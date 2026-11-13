import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.coordinates import AltAz, SkyCoord

from targetlist_generator import cli, io, make_targetlist, sky, standards


def test_spt_from_num():
    assert [io.spt_from_num(n) for n in (60, 72.5, 75.2, 89.9, 90)] == ['M0', 'L2.5', 'L5', 'Y0', 'Y0']
    for bad in (59, 99.8, 105):
        with pytest.raises(ValueError):
            io.spt_from_num(bad)


def test_read_targets_mixed_formats(target_csv):
    df = io.read_targets(target_csv)
    assert list(df.columns) == ['shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', 'tmassk']
    assert df.loc[0, 'dec'] == pytest.approx(-0.5)
    assert ' ' not in ''.join(df['shortname'])
    assert set(df['spt']) <= {'M5', 'L2.5', 'T6', 'L1'}
    assert df['pmra'].eq(0).all() and df['epoch'].eq(2000).all()


def test_read_targets_missing_columns(tmp_path):
    path = tmp_path / 'bad.csv'
    pd.DataFrame({'name': ['a'], 'ra': [1.0], 'dec': [1.0]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match=r"\['spt', 'tmassk'\]"):
        io.read_targets(path)


def test_formatted_sign_and_wrap():
    """Old code lost the sign of -0.5 deg and wrote 60 s / 24 h."""
    df = pd.DataFrame({'shortname': ['a', 'b'], 'ra': [359.9999999, 10.0], 'dec': [-0.5, 89.99999999],
                       'epoch': [2000.0] * 2, 'pmra': [0.0] * 2, 'pmdec': [0.0] * 2, 'spt': ['L1'] * 2, 'k': [1.0] * 2})
    out = io.formatted(df, 'k')
    assert list(out['ra']) == ['00:00:00.00', '00:40:00.00']
    assert list(out['dec']) == ['-00:30:00.0', '90:00:00.0']


def test_dec_box_and_ra_wrap():
    ra, dec = np.array([350.0, 10.0, 180.0]), np.array([-80.0, 0.0, 85.0])
    assert list(sky.in_box(ra, dec, ('23h', '1h'))) == [True, True, False]
    assert list(sky.in_box(ra, dec, dec_range=('-90', '10'))) == [True, True, False]


def test_airmass_matches_astropy():
    loc = sky.site_location(lat=28.76, lon=-17.88)
    times = sky.night_times('2026-11-14', loc, -18)
    ra, dec = np.array([30.0, 300.0]), np.array([20.0, -10.0])
    ours = sky.airmass(ra, dec, times[::12], loc)
    alt = SkyCoord(ra * u.deg, dec * u.deg)[None, :].transform_to(AltAz(obstime=times[::12, None], location=loc)).alt
    up = alt > 20 * u.deg
    assert np.allclose(ours[up], 1 / np.sin(alt[up]), rtol=0.002)  # no refraction in either


def test_night_southern_and_polar():
    assert len(sky.night_times('2026-06-21', sky.site_location(lat=-24.6, lon=-70.4))) > 100
    with pytest.raises(ValueError, match='never sets'):
        sky.night_times('2026-06-21', sky.site_location(lat=78.2, lon=15.6))


def test_nearest_standards_brute_force(target_csv):
    targets = io.read_targets(target_csv)
    stds = standards.load_standards('tmassk', spt='A0V')
    got = standards.nearest_standards(targets, stds, offset=30)
    assert len(got) == len(targets) and set(got['spt']) == {'A0V'}
    for i in range(len(targets)):
        c = SkyCoord((targets.ra[i] + 7.5) * u.deg, targets.dec[i] * u.deg)
        sep = c.separation(SkyCoord(stds.ra.to_numpy() * u.deg, stds.dec.to_numpy() * u.deg))
        assert got.shortname[i] == stds.shortname.iloc[np.argmin(sep)]


def test_make_targetlist_outputs(target_csv, tmp_path):
    """Old code wrote 16 misaligned columns and overlapping .cat chunks."""
    kept, paths = make_targetlist(target_csv, tmp_path / 'out', log=lambda s: None)
    irtf = (tmp_path / 'out' / 'targets_irtf.txt').read_text().splitlines()
    body = [line.split('\t') for line in irtf[3:]]
    assert len(body) == 2 * len(kept) == 50 and {len(row) for row in body} == {9}
    assert [row[0] for row in body] == [str(i) for i in range(1, 51)]
    assert all(row[7] in ('B9V', 'A0V', 'A1V') for row in body[1::2])
    cats = sorted((tmp_path / 'out' / 'cats').glob('*.cat'))
    names = [line.split()[0] for p in cats for line in p.read_text().splitlines()]
    assert len(cats) == 3 and names == list(kept['shortname'])


def test_night_cut_keeps_only_visible(target_csv, tmp_path):
    kept, _ = make_targetlist(target_csv, tmp_path, with_standards=False, log=lambda s: None,
                              date='2026-11-14', lat=19.83, lon=-155.47, max_airmass=1.5, min_time=60)
    assert 0 < len(kept) < 25
    loc = sky.site_location(lat=19.83, lon=-155.47)
    am = sky.airmass(kept.ra.to_numpy(), kept.dec.to_numpy(), sky.night_times('2026-11-14', loc), loc)
    assert ((am <= 1.5).sum(axis=0) * 5 >= 60).all()


def test_cli_errors_cleanly(target_csv, tmp_path):
    with pytest.raises(SystemExit, match='without a --date'):
        cli.main(['make', str(target_csv), '-o', str(tmp_path), '--lat', '19'])
    with pytest.raises(SystemExit, match='Not in the target list'):
        cli.main(['finder', str(target_csv), 'nobody', '-o', str(tmp_path)])


@pytest.mark.network
def test_finder_chart(target_csv, tmp_path):
    cli.main(['finder', str(target_csv), 'J0003 test', '-o', str(tmp_path), '--fov', '1'])
    assert (tmp_path / 'fcharts' / 'J0003_test_1arcmin.png').stat().st_size > 10_000


def test_plot_finder_offline(tmp_path):
    from astropy.io import fits
    from astropy.wcs import WCS

    from targetlist_generator import finder
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [150.0, 22.0]
    wcs.wcs.crpix = [finder.PIXELS / 2] * 2
    wcs.wcs.cdelt = [-2 / 60 / finder.PIXELS, 2 / 60 / finder.PIXELS]
    data = np.random.default_rng(0).normal(size=(finder.PIXELS, finder.PIXELS))
    hdu = fits.PrimaryHDU(data, header=wcs.to_header())
    path = finder.plot_finder(hdu, 'J1000+2200', 'CDS/P/2MASS/K', 2.0, tmp_path / 'f.png')
    assert path.stat().st_size > 10_000


IRTF = {'lat': 19.83, 'lon': -155.47}


def test_matched_standards_brute_force(target_csv):
    loc = sky.site_location(**IRTF)
    times = sky.night_times('2026-11-14', loc)
    targets = io.read_targets(target_csv).iloc[:6]
    stds = standards.load_standards('tmassk', spt='A0V')
    got, _ = standards.matched_standards(targets, stds, times, loc, offset=30, tol=0.1)
    t_am = sky.airmass(targets.ra.to_numpy(), targets.dec.to_numpy(), times, loc)
    s_am = sky.airmass(stds.ra.to_numpy(), stds.dec.to_numpy(), times, loc)
    for i in range(len(targets)):
        best = t_am[:, i].argmin()
        with np.errstate(invalid='ignore'):
            diff = np.nan_to_num(np.abs(s_am[min(best + 6, len(times) - 1)] - t_am[best, i]), nan=np.inf)
        sep = SkyCoord(targets.ra[i] * u.deg, targets.dec[i] * u.deg).separation(
            SkyCoord(stds.ra.to_numpy() * u.deg, stds.dec.to_numpy() * u.deg)).deg
        want = np.argmin(np.where(diff <= 0.1, sep, np.inf)) if (diff <= 0.1).any() else np.argmin(diff)
        assert got.shortname[i] == stds.shortname.iloc[want]


def test_schedule_order_and_files(target_csv, tmp_path):
    kept, paths = make_targetlist(target_csv, tmp_path, log=lambda s: None, date='2026-11-14', **IRTF)
    assert list(kept['best_utc']) == sorted(kept['best_utc'])
    sched = pd.read_csv(tmp_path / 'targets_schedule.csv')
    assert list(sched.columns) == ['shortname', 'best_utc', 'min_airmass', 'start_utc', 'end_utc', 'moon_sep',
                                   'standard']
    assert (sched['start_utc'] <= sched['best_utc']).all() and (sched['best_utc'] <= sched['end_utc']).all()
    assert (sched['min_airmass'] <= 2).all()
    irtf = (tmp_path / 'targets_irtf.txt').read_text().splitlines()[3::2]
    assert [line.split('\t')[1] for line in irtf] == list(kept['shortname'])
    assert (tmp_path / 'targets_airmass.png').stat().st_size > 20_000
    kept_in, _ = make_targetlist(target_csv, tmp_path, log=lambda s: None, date='2026-11-14', order='input', **IRTF)
    names = list(io.read_targets(target_csv)['shortname'])
    assert sorted(kept_in['shortname'], key=names.index) == list(kept_in['shortname'])


def test_moon_separation_and_cut(target_csv, tmp_path):
    from astropy.coordinates import get_body

    from targetlist_generator import schedule
    loc = sky.site_location(**IRTF)
    times = sky.night_times('2026-11-14', loc)
    targets = io.read_targets(target_csv)
    body, lit = schedule.moon(times, loc)
    assert 0 <= lit <= 1
    ref = get_body('moon', times[len(times) // 2], loc)
    stars = SkyCoord(targets.ra.to_numpy() * u.deg, targets.dec.to_numpy() * u.deg).transform_to(ref.frame)
    assert np.allclose(schedule.moon_separation(targets, body), ref.separation(stars).deg, atol=1e-6)
    kept, _ = make_targetlist(target_csv, tmp_path, log=lambda s: None, date='2026-11-14', min_moon_sep=60, **IRTF)
    assert (kept['moon_sep'] >= 60).all()


@pytest.mark.network
def test_aladin_batch(tmp_path):
    import shutil

    from targetlist_generator import finder
    if shutil.which('aladin') is None:
        pytest.skip('aladin not installed')
    items = [('A', 150.0, 22.0), ('B', 151.0, 23.0)]
    paths, failed = finder.make_finders(items, tmp_path, urls=('http://127.0.0.1:9/none',), log=lambda s: None)
    assert failed == [] and [p.name for p in paths] == ['A_2arcmin.png', 'B_2arcmin.png']
