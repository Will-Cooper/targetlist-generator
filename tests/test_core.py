import re

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.coordinates import AltAz, SkyCoord
from astropy.table import Table

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
    kept, paths = make_targetlist(target_csv, tmp_path / 'out', outputs=['irtf', 'cat'], log=lambda s: None)
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
                              date='2026-11-14', lat=19.83, lon=-155.47, max_airmass=1.5, min_time=60,
                              min_moon_sep='off')
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
    path = tmp_path / 'f.png'
    finder.finder_figure(hdu, 'J1000+2200', 'CDS/P/2MASS/K', 2.0).savefig(path)
    assert path.stat().st_size > 10_000


IRTF = {'lat': 19.83, 'lon': -155.47}


def test_site_lookup():
    assert sky.site_location('IRTF').lat.deg == pytest.approx(19.826, abs=1e-3)
    assert sky.site_location('lapalma').lon.deg == pytest.approx(-17.88, abs=1e-2)
    with pytest.raises(ValueError, match="did you mean.*paranal"):
        sky.site_location('Paranel')


def test_moon_limit():
    from targetlist_generator import schedule
    assert schedule.moon_limit('auto', 0) == 10 and schedule.moon_limit('auto', 1) == 40
    assert schedule.moon_limit('off', 0.5) == 0 and schedule.moon_limit('25', 0.5) == 25
    with pytest.raises(ValueError):
        schedule.moon_limit('200', 0.5)


def test_standard_never_across_sky():
    """A standard at the target's airmass 180 deg away must lose to a nearby one slightly off in airmass."""
    loc = sky.site_location(lat=0.0, lon=0.0)
    times = sky.night_times('2026-03-20', loc)
    k = len(times) // 2
    lst = times[k].sidereal_time('apparent', longitude=loc.lon).deg
    target = pd.DataFrame({'shortname': ['T'], 'ra': [lst], 'dec': [10.0]})  # transiting, altitude 80 deg
    stds = pd.DataFrame({'shortname': ['far', 'near'], 'ra': [(lst + 180) % 360, lst + 5], 'dec': [-10.0, 25.0]})
    am_t = sky.airmass(target.ra.to_numpy(), target.dec.to_numpy(), times[k:k + 1], loc)[0]
    got, notes, _, sep = standards.pick(target, stds, times, loc, np.array([k]), am_t, tol=0.01, max_sep=20)
    assert got.shortname[0] == 'near' and sep[0] < 20 and notes[0] == 'no standard within 0.01 airmass'


def test_planner_no_overlap_no_drop(target_csv):
    from targetlist_generator import planner
    loc = sky.site_location(**IRTF)
    times = sky.night_times('2026-11-30', loc)
    targets = io.read_targets(target_csv)
    targets = targets[sky.visible(targets.ra.to_numpy(), targets.dec.to_numpy(), times, loc)].reset_index(drop=True)
    minutes = np.full(len(targets), 45.0)
    plan, paired = planner.plan_night(targets, times, loc, 2.0, minutes, 5, 10, standards.load_standards('tmassk'))
    assert sorted(plan.shortname) == sorted(targets.shortname)  # nothing dropped
    assert list(plan.target_start) == sorted(plan.target_start) and len(paired) == len(plan)
    blocks = plan[plan.status == 'planned']
    ends = pd.to_datetime(blocks.target_start) + pd.Timedelta(minutes=60)  # 45 + 5 + 10
    assert (pd.to_datetime(blocks.target_start).iloc[1:].to_numpy() >= ends.iloc[:-1].to_numpy()).all()
    assert (blocks.target_start >= times[0].iso[:16]).all() and (ends <= pd.Timestamp(times[-1].iso) +
                                                                 pd.Timedelta(minutes=5)).all()
    assert (plan.mean_airmass[plan.status == 'planned'] <= 2).all()


def test_target_time_column_and_plan_files(target_csv, tmp_path):
    df = pd.read_csv(target_csv).assign(exptime=20)
    df.loc[0, 'exptime'] = -1
    path = tmp_path / 'timed.csv'
    df.to_csv(path, index=False)
    with pytest.raises(ValueError, match='positive number of minutes'):
        make_targetlist(path, tmp_path, target_time_col='exptime', min_moon_sep='off', log=lambda s: None,
                        date='2026-11-30', ra_range=('0', '360'), **IRTF)
    df.loc[0, 'exptime'] = 20
    df.to_csv(path, index=False)
    plan, paths = make_targetlist(path, tmp_path, target_time_col='exptime', log=lambda s: None,
                                  date='2026-11-30', **IRTF)
    starts, ends = pd.to_datetime(plan.target_start), pd.to_datetime(plan.target_end)
    assert ((ends - starts) == pd.Timedelta(minutes=20)).all()
    assert (tmp_path / 'timed_plan.csv').exists() and (tmp_path / 'timed_visibility_1.png').stat().st_size > 30_000


def test_every_format(target_csv, tmp_path):
    from targetlist_generator import formats
    plan, paths = make_targetlist(target_csv, tmp_path, outputs=list(formats.FORMATS), log=lambda s: None)
    names = {p.name for p in paths}
    for expected in ('targets_irtf.txt', 'targets_irtf_nostds.txt', 'targets_list.txt', 'targets_list.csv',
                     'targets_list.vot', 'targets_keck.txt', 'targets_lick.txt', 'targets_eso_p2.txt',
                     'targets_subaru_ope.txt', 'targets_mmt.cat', 'targets_magellan.cat', 'targets_cat_1.cat'):
        assert expected in names
    n = 2 * len(plan)  # each target followed by its standard
    keck = (tmp_path / 'targets_keck.txt').read_text().splitlines()
    assert len(keck) == n and all(len(line.split()) >= 8 and line[15] == ' ' for line in keck)
    assert keck[0][16:27] == plan.pipe(io.sexagesimal, sep=' ', sign=True)[0][0]
    assert all(len(line.split(', ')) == 5 for line in (tmp_path / 'targets_eso_p2.txt').read_text().splitlines())
    subaru = (tmp_path / 'targets_subaru_ope.txt').read_text().splitlines()
    assert all(re.fullmatch(r'\w+=OBJECT="[^"]+" RA=\d{6}\.\d{3} DEC=[+-]\d{6}\.\d{2} EQUINOX=2000\.0', s)
               for s in subaru)
    assert all(len(line.split()) == 8 for line in (tmp_path / 'targets_mmt.cat').read_text().splitlines())
    magellan = (tmp_path / 'targets_magellan.cat').read_text().splitlines()
    assert magellan[0].startswith('#') and all(len(line.split()) == 16 for line in magellan[1:])
    assert len(Table.read(tmp_path / 'targets_list.vot')) == n
    with pytest.raises(SystemExit):
        cli.main(['make', str(target_csv), '-o', str(tmp_path), '--formats', 'irtf,nonsense'])


def test_keck_proper_motion_units():
    from targetlist_generator import formats
    rows = pd.DataFrame({'dec': [60.0], 'pmra': [150.0], 'pmdec': [-300.0]})  # mas/yr, pmra with cos dec
    pmra, pmdec = formats.time_pm(rows)
    assert pmra[0] == pytest.approx(0.15 / 15 / 0.5) and pmdec[0] == pytest.approx(-0.3)


@pytest.mark.network
def test_aladin_batch(tmp_path):
    import shutil

    from targetlist_generator import finder
    if shutil.which('aladin') is None:
        pytest.skip('aladin not installed')
    items = [('A', 150.0, 22.0), ('B', 151.0, 23.0)]
    paths, failed = finder.make_finders(items, tmp_path, urls=('http://127.0.0.1:9/none',), log=lambda s: None)
    assert failed == [] and [p.name for p in paths] == ['A_2arcmin.png', 'B_2arcmin.png']


def spec_row(name, ra, dec, pmra=0.0, pmdec=0.0, epoch=2000.0, mag=0.0, spt=''):
    return pd.DataFrame({'shortname': [name], 'ra': [sky.parse_angle(ra, u.hourangle)],
                         'dec': [sky.parse_angle(dec, u.deg)], 'epoch': [epoch], 'pmra': [pmra], 'pmdec': [pmdec],
                         'spt': [spt], 'k': [mag]})


def test_formats_match_spec_examples(tmp_path):
    """Lines rebuilt from the examples in each observatory's format documentation."""
    from targetlist_generator import formats
    formats.subaru(tmp_path / 's', spec_row('M51', '13:29:55', '+47:18:55'), None, 'k')
    assert (tmp_path / 's_subaru_ope.txt').read_text() == \
        'M51=OBJECT="M51" RA=132955.000 DEC=+471855.00 EQUINOX=2000.0\n'
    formats.eso(tmp_path / 'e', spec_row('Sz100', '16:08:25.8', '-39:06:01', epoch=2013.5), None, 'k')
    assert (tmp_path / 'e_eso_p2.txt').read_text() == 'Sz100, 16:08:25.800, -39:06:01.00, J2000, 2013.5\n'
    # MMT example 'FK5-0001 0:08:23.265 +29:05:25.58 +01.039 -016.33 2.06 A0 J2000.0' (pm per century)
    cos_dec = np.cos(np.radians(sky.parse_angle('+29:05:25.58', u.deg)))
    row = spec_row('FK5-0001', '0:08:23.265', '+29:05:25.58', pmra=1.039 / 100 * 15 * 1000 * cos_dec,
                   pmdec=-16.33 / 100 * 1000, mag=2.06, spt='A0')
    # Keck LGS-page example 'NIR standard 21 30 31.700 +12 37 51.30 2000.0': name in columns 1-16
    formats.keck(tmp_path / 'k', spec_row('NIR_standard', '21:30:31.70', '+12:37:51.3'), None, 'k')
    assert (tmp_path / 'k_keck.txt').read_text() == 'NIR_standard    21 30 31.70 +12 37 51.3 2000.0\n'
    # Lick example 'obj1a  12 34 56  1 2 3  2000.0' (free-form whitespace, so padding and precision may differ)
    formats.lick(tmp_path / 'l', spec_row('obj1a', '12:34:56', '1:02:03'), None, 'k')
    assert (tmp_path / 'l_lick.txt').read_text() == 'obj1a 12 34 56.00 +01 02 03.0 2000.0\n'
    formats.mmt(tmp_path / 'm', row, None, 'k')
    assert (tmp_path / 'm_mmt.cat').read_text() == 'FK5-0001 00:08:23.265 +29:05:25.58 +1.039 -16.33 2.06 A0 J2000.0\n'


def test_default_formats_follow_observatory(target_csv, tmp_path):
    from targetlist_generator import formats
    assert formats.defaults('irtf') == ('irtf',) and formats.defaults('keck') == ('keck',)
    assert formats.defaults('paranal') == ('eso',) and formats.defaults('lapalma') == ('cat',)
    assert formats.defaults('lco') == ('magellan',) and formats.defaults('Subaru') == ('subaru',)
    assert formats.defaults('ALMA') == ('csv',) and formats.defaults(None) == ('csv',)
    assert formats.defaults(None, from_samp=True) == ('votable',)
    assert formats.defaults('mmt', from_samp=True) == ('mmt',)  # the observatory wins
    cli.main(['make', str(target_csv), '-o', str(tmp_path), '--date', '2026-11-30', '--site', 'keck'])
    assert (tmp_path / 'targets_keck.txt').exists() and not (tmp_path / 'targets_list.csv').exists()
