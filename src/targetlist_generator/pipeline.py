"""One entry point for the CLI and GUI: cut a target list, pair standards, write lists and finder charts."""
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import NamedTuple

import pandas as pd
from astropy.coordinates import EarthLocation
from astropy.time import Time

from . import finder, io, schedule, sky, standards

Log = Callable[[str], None]


class Night(NamedTuple):
    times: Time
    location: EarthLocation
    max_airmass: float


def cut_targets(targets: pd.DataFrame, log: Log, date: str | None = None, site: str | None = None,
                lat: float | None = None, lon: float | None = None, twilight: float = -12.0,
                max_airmass: float = 2.0, min_time: float = 30.0, min_moon_sep: float = 0.0,
                ra_range: tuple | None = None, dec_range: tuple | None = None) -> tuple[pd.DataFrame, Night | None]:
    """Apply the RA/Dec box and, when ``date`` is given, the night-visibility and Moon cuts.
    Returns the kept targets (with a moon_sep column when a night is given) and the night."""
    if date is None and (site or lat is not None or lon is not None or min_moon_sep):
        raise ValueError('A site, latitude/longitude or Moon separation was given without a --date for the night')
    ra, dec = targets['ra'].to_numpy(), targets['dec'].to_numpy()
    keep = sky.in_box(ra, dec, ra_range, dec_range)
    if ra_range or dec_range:
        log(f'RA/Dec box keeps {keep.sum()} of {len(targets)} targets')
    night = None
    if date is not None:
        location = sky.site_location(site, lat, lon)
        times = sky.night_times(date, location, twilight)
        night = Night(times, location, max_airmass)
        log(f'Night of {date}: Sun below {twilight:g} deg from {times[0].iso[:16]} to {times[-1].iso[:16]} UTC')
        seen = sky.visible(ra, dec, times, location, max_airmass, min_time)
        log(f'{seen.sum()} of {len(targets)} targets reach airmass <= {max_airmass:g} for >= {min_time:g} min')
        body, lit = schedule.moon(times, location)
        targets = targets.assign(moon_sep=schedule.moon_separation(targets, body).round(1))
        far = targets['moon_sep'].to_numpy() >= min_moon_sep
        log(f'Moon {100 * lit:.0f} per cent illuminated'
            + (f'; {far.sum()} targets at least {min_moon_sep:g} deg from it' if min_moon_sep else ''))
        keep &= seen & far
    if not (ra_range or dec_range or date):
        log('No RA/Dec box or night given: keeping every target')
    if not keep.any():
        raise ValueError('No targets survive the cuts; check the date, site and limits')
    return targets[keep].reset_index(drop=True), night


def make_targetlist(target_list: str | Path, out_dir: str | Path = 'targetlist_out', mag_col: str = 'tmassk',
                    with_standards: bool = True, standards_file: str | Path | None = None, std_band: str = 'K',
                    std_spt: str | None = None, std_offset: float = 30.0, std_airmass_tol: float = 0.1,
                    order: str = 'time', finders: bool = False, survey: str = finder.SURVEYS[0], fov: float = 2.0,
                    urls: tuple[str, ...] = finder.HIPS2FITS, log: Log = print,
                    **cuts) -> tuple[pd.DataFrame, list[Path]]:
    """Write into ``out_dir``: ``<stem>_irtf.txt`` (targets interleaved with standards), ``<stem>_irtf_nostds.txt``,
    ``cats/<stem>_cat_N.cat``, with a night also ``<stem>_schedule.csv`` and ``<stem>_airmass.png``, and optionally
    ``fcharts/*.png``. ``cuts`` go to ``cut_targets``; ``order`` 'time' sorts by best time when a night is given,
    'input' keeps the input order. Returns the kept targets (plus any timing columns) and the written paths."""
    if order not in ('time', 'input'):
        raise ValueError(f"order must be 'time' or 'input', not {order!r}")
    out_dir = Path(out_dir)
    mag_col = mag_col.lower()
    targets = io.read_targets(target_list, mag_col)
    log(f'Read {len(targets)} targets from {target_list}')
    kept, night = cut_targets(targets, log, **cuts)
    if night is not None:
        kept = pd.concat([kept, schedule.timing(kept, *night)], axis=1)
        if order == 'time':
            kept = kept.sort_values('best_utc', kind='stable').reset_index(drop=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'cats').mkdir(exist_ok=True)
    stem = Path(target_list).stem
    paths = [io.write_irtf(out_dir / f'{stem}_irtf_nostds.txt', kept, mag_col)]
    if with_standards:
        stds = standards.load_standards(mag_col, std_band, std_spt, standards_file)
        if night is None:
            paired = standards.nearest_standards(kept, stds, std_offset)
            log(f'Paired each target with the nearest of {len(stds)} standards, {std_offset:g} min later in RA')
        else:
            paired, unmatched = standards.matched_standards(kept, stds, night.times, night.location,
                                                            std_offset, std_airmass_tol)
            log(f'Paired each target with the nearest of {len(stds)} standards within {std_airmass_tol:g} airmass '
                f'{std_offset:g} min after its best time' + (f'; {unmatched} had none, so got the closest '
                                                              'in airmass' if unmatched else ''))
        kept = kept.assign(standard=paired['shortname'].to_numpy())
        paths.append(io.write_irtf(out_dir / f'{stem}_irtf.txt', kept, mag_col, paired))
    paths += io.write_cats(out_dir / 'cats' / stem, kept, mag_col)
    if night is not None:
        columns = [c for c in ('shortname', 'best_utc', 'min_airmass', 'start_utc', 'end_utc', 'moon_sep',
                               'standard') if c in kept.columns]
        kept[columns].to_csv(out_dir / f'{stem}_schedule.csv', index=False)
        paths += [out_dir / f'{stem}_schedule.csv',
                  schedule.airmass_plot(kept, *night, out_dir / f'{stem}_airmass.png')]
    if finders:
        paths += make_finders(kept, kept['shortname'], out_dir, survey, fov, urls, log)
    log(f'Wrote {len(paths)} files to {out_dir}')
    return kept, paths


def make_finders(targets: pd.DataFrame, names: Iterable[str], out_dir: str | Path,
                 survey: str = finder.SURVEYS[0], fov: float = 2.0, urls: tuple[str, ...] = finder.HIPS2FITS,
                 log: Log = print) -> list[Path]:
    """Finder charts for ``names`` (matched to shortname, spaces as underscores) into ``out_dir/fcharts``;
    raises listing every chart that could not be made."""
    names = [str(n).strip().replace(' ', '_') for n in names]
    unknown = sorted(set(names) - set(targets['shortname']))
    if unknown:
        raise ValueError(f'Not in the target list: {unknown}')
    rows = targets.set_index('shortname')
    items = [(n, rows.at[n, 'ra'], rows.at[n, 'dec']) for n in names]
    paths, failed = finder.make_finders(items, Path(out_dir) / 'fcharts', survey, fov, urls, log)
    if failed:
        raise RuntimeError(f'{len(failed)} finder charts failed: {failed}')
    return paths
