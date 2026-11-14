"""One entry point for the CLI and GUI: cut a target list, plan the night, pair standards, write lists and charts."""
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
from astropy.coordinates import EarthLocation
from astropy.time import Time

from . import finder, formats, io, planner, schedule, sky, standards

Log = Callable[[str], None]
PLAN_COLUMNS = ('shortname', 'status', 'target_start', 'target_end', 'mean_airmass', 'standard', 'std_start',
                'std_airmass', 'std_sep', 'std_note', 'moon_sep', 'ra', 'dec')


class Night(NamedTuple):
    date: str
    title: str
    times: Time
    location: EarthLocation
    max_airmass: float


def cut_targets(targets: pd.DataFrame, log: Log, date: str | None = None, site: str | None = None,
                lat: float | None = None, lon: float | None = None, twilight: float = -12.0,
                max_airmass: float = 2.0, min_time: float = 30.0, min_moon_sep: str | float = 'auto',
                ra_range: tuple | None = None, dec_range: tuple | None = None) -> tuple[pd.DataFrame, Night | None]:
    """Apply the RA/Dec box and, when ``date`` is given, the night-visibility and Moon cuts
    (``min_moon_sep`` 'auto', 'off' or degrees; see ``schedule.moon_limit``).
    Returns the kept targets (with a moon_sep column when a night is given) and the night."""
    if date is None and (site or lat is not None or lon is not None):
        raise ValueError('A site or latitude/longitude was given without a --date for the night')
    ra, dec = targets['ra'].to_numpy(), targets['dec'].to_numpy()
    keep = sky.in_box(ra, dec, ra_range, dec_range)
    if ra_range or dec_range:
        log(f'RA/Dec box keeps {keep.sum()} of {len(targets)} targets')
    night = None
    if date is not None:
        location = sky.site_location(site, lat, lon)
        times = sky.night_times(date, location, twilight)
        night = Night(date, site or f'{lat:.4f}, {lon:.4f}', times, location, max_airmass)
        log(f'Night of {date}: Sun below {twilight:g} deg from {times[0].iso[:16]} to {times[-1].iso[:16]} UTC')
        seen = sky.visible(ra, dec, times, location, max_airmass, min_time)
        log(f'{seen.sum()} of {len(targets)} targets reach airmass <= {max_airmass:g} for >= {min_time:g} min')
        body, lit = schedule.moon(times, location)
        limit = schedule.moon_limit(min_moon_sep, lit)
        targets = targets.assign(moon_sep=schedule.moon_separation(targets, body).round(1))
        far = targets['moon_sep'].to_numpy() >= limit
        log(f'Moon {100 * lit:.0f} per cent illuminated'
            + (f'; {far.sum()} targets at least {limit:.0f} deg from it' if limit else ''))
        keep &= seen & far
    if not (ra_range or dec_range or date):
        log('No RA/Dec box or night given: keeping every target')
    if not keep.any():
        raise ValueError('No targets survive the cuts; check the date, site and limits')
    return targets[keep].reset_index(drop=True), night


def target_minutes(targets: pd.DataFrame, target_time: float, target_time_col: str | None) -> np.ndarray:
    """Observing time per target [min]: the ``target_time_col`` column when given, else ``target_time``."""
    if target_time_col is None:
        return np.full(len(targets), float(target_time))
    minutes = pd.to_numeric(targets[target_time_col.lower()], errors='coerce').to_numpy(dtype=float)
    bad = targets['shortname'][~(np.isfinite(minutes) & (minutes > 0))].tolist()
    if bad:
        raise ValueError(f'Column {target_time_col} needs a positive number of minutes for {bad}')
    return minutes


def make_targetlist(target_list: str | Path, out_dir: str | Path, mag_col: str = 'tmassk',
                    target_time: float = 30.0, target_time_col: str | None = None, cal_time: float = 5.0,
                    with_standards: bool = True, standards_file: str | Path | None = None, std_band: str = 'K',
                    std_spt: str | None = None, std_time: float = 10.0, std_airmass_tol: float = 0.1,
                    std_max_sep: float = 20.0, outputs: Iterable[str] | None = None, finders: bool = False,
                    survey: str = finder.SURVEYS[0], fov: float = 2.0, urls: tuple[str, ...] = finder.HIPS2FITS,
                    log: Log = print, show=None, progress=None, **cuts) -> tuple[pd.DataFrame, list[Path]]:
    """Cut ``target_list`` (``cuts`` go to ``cut_targets``), plan the night when a date is given, pair standards
    and write into ``out_dir``: the chosen ``outputs`` (see ``formats.FORMATS``; None = ``formats.defaults`` for the
    site), with a night also
    ``<stem>_plan.csv`` and ``<stem>_visibility_N.png``, and optionally ``fcharts/*.png``.
    ``show(kind, name, item)`` receives each visibility figure (kind "visibility") and finder chart (kind
    "finder": a Figure, or a PNG path when drawn by Aladin) as it is made; ``progress(message, done, total)``
    reports the current stage (total 0 while busy). Returns the plan (or the kept targets without a night)
    and the written paths."""
    progress = progress or (lambda message, done, total: None)
    progress('Planning the night', 0, 0)
    out_dir = Path(out_dir)
    mag_col = mag_col.lower()
    targets = io.read_targets(target_list, mag_col, keep=(target_time_col,) if target_time_col else ())
    log(f'Read {len(targets)} targets from {target_list}')
    kept, night = cut_targets(targets, log, **cuts)
    stds = standards.load_standards(mag_col, std_band, std_spt, standards_file) if with_standards else None
    if night is None:
        rows, paired = kept, None
        if stds is not None:
            paired = standards.nearest_standards(kept, stds, target_time + cal_time)
            log(f'No night given: paired each target with the nearest of {len(stds)} standards on the sky')
    else:
        minutes = target_minutes(kept, target_time, target_time_col)
        rows, paired = planner.plan_night(kept, night.times, night.location, night.max_airmass, minutes, cal_time,
                                          std_time if stds is not None else 0.0, stds, std_airmass_tol, std_max_sep)
        extra = int((rows['status'] == 'extra').sum())
        log(f'Planned {len(rows) - extra} blocks back to back' + (f'; {extra} more did not fit and are '
                                                                  'listed at their best time' if extra else ''))
        if paired is not None and (rows['std_note'] != '').any():
            log(f'{(rows["std_note"] != "").sum()} standards are compromises: see std_note in the plan')
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / Path(target_list).stem
    keys = list(outputs) if outputs is not None else list(formats.defaults(cuts.get('site')))
    paths = formats.write(keys, stem, rows, paired, mag_col)
    if night is not None:
        plan_path = stem.with_name(f'{stem.name}_plan.csv')
        rows[[c for c in PLAN_COLUMNS if c in rows.columns]].round({'ra': 6, 'dec': 6}).to_csv(plan_path, index=False)
        progress('Drawing visibility plots', 0, 0)
        paths += [plan_path, *schedule.visibility_plots(
            rows, night.date, night.location, night.max_airmass, night.title, stem,
            None if show is None else lambda name, fig: show('visibility', name, fig))]
    if finders:
        paths += make_finders(rows, rows['shortname'], out_dir, survey, fov, urls, log,
                              None if show is None else lambda name, item: show('finder', name, item), progress)
    log(f'Wrote {len(paths)} files to {out_dir}')
    return rows, paths


def make_finders(targets: pd.DataFrame, names: Iterable[str], out_dir: str | Path,
                 survey: str = finder.SURVEYS[0], fov: float = 2.0, urls: tuple[str, ...] = finder.HIPS2FITS,
                 log: Log = print, show=None, progress=None) -> list[Path]:
    """Finder charts for ``names`` (matched to shortname, spaces as underscores) into ``out_dir/fcharts``;
    raises listing every chart that could not be made. ``show``/``progress`` as in ``finder.make_finders``."""
    names = [str(n).strip().replace(' ', '_') for n in names]
    unknown = sorted(set(names) - set(targets['shortname']))
    if unknown:
        raise ValueError(f'Not in the target list: {unknown}')
    rows = targets.set_index('shortname')
    items = [(n, rows.at[n, 'ra'], rows.at[n, 'dec']) for n in names]
    paths, failed = finder.make_finders(items, Path(out_dir) / 'fcharts', survey, fov, urls, log, show, progress)
    if failed:
        raise RuntimeError(f'{len(failed)} finder charts failed: {failed}')
    return paths
