"""One entry point for the CLI and GUI: cut a target list, pair standards, write lists and finder charts."""
from collections.abc import Callable, Iterable
from pathlib import Path

import pandas as pd

from . import finder, io, sky, standards

Log = Callable[[str], None]


def cut_targets(targets: pd.DataFrame, log: Log, date: str | None = None, site: str | None = None,
                lat: float | None = None, lon: float | None = None, twilight: float = -12.0,
                max_airmass: float = 2.0, min_time: float = 30.0,
                ra_range: tuple | None = None, dec_range: tuple | None = None) -> pd.DataFrame:
    """Apply the RA/Dec box and, when ``date`` is given, the night-visibility cut."""
    if date is None and (site or lat is not None or lon is not None):
        raise ValueError('A site or latitude/longitude was given without a --date for the night')
    keep = sky.in_box(targets['ra'].to_numpy(), targets['dec'].to_numpy(), ra_range, dec_range)
    if ra_range or dec_range:
        log(f'RA/Dec box keeps {keep.sum()} of {len(targets)} targets')
    if date is not None:
        location = sky.site_location(site, lat, lon)
        times = sky.night_times(date, location, twilight)
        log(f'Night of {date}: Sun below {twilight:g} deg from {times[0].iso[:16]} to {times[-1].iso[:16]} UTC')
        seen = sky.visible(targets['ra'].to_numpy(), targets['dec'].to_numpy(), times, location,
                           max_airmass, min_time)
        log(f'{seen.sum()} of {len(targets)} targets reach airmass <= {max_airmass:g} for >= {min_time:g} min')
        keep &= seen
    if not (ra_range or dec_range or date):
        log('No RA/Dec box or night given: keeping every target')
    if not keep.any():
        raise ValueError('No targets survive the cuts; check the date, site and limits')
    return targets[keep].reset_index(drop=True)


def make_targetlist(target_list: str | Path, out_dir: str | Path = 'targetlist_out', mag_col: str = 'tmassk',
                    with_standards: bool = True, standards_file: str | Path | None = None, std_band: str = 'K',
                    std_spt: str | None = None, std_offset: float = 30.0, finders: bool = False,
                    survey: str = finder.SURVEYS[0], fov: float = 2.0, urls: tuple[str, ...] = finder.HIPS2FITS,
                    log: Log = print, **cuts) -> tuple[pd.DataFrame, list[Path]]:
    """Write ``<stem>_irtf.txt`` (targets interleaved with standards), ``<stem>_irtf_nostds.txt``,
    ``cats/<stem>_cat_N.cat`` and optionally ``fcharts/*.png`` into ``out_dir``; ``cuts`` go to ``cut_targets``.
    Returns the kept targets and the written paths."""
    out_dir = Path(out_dir)
    mag_col = mag_col.lower()
    targets = io.read_targets(target_list, mag_col)
    log(f'Read {len(targets)} targets from {target_list}')
    kept = cut_targets(targets, log, **cuts)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'cats').mkdir(exist_ok=True)
    stem = Path(target_list).stem
    paths = [io.write_irtf(out_dir / f'{stem}_irtf_nostds.txt', kept, mag_col)]
    if with_standards:
        stds = standards.load_standards(mag_col, std_band, std_spt, standards_file)
        paired = standards.nearest_standards(kept, stds, std_offset)
        paths.append(io.write_irtf(out_dir / f'{stem}_irtf.txt', kept, mag_col, paired))
        log(f'Paired each target with the nearest of {len(stds)} standards, {std_offset:g} min later in RA')
    paths += io.write_cats(out_dir / 'cats' / stem, kept, mag_col)
    if finders:
        paths += make_finders(kept, kept['shortname'], out_dir, survey, fov, urls, log)
    log(f'Wrote {len(paths)} files to {out_dir}')
    return kept, paths


def make_finders(targets: pd.DataFrame, names: Iterable[str], out_dir: str | Path,
                 survey: str = finder.SURVEYS[0], fov: float = 2.0, urls: tuple[str, ...] = finder.HIPS2FITS,
                 log: Log = print) -> list[Path]:
    """Finder charts for ``names`` (matched to shortname, spaces as underscores) into ``out_dir/fcharts``.
    A failed download is logged and skipped; an error listing every failure is raised at the end."""
    names = [str(n).strip().replace(' ', '_') for n in names]
    unknown = sorted(set(names) - set(targets['shortname']))
    if unknown:
        raise ValueError(f'Not in the target list: {unknown}')
    rows = targets.set_index('shortname')
    paths, failed = [], []
    for name in names:
        try:
            paths.append(finder.make_finder(name, rows.at[name, 'ra'], rows.at[name, 'dec'],
                                            Path(out_dir) / 'fcharts', survey, fov, urls))
            log(f'Finder chart: {paths[-1].name}')
        except OSError as err:  # network/HTTP errors; one bad target should not cost the rest
            failed.append(name)
            log(f'Finder chart FAILED for {name}: {err}')
    if failed:
        raise RuntimeError(f'{len(failed)} finder charts failed: {failed}')
    return paths
