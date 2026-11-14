"""Target-list writers for telescopes and tools. Each takes the targets in observing order, the paired standards
(one row per target, or None) and the magnitude column, and returns the written paths.

Input proper motions are mas/yr with pmra including cos(dec), as in Gaia. Keck, MMT and Magellan convert them to
their documented units; Lick and ESO omit them (cos dec convention unstated); IRTF copies them unchanged (units not
confirmed, see ``irtf``). Coordinates are written at the catalogue epoch with equinox J2000."""
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

from .io import formatted, sexagesimal
from .sky import find_site

CAT_CHUNK = 10  # targets per ING .cat file
RULE = '#+-----+-----+-----+-----+-----+-----+-----+-----+-------+'
Writer = Callable[[Path, pd.DataFrame, pd.DataFrame | None, str], list[Path]]


def interleave(targets: pd.DataFrame, stds: pd.DataFrame | None) -> pd.DataFrame:
    """Targets each followed by their standard, with a ``kind`` column (target/standard)."""
    rows = targets.assign(kind='target').reset_index(drop=True)
    if stds is None:
        return rows
    std_rows = stds.assign(kind='standard').reset_index(drop=True)
    return pd.concat([rows, std_rows]).sort_index(kind='stable').reset_index(drop=True)  # t0, s0, t1, s1, ...


def time_pm(rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Proper motion as seconds of RA per year (no cos dec) and arcsec per year."""
    cos_dec = np.cos(np.radians(rows['dec'].to_numpy()))
    return rows['pmra'].to_numpy() / 1000 / 15 / cos_dec, rows['pmdec'].to_numpy() / 1000


def clean_names(names: pd.Series, disallowed: str | None, fmt: str, max_len: int | None = None) -> pd.Series:
    """Names with characters matching ``disallowed`` removed (spaces collapsed); raises when a name ends up empty
    or longer than ``max_len``, or two different objects end up with the same name."""
    out = names if disallowed is None else names.str.replace(disallowed, '', regex=True)
    out = out.str.replace(r'\s+', ' ', regex=True).str.strip()
    too_long = out.str.len() > max_len if max_len else pd.Series(False, index=out.index)
    bad = sorted(set(names[(out == '') | too_long]))
    if bad:
        raise ValueError(f'{fmt} names must be 1-{max_len or "any"} characters after cleaning; rename {bad}')
    pairs = pd.DataFrame({'old': names, 'new': out}).drop_duplicates()
    clash = sorted(pairs.loc[pairs['new'].duplicated(keep=False), 'old'])
    if clash:
        raise ValueError(f'{fmt} names clash after removing disallowed characters; rename {clash}')
    return out


def irtf(stem: Path, targets, stds, mag_col) -> list[Path]:
    """IRTF starcat lists: ``_irtf.txt`` with standards (when given) and ``_irtf_nostds.txt``, in the version 1
    layout. pmra/pmdec are copied from the table unchanged: the units the StarCat expects (IRTF TCS3 StarCat user
    guide, 'Userlist v3 text format') could not be confirmed, so check them before relying on proper motions."""
    def write(path: Path, rows: pd.DataFrame) -> Path:
        out = formatted(rows, mag_col)
        out.index = np.arange(1, len(out) + 1)
        header = f'#|Index| Name|  RA | DEC |Epoch| PMRA|PMDEC|SpT  | {mag_col.upper():>6}|'
        path.write_text('\n'.join((RULE, header, RULE, '')) + out.to_csv(sep='\t', header=False, lineterminator='\n'))
        return path
    paths = [write(stem.with_name(f'{stem.name}_irtf_nostds.txt'), targets)]
    if stds is not None:
        paths.append(write(stem.with_name(f'{stem.name}_irtf.txt'), interleave(targets, stds)))
    return paths


def ing_cat(stem: Path, targets, stds, mag_col) -> list[Path]:
    """ING object-visibility ``cats/<stem>_cat_N.cat`` files, ``CAT_CHUNK`` targets each: name ra dec epoch."""
    folder = stem.parent / 'cats'
    folder.mkdir(exist_ok=True)
    rows = formatted(targets, mag_col)[['shortname', 'ra', 'dec', 'epoch']]
    paths = []
    for k, start in enumerate(range(0, len(rows), CAT_CHUNK), 1):
        path = folder / f'{stem.name}_cat_{k}.cat'
        rows.iloc[start:start + CAT_CHUNK].to_csv(path, sep=' ', header=False, index=False, lineterminator='\n')
        paths.append(path)
    return paths


def table(rows: pd.DataFrame, mag_col: str) -> pd.DataFrame:
    """Decimal and sexagesimal positions with kind, for the plain/CSV/VOTable lists."""
    out = formatted(rows, mag_col).rename(columns={'ra': 'ra_hms', 'dec': 'dec_dms'})
    return out.assign(kind=rows['kind'].to_numpy(), ra=rows['ra'].round(6).to_numpy(),
                      dec=rows['dec'].round(6).to_numpy())


def plain(stem: Path, targets, stds, mag_col) -> list[Path]:
    """Whitespace-aligned text table of targets and standards in observing order."""
    path = stem.with_name(f'{stem.name}_list.txt')
    path.write_text(table(interleave(targets, stds), mag_col).to_string(index=False) + '\n')
    return [path]


def csv(stem: Path, targets, stds, mag_col) -> list[Path]:
    path = stem.with_name(f'{stem.name}_list.csv')
    table(interleave(targets, stds), mag_col).to_csv(path, index=False)
    return [path]


def votable(stem: Path, targets, stds, mag_col) -> list[Path]:
    """VOTable, e.g. to open in TOPCAT."""
    path = stem.with_name(f'{stem.name}_list.vot')
    Table.from_pandas(table(interleave(targets, stds), mag_col)).write(path, format='votable', overwrite=True)
    return [path]


def keck(stem: Path, targets, stds, mag_col) -> list[Path]:
    """Keck starlist (www2.keck.hawaii.edu/observing/starlist.html): name in columns 1-16 (max 15 characters),
    'hh mm ss.ss +dd mm ss.s 2000.0', then pmra= (s of RA/yr) and pmdec= (arcsec/yr) when non-zero.
    No comment or blank lines."""
    rows = interleave(targets, stds)
    names = clean_names(rows['shortname'], None, 'Keck', 15)
    ra, dec = sexagesimal(rows, sep=' ', ra_precision=2, dec_precision=1, sign=True)
    pmra, pmdec = time_pm(rows)
    lines = [f'{n:<16}{r} {d} 2000.0' + (f' pmra={a:+.5f} pmdec={b:+.4f}' if a or b else '')
             for n, r, d, a, b in zip(names, ra, dec, pmra, pmdec, strict=True)]
    path = stem.with_name(f'{stem.name}_keck.txt')
    path.write_text('\n'.join(lines) + '\n')
    return [path]


def lick(stem: Path, targets, stds, mag_col) -> list[Path]:
    """Lick (Shane) starlist (mthamilton.ucolick.org/techdocs/telescopes/starlistFormat.html): whitespace-separated
    'name h m s d m s 2000.0'; names without spaces. Proper motions omitted (the cos dec convention is unstated)."""
    rows = interleave(targets, stds)
    names = clean_names(rows['shortname'], None, 'Lick')
    ra, dec = sexagesimal(rows, sep=' ', ra_precision=2, dec_precision=1, sign=True)
    path = stem.with_name(f'{stem.name}_lick.txt')
    path.write_text('\n'.join(f'{n} {r} {d} 2000.0' for n, r, d in zip(names, ra, dec, strict=True))
                    + '\n')
    return [path]


def eso(stem: Path, targets, stds, mag_col) -> list[Path]:
    """ESO p2 'Import Target List' (eso.org/sci/observing/phase2/p2intro/p2-tutorials/p2-ImportTargetList.html):
    'name, HH:MM:SS.SSS, +DD:MM:SS.SS, J2000, epoch'; names of [a-zA-Z0-9.+- ] (underscores become spaces).
    Proper motions omitted (the cos dec convention is unstated)."""
    rows = interleave(targets, stds)
    names = clean_names(rows['shortname'].str.replace('_', ' '), r'[^a-zA-Z0-9.+\- ]', 'ESO p2')
    ra, dec = sexagesimal(rows, ra_precision=3, dec_precision=2, sign=True)
    path = stem.with_name(f'{stem.name}_eso_p2.txt')
    path.write_text('\n'.join(f'{n}, {r}, {d}, J2000, {e:.1f}'
                              for n, r, d, e in zip(names, ra, dec, rows['epoch'], strict=True)) + '\n')
    return [path]


def subaru(stem: Path, targets, stds, mag_col) -> list[Path]:
    """Subaru OPE target definitions for the <PARAMETER_LIST> section (subarutelescope.org/Instruments/HSC/ope.html):
    'ALIAS=OBJECT="name" RA=hhmmss.sss DEC=+ddmmss.ss EQUINOX=2000.0'; aliases of letters, digits and underscores."""
    rows = interleave(targets, stds)
    aliases = clean_names(rows['shortname'].str.replace(r'[+\-.]', '_', regex=True), r'[^A-Za-z0-9_]', 'Subaru')
    ra, dec = sexagesimal(rows, sep='', ra_precision=3, dec_precision=2, sign=True)
    path = stem.with_name(f'{stem.name}_subaru_ope.txt')
    path.write_text('\n'.join(f'{a}=OBJECT="{n}" RA={r} DEC={d} EQUINOX=2000.0'
                              for a, n, r, d in zip(aliases, rows['shortname'], ra, dec, strict=True)) + '\n')
    return [path]


def mmt(stem: Path, targets, stds, mag_col) -> list[Path]:
    """MMT catalogue (mmto.org/?p=544): 'ID hh:mm:ss.sss +dd:mm:ss.ss pmRA pmDec mag SpT epoch', blank-separated,
    no empty fields; proper motions in s of RA and arcsec per century."""
    rows = interleave(targets, stds)
    names = clean_names(rows['shortname'], None, 'MMT')
    ra, dec = sexagesimal(rows, ra_precision=3, dec_precision=2, sign=True)
    pmra, pmdec = time_pm(rows)
    mags = [f'{m:.2f}' if np.isfinite(m) else '-' for m in rows[mag_col]]
    spts = [s.replace(' ', '') or '-' for s in rows['spt']]
    lines = [f'{n} {r} {d} {a * 100:+.3f} {b * 100:+.2f} {m} {s} J{e:.1f}' for n, r, d, a, b, m, s, e in
             zip(names, ra, dec, pmra, pmdec, mags, spts, rows['epoch'], strict=True)]
    path = stem.with_name(f'{stem.name}_mmt.cat')
    path.write_text('\n'.join(lines) + '\n')
    return [path]


def magellan(stem: Path, targets, stds, mag_col) -> list[Path]:
    """Magellan catalogue (lco.cl/observing-catalogs-format/): 16 whitespace-separated fields: number, name, RA,
    Dec, equinox, pmRA (s/yr), pmDec (arcsec/yr), rotator offset angle and mode (0.0 OFF here: set them for your
    instrument), two guide-probe RA/Dec/equinox triples (zeros: none chosen), proper-motion epoch."""
    rows = interleave(targets, stds)
    names = clean_names(rows['shortname'], r'[^A-Za-z0-9+\-._]', 'Magellan')
    ra, dec = sexagesimal(rows, ra_precision=2, dec_precision=1, sign=True)
    pmra, pmdec = time_pm(rows)
    probes = '00:00:00.00 +00:00:00.0 2000.0 00:00:00.00 +00:00:00.0 2000.0'
    lines = ['# No. Name RA Dec Equinox pmRA pmDec RotOff RotMode GP1RA GP1Dec GP1Eq GP2RA GP2Dec GP2Eq pmEpoch']
    lines += [f'{i} {n} {r} {d} 2000.0 {a:.4f} {b:.3f} 0.0 OFF {probes} {e:.1f}' for i, (n, r, d, a, b, e) in
              enumerate(zip(names, ra, dec, pmra, pmdec, rows['epoch'], strict=True), 1)]
    path = stem.with_name(f'{stem.name}_magellan.cat')
    path.write_text('\n'.join(lines) + '\n')
    return [path]


FORMATS: dict[str, tuple[str, Writer]] = {
    'irtf': ('IRTF starcat', irtf),
    'cat': ('ING object visibility (.cat)', ing_cat),
    'txt': ('Plain text table', plain),
    'csv': ('CSV', csv),
    'votable': ('VOTable (TOPCAT)', votable),
    'keck': ('Keck starlist', keck),
    'lick': ('Lick starlist', lick),
    'eso': ('ESO p2 target import', eso),
    'subaru': ('Subaru OPE targets', subaru),
    'mmt': ('MMT catalogue', mmt),
    'magellan': ('Magellan catalogue', magellan),
}
SITE_FORMATS = {  # bundled observatory name -> the list format its telescopes read
    'NASA Infrared Telescope Facility': ('irtf',),
    'W. M. Keck Observatory': ('keck',),
    'Lick Observatory': ('lick',),
    'Paranal Observatory (ESO)': ('eso',),
    'La Silla Observatory (ESO)': ('eso',),
    'Subaru Telescope': ('subaru',),
    'Multiple Mirror Telescope': ('mmt',),
    'Las Campanas Observatory': ('magellan',),
    'Roque de los Muchachos, La Palma': ('cat',),
}


def defaults(site: str | None = None, from_samp: bool = False) -> tuple[str, ...]:
    """Formats to write when none are chosen: the observatory's own list format, else a VOTable for a table
    that came from TOPCAT (SAMP), else CSV. ``site`` is a bundled observatory name or alias."""
    if site:
        name = find_site(site)['name']
        if name in SITE_FORMATS:
            return SITE_FORMATS[name]
    return ('votable',) if from_samp else ('csv',)


def write(keys, stem: Path, targets: pd.DataFrame, stds: pd.DataFrame | None, mag_col: str) -> list[Path]:
    """Write the formats named in ``keys`` for ``stem`` (an output folder plus file stem)."""
    unknown = sorted(set(keys) - set(FORMATS))
    if unknown:
        raise ValueError(f'Unknown formats {unknown}; choose from {list(FORMATS)}')
    return [p for key in keys for p in FORMATS[key][1](stem, targets, stds, mag_col)]
