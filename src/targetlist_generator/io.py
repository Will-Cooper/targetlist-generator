"""Read target tables and write IRTF, no-standards and ING ``.cat`` target lists."""
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

REQUIRED = ('shortname', 'ra', 'dec', 'spt')
RENAME = {'name': 'shortname', 'median_sptname': 'spt', 'gaiag_spt': 'spt'}  # first match wins
TXT_NAMES = ['index', 'shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt']  # IRTF list columns, mag appended
CAT_CHUNK = 10  # targets per ING .cat file
RULE = '#+-----+-----+-----+-----+-----+-----+-----+-----+-------+'


def read_table(path: str | Path, mag_col: str = 'tmassk') -> pd.DataFrame:
    """Read a CSV, IRTF ``.txt`` list or astropy-readable table (FITS, ECSV, VOTable); lower-case names."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Cannot find target list {path}')
    suffix = path.suffix.lower()
    if suffix == '.csv':
        df = pd.read_csv(path)
    elif suffix in ('.txt', '.tsv'):  # our own IRTF output, so lists can be re-cut
        df = pd.read_csv(path, sep='\t', comment='#', header=None, names=TXT_NAMES + [mag_col.lower()])
        df = df.drop(columns='index')
    else:
        df = Table.read(path).to_pandas()
        for col in df.select_dtypes(object).columns:  # FITS strings arrive as bytes
            df[col] = [v.decode().strip() if isinstance(v, bytes) else v for v in df[col]]
    df = df.rename(columns=str.lower)
    for src, dst in RENAME.items():
        if src in df.columns and dst not in df.columns:
            df = df.rename(columns={src: dst})
    return df


def spt_from_num(num: float) -> str:
    """Spectral type number (M0 = 60, L0 = 70, T0 = 80, Y0 = 90) to string, rounded to half a subtype."""
    half = 0.5 * round(num / 0.5)
    if not 60 <= half < 100:
        raise ValueError(f'Spectral type number {num} outside M0-Y9.5 (60-99.5)')
    return f"{'MLTY'[int(half // 10) - 6]}{half % 10:g}"


def read_targets(path: str | Path, mag_col: str = 'tmassk') -> pd.DataFrame:
    """Read a target list into columns shortname, ra, dec (deg), epoch, pmra, pmdec, spt and ``mag_col``."""
    return normalise(read_table(path, mag_col), mag_col, source=str(path))


def normalise(df: pd.DataFrame, mag_col: str, source: str = 'table') -> pd.DataFrame:
    """Check required columns and convert to the standard column set (see ``read_targets``)."""
    df = df.rename(columns=str.lower)
    mag_col = mag_col.lower()
    missing = [c for c in (*REQUIRED, mag_col) if c not in df.columns]
    if missing:
        raise ValueError(f'{source} is missing columns {missing} (case-insensitive; "name" also accepted for '
                         f'"shortname"); it has {list(df.columns)}')
    df['shortname'] = df['shortname'].astype(str).str.strip().str.replace(r'\s+', '_', regex=True)  # split on spaces
    for col, unit in (('ra', u.hourangle), ('dec', u.deg)):  # numbers are degrees; sexagesimal strings in unit
        deg = pd.to_numeric(df[col], errors='coerce')
        text = deg.isna() & df[col].notna()
        if text.any():
            deg[text] = Angle(df.loc[text, col].astype(str).str.strip().to_numpy(), unit=unit).deg
        if deg.isna().any():
            raise ValueError(f'{source}: no {col} for {list(df.loc[deg.isna(), "shortname"])}')
        df[col] = deg
    for col, default in (('epoch', 2000.0), ('pmra', 0.0), ('pmdec', 0.0)):
        df[col] = pd.to_numeric(df[col]).fillna(default) if col in df.columns else default
    sptnum = pd.to_numeric(df['spt'], errors='coerce')  # numeric types (L0 = 70) may be mixed with strings
    df['spt'] = [spt_from_num(n) if np.isfinite(n) else ('' if pd.isna(s) else str(s).strip())
                 for n, s in zip(sptnum, df['spt'], strict=True)]
    return df[['shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', mag_col]].reset_index(drop=True)


def formatted(df: pd.DataFrame, mag_col: str) -> pd.DataFrame:
    """Sexagesimal, rounded copy of a target table for writing."""
    # round to the written precision first so 23:59:59.999 wraps to 00:00:00.00, never 24:00 or 60 s
    ra = (np.round(df['ra'].to_numpy() * 240, 2) % 86400) / 240  # 240 s of time per degree
    dec = np.round(df['dec'].to_numpy() * 3600, 1) / 3600
    coords = SkyCoord(ra * u.deg, dec * u.deg)
    return pd.DataFrame({
        'shortname': df['shortname'].to_numpy(),
        'ra': coords.ra.to_string(u.hourangle, sep=':', precision=2, pad=True),
        'dec': coords.dec.to_string(u.deg, sep=':', precision=1, pad=True),
        'epoch': [f'{e:g}' for e in df['epoch']],
        'pmra': df['pmra'].round(2).to_numpy(),
        'pmdec': df['pmdec'].round(2).to_numpy(),
        'spt': df['spt'].to_numpy(),
        mag_col: df[mag_col].fillna(0).round(2).to_numpy(),  # IRTF list needs a number; 0 marks unknown
    })


def write_irtf(path: Path, targets: pd.DataFrame, mag_col: str, standards: pd.DataFrame | None = None) -> Path:
    """IRTF starcat list; with ``standards`` (one row per target) each target is followed by its standard."""
    rows = formatted(targets, mag_col)
    if standards is not None:
        rows = pd.concat([rows, formatted(standards, mag_col)]).sort_index(kind='stable')  # t0, s0, t1, s1, ...
    rows.index = np.arange(1, len(rows) + 1)
    header = f'#|Index| Name|  RA | DEC |Epoch| PMRA|PMDEC|SpT  | {mag_col.upper():>6}|'
    path.write_text('\n'.join((RULE, header, RULE, '')) + rows.to_csv(sep='\t', header=False, lineterminator='\n'))
    return path


def write_cats(stem: Path, targets: pd.DataFrame, mag_col: str) -> list[Path]:
    """ING object-visibility ``.cat`` files, ``CAT_CHUNK`` targets each: name ra dec epoch."""
    rows = formatted(targets, mag_col)[['shortname', 'ra', 'dec', 'epoch']]
    paths = []
    for k, start in enumerate(range(0, len(rows), CAT_CHUNK), 1):
        path = stem.with_name(f'{stem.name}_cat_{k}.cat')
        rows.iloc[start:start + CAT_CHUNK].to_csv(path, sep=' ', header=False, index=False, lineterminator='\n')
        paths.append(path)
    return paths
