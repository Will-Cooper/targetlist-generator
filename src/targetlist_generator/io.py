"""Read target tables (CSV, TSV, text, FITS, VOTable, ECSV, IRTF lists) into one standard set of columns."""
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

REQUIRED = ('shortname', 'ra', 'dec')
TXT_NAMES = ['index', 'shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt']  # IRTF list columns, mag appended
IRTF_RULE = '#+-----'
TEXT_SUFFIXES = ('.txt', '.dat', '.tbl', '.ascii', '.text')
FILE_TYPES = ('CSV, TSV, text tables with a header line (whitespace or IPAC .tbl), FITS (also gzipped), '
              'VOTable (.vot/.xml), ECSV, and this tool\'s own IRTF lists')
# Standard column -> names recognised for it (lower case, first match wins); the user can always choose others.
ALIASES = {
    'shortname': ('shortname', 'name', 'designation', 'object', 'target', 'main_id', 'id', 'source_id', 'star'),
    'ra': ('ra', 'raj2000', '_raj2000', 'ra_icrs', 'ra_deg', 'radeg', 'ra2000', 'ra_j2000', 'alpha'),
    'dec': ('dec', 'dej2000', '_dej2000', 'de_icrs', 'dec_deg', 'decdeg', 'decj2000', 'dec_j2000', 'de', 'delta'),
    'spt': ('spt', 'sptype', 'sp_type', 'spectral_type', 'spectype', 'median_sptname', 'gaiag_spt', 'spt_opt',
            'spt_ir'),
    'pmra': ('pmra', 'pm_ra', 'pmra_cosdec', 'pmracosdec'),
    'pmdec': ('pmdec', 'pm_dec', 'pmde'),
    'epoch': ('epoch', 'ref_epoch'),
}
MAGNITUDES = ('tmassk', 'kmag', 'ks', 'ksmag', 'k', 'k_mag', 'hmag', 'h', 'jmag', 'j', 'phot_g_mean_mag', 'gmag',
              'g', 'vmag', 'v')


def read_table(path: str | Path) -> pd.DataFrame:
    """Read any table in ``FILE_TYPES`` with lower-case column names, unchanged otherwise."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Cannot find table {path}')
    suffix = path.suffix.lower()
    if suffix == '.csv':
        df = pd.read_csv(path)
    elif suffix == '.tsv':
        df = pd.read_csv(path, sep='\t')
    elif suffix in TEXT_SUFFIXES and path.read_text(errors='replace').startswith(IRTF_RULE):  # an IRTF list
        mag = path.read_text().splitlines()[1].split('|')[-2].strip().lower()
        df = pd.read_csv(path, sep='\t', comment='#', header=None, names=TXT_NAMES + [mag]).drop(columns='index')
    else:
        df = Table.read(path, format='ascii') if suffix in TEXT_SUFFIXES else Table.read(path)
        df = df.to_pandas()
        for col in [c for c in df.columns if pd.api.types.is_object_dtype(df[c])]:  # FITS strings arrive as bytes
            df[col] = [v.decode().strip() if isinstance(v, bytes) else v for v in df[col]]
    df = df.rename(columns=lambda c: str(c).strip().lower())
    clash = df.columns[df.columns.duplicated()].tolist()
    if clash:
        raise ValueError(f'{path.name} has columns that differ only in case: {clash}')
    return df


def guess_columns(columns) -> dict[str, str | None]:
    """Best guess at which column holds each standard field (``ALIASES``), or None."""
    names = list(columns)
    return {field: next((n for n in aliases if n in names), None) for field, aliases in ALIASES.items()}


def magnitude_column(columns, mag_col: str | None = None) -> str:
    """The chosen magnitude column, else the first of ``MAGNITUDES`` present."""
    if mag_col:
        return mag_col.lower()
    guess = next((n for n in MAGNITUDES if n in list(columns)), None)
    if guess is None:
        raise ValueError(f'No magnitude column recognised: choose one (--mag-col, or in the window) '
                         f'from {list(columns)}')
    return guess


def apply_columns(df: pd.DataFrame, columns: dict[str, str | None] | None = None) -> pd.DataFrame:
    """Rename the columns holding each standard field to its standard name: the best guess, overridden by
    ``columns`` ({field: column}, with 'name' accepted for 'shortname' and None to ignore a field)."""
    choice = guess_columns(df.columns)
    for field, col in (columns or {}).items():
        field = 'shortname' if field == 'name' else field
        if field not in ALIASES:
            raise ValueError(f'Unknown field {field!r}; choose from {["name", *list(ALIASES)[1:]]}')
        if col and col.lower() not in df.columns:
            raise ValueError(f'No column {col!r} for {field}; columns are {list(df.columns)}')
        choice[field] = col.lower() if col else None
    used = {c for c in choice.values() if c}
    df = df.drop(columns=[f for f, c in choice.items() if f in df.columns and c != f and f not in used])
    return df.rename(columns={c: f for f, c in choice.items() if c and c != f})


def spt_from_num(num: float) -> str:
    """Spectral type number (M0 = 60, L0 = 70, T0 = 80, Y0 = 90) to string, rounded to half a subtype."""
    half = 0.5 * round(num / 0.5)
    if not 60 <= half < 100:
        raise ValueError(f'Spectral type number {num} outside M0-Y9.5 (60-99.5)')
    return f"{'MLTY'[int(half // 10) - 6]}{half % 10:g}"


def read_targets(path: str | Path, mag_col: str | None = None, keep: tuple[str, ...] = (),
                 columns: dict[str, str | None] | None = None) -> pd.DataFrame:
    """Read a target list into columns shortname, ra, dec (deg), epoch, pmra, pmdec (mas/yr, pmra including
    cos dec as in Gaia), spt, the magnitude column and any ``keep`` columns. Fields come from the best-guess
    columns unless ``columns`` says otherwise (see ``apply_columns``); the magnitude column (guessed when
    ``mag_col`` is None) is named in ``df.attrs['mag_col']``."""
    raw = read_table(path)
    mag_col = magnitude_column(raw.columns, mag_col)
    df = normalise(apply_columns(raw, columns), mag_col, source=str(path), keep=keep)
    df.attrs['mag_col'] = mag_col
    return df


def normalise(df: pd.DataFrame, mag_col: str, source: str = 'table', keep: tuple[str, ...] = ()) -> pd.DataFrame:
    """Check required columns and convert to the standard column set (see ``read_targets``)."""
    df = df.rename(columns=str.lower)
    mag_col = mag_col.lower()
    keep = tuple(c.lower() for c in keep)
    missing = [c for c in (*REQUIRED, mag_col, *keep) if c not in df.columns]
    if missing:
        raise ValueError(f'{source} has no column for {missing}: choose the right columns, or add them; it has '
                         f'{list(df.columns)}')
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
    if 'spt' not in df.columns:  # optional: written into the lists when known
        df['spt'] = ''
    sptnum = pd.to_numeric(df['spt'], errors='coerce')  # numeric types (L0 = 70) may be mixed with strings
    df['spt'] = [spt_from_num(n) if np.isfinite(n) else ('' if pd.isna(s) else str(s).strip())
                 for n, s in zip(sptnum, df['spt'], strict=True)]
    columns = ['shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', mag_col]
    if 'priority' in df.columns:  # optional: 1 (highest) to 5 (lowest), used by the planner
        columns.append('priority')
    return df[columns + [c for c in keep if c not in columns]].reset_index(drop=True)


def sexagesimal(df: pd.DataFrame, sep: str = ':', ra_precision: int = 2, dec_precision: int = 1,
                sign: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """RA (hours) and Dec strings, rounded to the written precision first so 23:59:59.999 wraps to 00:00:00,
    never 24:00 or 60 s; ``sign`` forces a '+' on positive declinations."""
    ra = (np.round(df['ra'].to_numpy() * 240, ra_precision) % 86400) / 240  # 240 s of time per degree
    dec = np.round(df['dec'].to_numpy() * 3600, dec_precision) / 3600
    coords = SkyCoord(ra * u.deg, dec * u.deg)
    return (coords.ra.to_string(u.hourangle, sep=sep, precision=ra_precision, pad=True),
            coords.dec.to_string(u.deg, sep=sep, precision=dec_precision, pad=True, alwayssign=sign))


def formatted(df: pd.DataFrame, mag_col: str) -> pd.DataFrame:
    """Sexagesimal, rounded copy of a target table for writing."""
    ra, dec = sexagesimal(df)
    return pd.DataFrame({
        'shortname': df['shortname'].to_numpy(),
        'ra': ra,
        'dec': dec,
        'epoch': [f'{e:g}' for e in df['epoch']],
        'pmra': df['pmra'].round(2).to_numpy(),
        'pmdec': df['pmdec'].round(2).to_numpy(),
        'spt': df['spt'].to_numpy(),
        mag_col: df[mag_col].fillna(0).round(2).to_numpy(),  # IRTF list needs a number; 0 marks unknown
    })
