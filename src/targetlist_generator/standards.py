"""Telluric standards: load the bundled (or a user) list and pick one near each target on the sky and in airmass."""
from importlib.resources import files
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from . import sky
from .io import normalise

BUNDLED = files(__package__) / 'data' / 'large_standards_bright.csv'  # B9V-A1V, K = 5-9 mag, 4969 stars
BANDS = ('B', 'V', 'G', 'J', 'H', 'K')


def load_standards(mag_col: str, band: str = 'K', spt: str | None = None, path: str | Path | None = None):
    """Standards with ``band`` renamed to ``mag_col`` (to sit in the target list), optionally one spectral type."""
    raw = pd.read_csv(path or BUNDLED)
    if band not in raw.columns:
        raise ValueError(f'Standards have no band {band}; available: {[c for c in BANDS if c in raw.columns]}')
    raw = raw[['shortname', 'ra', 'dec', 'spt', band]].rename(columns={band: mag_col})
    if spt:
        raw = raw[raw['spt'].str.strip() == spt]
        if raw.empty:
            raise ValueError(f'No standards of spectral type {spt}')
    return normalise(raw, mag_col, source=str(path or BUNDLED))


def pick(targets: pd.DataFrame, stds: pd.DataFrame, times, location, std_at: np.ndarray, target_am: np.ndarray,
         tol: float = 0.1, max_sep: float = 20.0) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Standard per target, observed at grid index ``std_at``: nearest on the sky among those within ``max_sep``
    deg and within ``tol`` of the target's airmass ``target_am``; else the closest in airmass within ``max_sep``;
    else (none within ``max_sep`` above the horizon) the nearest on the sky, noted. Small slews beat perfect
    airmass matches, so a standard is never taken from across the sky while a nearby one exists.
    Returns (standards, notes, standard airmass, separation [deg])."""
    s_am = sky.airmass(stds['ra'].to_numpy(), stds['dec'].to_numpy(), times, location)[std_at]  # (N, S)
    with np.errstate(invalid='ignore'):  # inf - inf below the horizon -> nan, never a candidate
        diff = np.nan_to_num(np.abs(s_am - np.asarray(target_am)[:, None]), nan=np.inf)
    sep = coords(targets)[:, None].separation(coords(stds)[None, :]).deg
    near = (sep <= max_sep) & np.isfinite(diff)
    match = near & (diff <= tol)
    choice = np.argmin(np.where(match, sep, np.inf), axis=1)
    loose = ~match.any(axis=1)
    choice[loose] = np.argmin(np.where(near, diff, np.inf), axis=1)[loose]
    far = ~near.any(axis=1)
    choice[far] = np.argmin(sep, axis=1)[far]
    notes = np.where(far, f'no standard above the horizon within {max_sep:g} deg',
                     np.where(loose, f'no standard within {tol:g} airmass', ''))
    rows = np.arange(len(targets))
    return stds.iloc[choice].reset_index(drop=True), notes, s_am[rows, choice], sep[rows, choice]


def coords(df: pd.DataFrame) -> SkyCoord:
    return SkyCoord(df['ra'].to_numpy() * u.deg, df['dec'].to_numpy() * u.deg)


def nearest_standards(targets: pd.DataFrame, standards: pd.DataFrame, offset: float = 30.0) -> pd.DataFrame:
    """Row per target: the standard nearest on the sky to the target's position ``offset`` minutes later in RA,
    i.e. one that will be at a similar airmass when observed straight after the target."""
    shifted = SkyCoord((targets['ra'].to_numpy() + offset / 4) * u.deg, targets['dec'].to_numpy() * u.deg)
    stds = SkyCoord(standards['ra'].to_numpy() * u.deg, standards['dec'].to_numpy() * u.deg)
    nearest = np.argmin(shifted[:, None].separation(stds[None, :]), axis=1)  # no scipy needed for a few 1e6 pairs
    return standards.iloc[nearest].reset_index(drop=True)
