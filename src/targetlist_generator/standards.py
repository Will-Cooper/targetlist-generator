"""Telluric standards: load the bundled (or a user) list and pair each target with its nearest standard."""
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


def matched_standards(targets: pd.DataFrame, standards: pd.DataFrame, times, location, offset: float = 30.0,
                      tol: float = 0.1) -> tuple[pd.DataFrame, int]:
    """Row per target: among standards whose airmass ``offset`` minutes after the target's best (minimum-airmass)
    time is within ``tol`` of the target's, the nearest on the sky; with none within ``tol``, the closest in
    airmass. Also returns how many targets had no standard within ``tol``."""
    t_am = sky.airmass(targets['ra'].to_numpy(), targets['dec'].to_numpy(), times, location)  # (T, N)
    s_am = sky.airmass(standards['ra'].to_numpy(), standards['dec'].to_numpy(), times, location)  # (T, S)
    best = t_am.argmin(axis=0)
    later = np.clip(best + round(offset / sky.STEP.to_value(u.min)), 0, len(times) - 1)  # stays inside the night
    with np.errstate(invalid='ignore'):  # inf - inf below the horizon -> nan, never a candidate
        diff = np.abs(s_am[later] - t_am[best, np.arange(len(targets))][:, None])  # (N, S)
    diff = np.nan_to_num(diff, nan=np.inf)
    sep = coords(targets)[:, None].separation(coords(standards)[None, :]).deg
    choice = np.argmin(np.where(diff <= tol, sep, np.inf), axis=1)
    none = ~(diff <= tol).any(axis=1)
    choice[none] = np.argmin(diff[none], axis=1)
    return standards.iloc[choice].reset_index(drop=True), int(none.sum())


def coords(df: pd.DataFrame) -> SkyCoord:
    return SkyCoord(df['ra'].to_numpy() * u.deg, df['dec'].to_numpy() * u.deg)


def nearest_standards(targets: pd.DataFrame, standards: pd.DataFrame, offset: float = 30.0) -> pd.DataFrame:
    """Row per target: the standard nearest on the sky to the target's position ``offset`` minutes later in RA,
    i.e. one that will be at a similar airmass when observed straight after the target."""
    shifted = SkyCoord((targets['ra'].to_numpy() + offset / 4) * u.deg, targets['dec'].to_numpy() * u.deg)
    stds = SkyCoord(standards['ra'].to_numpy() * u.deg, standards['dec'].to_numpy() * u.deg)
    nearest = np.argmin(shifted[:, None].separation(stds[None, :]), axis=1)  # no scipy needed for a few 1e6 pairs
    return standards.iloc[nearest].reset_index(drop=True)
