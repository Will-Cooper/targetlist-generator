"""Telluric standards: load the bundled (or a user) list and pair each target with its nearest standard."""
from importlib.resources import files
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

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


def nearest_standards(targets: pd.DataFrame, standards: pd.DataFrame, offset: float = 30.0) -> pd.DataFrame:
    """Row per target: the standard nearest on the sky to the target's position ``offset`` minutes later in RA,
    i.e. one that will be at a similar airmass when observed straight after the target."""
    shifted = SkyCoord((targets['ra'].to_numpy() + offset / 4) * u.deg, targets['dec'].to_numpy() * u.deg)
    stds = SkyCoord(standards['ra'].to_numpy() * u.deg, standards['dec'].to_numpy() * u.deg)
    nearest = np.argmin(shifted[:, None].separation(stds[None, :]), axis=1)  # no scipy needed for a few 1e6 pairs
    return standards.iloc[nearest].reset_index(drop=True)
