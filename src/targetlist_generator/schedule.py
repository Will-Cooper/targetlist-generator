"""Per-target timing for one night, the Moon, and the airmass plot."""
from pathlib import Path

import astropy.units as u
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.utils import iers
from matplotlib import colormaps, rc_context
from matplotlib.figure import Figure

from . import sky


def moon(times: Time, location) -> tuple[SkyCoord, float]:
    """Moon position at mid-night and its illuminated fraction (from the Sun-Moon elongation)."""
    mid = times[len(times) // 2]
    with iers.conf.set_temp('auto_download', False):
        body = get_body('moon', mid, location)
        elongation = body.separation(get_sun(mid).transform_to(body.frame))
    return body, float((1 - np.cos(elongation)) / 2)


def moon_separation(targets: pd.DataFrame, body: SkyCoord) -> np.ndarray:
    """Target-Moon separations [deg], computed in the Moon's (topocentric) frame."""
    stars = SkyCoord(targets['ra'].to_numpy() * u.deg, targets['dec'].to_numpy() * u.deg)
    with iers.conf.set_temp('auto_download', False):
        return body.separation(stars.transform_to(body.frame)).deg


def timing(targets: pd.DataFrame, times: Time, location, max_airmass: float) -> pd.DataFrame:
    """best_utc (minimum airmass), min_airmass and the first/last UTC at airmass <= ``max_airmass``."""
    am = sky.airmass(targets['ra'].to_numpy(), targets['dec'].to_numpy(), times, location)
    up = am <= max_airmass
    iso = np.array([t[:16] for t in times.iso])
    first, last = up.argmax(axis=0), len(times) - 1 - up[::-1].argmax(axis=0)
    never = ~up.any(axis=0)
    return pd.DataFrame({'best_utc': iso[am.argmin(axis=0)], 'min_airmass': am.min(axis=0).round(3),
                         'start_utc': np.where(never, '', iso[first]), 'end_utc': np.where(never, '', iso[last])})


def airmass_plot(targets: pd.DataFrame, times: Time, location, max_airmass: float, path: Path) -> Path:
    """Airmass against UTC for every target, coloured by time of best airmass, names at their minima."""
    am = sky.airmass(targets['ra'].to_numpy(), targets['dec'].to_numpy(), times, location)
    rank = np.argsort(np.argsort(am.argmin(axis=0), kind='stable'), kind='stable')  # by time of best airmass
    colours = colormaps['viridis'](np.linspace(0, 0.9, am.shape[1]))[rank]
    am[am > max_airmass + 0.5] = np.nan
    when = times.datetime
    with rc_context({'font.family': 'serif'}):
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot()
        for i, name in enumerate(targets['shortname']):
            ax.plot(when, am[:, i], color=colours[i], lw=1)
            if np.isfinite(am[:, i]).any():
                j = np.nanargmin(am[:, i])
                ax.annotate(name, (when[j], am[j, i]), xytext=(0, -3), textcoords='offset points',
                            ha='center', va='top', fontsize=6, color=colours[i])
        ax.axhline(max_airmass, color='black', ls='--', lw=1)
        ax.set_ylim(max_airmass + 0.3, 0.95)
        ax.set_xlim(when[0], when[-1])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.set_xlabel(f'UTC from {times[0].iso[:10]} [h]')
        ax.set_ylabel('Airmass')
        fig.savefig(path, dpi=150, bbox_inches='tight')
    return path
