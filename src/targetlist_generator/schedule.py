"""The Moon and the night's visibility plots."""
from pathlib import Path

import astropy.units as u
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.utils import iers
from matplotlib import rc_context
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


def moon_limit(setting: str | float, illumination: float) -> float:
    """Minimum Moon separation [deg]: 'auto' = 10 + 30 x illumination (10 at new, 40 at full), 'off' = 0, or a
    number of degrees."""
    if setting == 'auto':
        return 10 + 30 * illumination
    if setting == 'off':
        return 0.0
    value = float(setting)
    if not 0 <= value <= 180:
        raise ValueError(f'Moon separation must be auto, off or 0-180 deg, not {setting!r}')
    return value


PER_PLOT = 10  # targets per visibility plot, as in the ING .cat lists
COLOURS = ('black', '#004488', '#DDAA33')  # colour-blind safe, then cycled through line styles
STYLES = ('-', '--', '-.', ':')


def visibility_plots(plan: pd.DataFrame, date: str, location, max_airmass: float, title: str,
                     stem: Path) -> list[Path]:
    """ING-style visibility plots, ``PER_PLOT`` targets each (``<stem>_visibility_N.png``): altitude against UTC
    from sunset to sunrise with an airmass axis, twilight shading (darker grey = brighter sky), the Moon (grey
    dotted) and each target's planned block drawn thick."""
    grid, sun = sky.day_grid(date, location)
    up = np.flatnonzero(sun < 0)
    span = slice(max(up[0] - 6, 0), min(up[-1] + 7, len(grid)))  # 30 min either side of sunset/sunrise
    times, sun = grid[span], sun[span]
    when = times.datetime
    iso = np.array([t[:16] for t in times.iso])
    with iers.conf.set_temp('auto_download', False):
        moon_alt = get_body('moon', times, location).transform_to(AltAz(obstime=times, location=location)).alt.deg
    alt_all = sky.altitude(plan['ra'].to_numpy(), plan['dec'].to_numpy(), times, location)
    limit_alt = np.degrees(np.arcsin(1 / max_airmass))
    paths = []
    for k, start in enumerate(range(0, len(plan), PER_PLOT), 1):
        chunk = plan.iloc[start:start + PER_PLOT]
        with rc_context({'font.family': 'serif'}):
            fig = Figure(figsize=(9, 5.5))
            ax = fig.add_subplot()
            for level in (-18, -12, -6, 0):
                ax.fill_between(when, 0, 90, where=sun > level, color='grey', alpha=0.12, lw=0)
            ax.plot(when, moon_alt, color='grey', ls=':', lw=1.5, label='Moon')
            ax.axhline(limit_alt, color='grey', ls='--', lw=0.8)
            for j, row in enumerate(chunk.itertuples()):
                style = {'color': COLOURS[j % 3], 'ls': STYLES[j // 3 % 4]}
                alt = alt_all[:, start + j]
                ax.plot(when, alt, lw=0.9, **style)
                block = (iso >= row.target_start) & (iso <= row.target_end)
                ax.plot(when[block], alt[block], lw=3.5, solid_capstyle='butt', **style)
                if block.any():
                    m = np.flatnonzero(block)[block.sum() // 2]
                    ax.annotate(row.shortname, (when[m], alt[m]), xytext=(0, 6), textcoords='offset points',
                                ha='center', fontsize=7, color=style['color'])
            ax.set_ylim(0, 90)
            ax.set_xlim(when[0], when[-1])
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            ax.minorticks_on()
            ax.tick_params(which='both', direction='in', top=True)
            to_am = (lambda a: 1 / np.sin(np.radians(np.clip(a, 1, 90))),
                     lambda x: np.degrees(np.arcsin(1 / np.clip(x, 1, None))))
            right = ax.secondary_yaxis('right', functions=to_am)
            right.set_ticks([1.0, 1.1, 1.2, 1.5, 2.0, 3.0])
            right.tick_params(direction='in')
            right.set_ylabel('Airmass')
            ax.set_xlabel(f'UTC from {times[0].iso[:10]} [h]')
            ax.set_ylabel('Altitude [deg]')
            ax.set_title(f'{title}, night of {date}' + (f' ({k})' if len(plan) > PER_PLOT else ''))
            ax.legend(loc='upper right', fontsize=7, frameon=False)
            path = stem.with_name(f'{stem.name}_visibility_{k}.png')
            fig.savefig(path, dpi=150, bbox_inches='tight')
        paths.append(path)
    return paths
