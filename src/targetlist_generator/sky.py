"""Night window and target visibility from an observatory site, plus plain RA/Dec box cuts."""
from difflib import get_close_matches
from functools import cache
from importlib.resources import files

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import FK5, AltAz, Angle, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astropy.utils import iers

STEP = 5 * u.min  # time resolution of the night grid
OBSERVATORIES = files(__package__) / 'data' / 'observatories.csv'  # snapshot of astropy's site registry


@cache
def observatories() -> pd.DataFrame:
    """Bundled observatories: name, aliases (';'-separated), lat, lon [deg, east +], height [m]."""
    return pd.read_csv(OBSERVATORIES, keep_default_na=False)


def find_site(site: str) -> pd.Series:
    """Observatory row by name or alias, case-insensitive; unknown names raise with the closest matches."""
    sites = observatories()
    keys = {n.lower(): i for i, row in sites.iterrows()
            for n in [row['name'], *filter(None, row['aliases'].split(';'))]}
    if site.strip().lower() not in keys:
        close = get_close_matches(site.strip().lower(), keys, n=5, cutoff=0.5)
        raise ValueError(f'Unknown observatory {site!r}' + (f'; did you mean {close}?' if close else '')
                         + ' Give --lat/--lon instead for a site not in the list.')
    return sites.loc[keys[site.strip().lower()]]


def site_location(site: str | None = None, lat: float | None = None, lon: float | None = None) -> EarthLocation:
    """Observatory from a bundled site name or alias (e.g. 'lapalma', 'irtf') or decimal-degree latitude/longitude."""
    if site and (lat is not None or lon is not None):
        raise ValueError('Give either a site name or latitude/longitude, not both')
    if site:
        row = find_site(site)
        return EarthLocation.from_geodetic(row['lon'] * u.deg, row['lat'] * u.deg, row['height'] * u.m)
    if lat is None or lon is None:
        raise ValueError('Need a site name, or both latitude and longitude (decimal degrees, east positive)')
    return EarthLocation.from_geodetic(lon * u.deg, lat * u.deg)


def day_grid(date: str, location: EarthLocation) -> tuple[Time, np.ndarray]:
    """``STEP`` grid over the 24 h from local noon of ``date``, and the Sun's altitude [deg] on it."""
    local_noon = Time(date) + (12 - location.lon.to_value(u.deg) / 15) * u.h
    grid = local_noon + np.arange(0, 24 * 60, STEP.to_value(u.min)) * u.min
    with iers.conf.set_temp('auto_download', False):  # minute precision is plenty; works offline at the telescope
        sun_alt = get_sun(grid).transform_to(AltAz(obstime=grid, location=location)).alt.deg
    return grid, sun_alt


def night_times(date: str, location: EarthLocation, twilight: float = -12.0) -> Time:
    """Times (``STEP`` apart) during the night starting on the evening of ``date`` with the Sun below ``twilight``."""
    grid, sun_alt = day_grid(date, location)
    night = grid[sun_alt < twilight]
    if len(night) == 0:
        raise ValueError(f'The Sun never sets below {twilight} deg on {date} at this site')
    return night


def altitude(ra: np.ndarray, dec: np.ndarray, times: Time, location: EarthLocation) -> np.ndarray:
    """Geometric altitude [deg], shape (len(times), len(ra)); ra, dec in degrees (J2000, precessed to date)."""
    with iers.conf.set_temp('auto_download', False):
        lst = times.sidereal_time('apparent', longitude=location.lon).deg
    now = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg).transform_to(FK5(equinox=times[0]))
    hour_angle = np.radians(lst[:, None] - now.ra.deg[None, :])
    lat, dec = location.lat.rad, now.dec.rad[None, :]
    return np.degrees(np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle)))


def airmass(ra: np.ndarray, dec: np.ndarray, times: Time, location: EarthLocation) -> np.ndarray:
    """Plane-parallel airmass (inf below the horizon), shape (len(times), len(ra)); ra, dec in degrees."""
    sin_alt = np.sin(np.radians(altitude(ra, dec, times, location)))
    with np.errstate(divide='ignore'):
        return np.where(sin_alt > 0, 1 / sin_alt, np.inf)


def visible(ra: np.ndarray, dec: np.ndarray, times: Time, location: EarthLocation,
            max_airmass: float = 2.0, min_time: float = 30.0) -> np.ndarray:
    """True where a target spends at least ``min_time`` minutes of ``times`` below ``max_airmass``."""
    minutes = (airmass(ra, dec, times, location) <= max_airmass).sum(axis=0) * STEP.to_value(u.min)
    return minutes >= min_time


def parse_angle(value: str | float, unit: u.Unit) -> float:
    """Degrees from a number (degrees) or a sexagesimal string ('20h', '20:30:00' in ``unit``, '-15d30m')."""
    try:
        return float(value)
    except ValueError:
        return Angle(str(value).strip(), unit=unit).deg


def in_box(ra: np.ndarray, dec: np.ndarray, ra_range: tuple | None = None, dec_range: tuple | None = None):
    """True inside RA (wrapping through 0 when lower > upper) and Dec limits, all inclusive."""
    keep = np.ones(len(ra), dtype=bool)
    if ra_range is not None:
        lo, hi = (parse_angle(v, u.hourangle) for v in ra_range)
        keep &= ((ra >= lo) & (ra <= hi)) if lo <= hi else ((ra >= lo) | (ra <= hi))
    if dec_range is not None:
        lo, hi = sorted(parse_angle(v, u.deg) for v in dec_range)
        keep &= (dec >= lo) & (dec <= hi)
    return keep
