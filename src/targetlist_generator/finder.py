"""Finder charts from CDS hips2fits cut-outs; Aladin on a HiPS mirror as the fallback when hips2fits is down."""
import shutil
import subprocess
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.figure import Figure  # pyplot-free, so safe in the GUI worker thread

HIPS2FITS = ('https://alasky.cds.unistra.fr/hips-image-services/hips2fits',  # tried in order
             'https://alaskybis.cds.unistra.fr/hips-image-services/hips2fits')
IRSA = 'https://irsa.ipac.caltech.edu/data/hips/CDS/'
MIRRORS = {'CDS/P/2MASS/K': IRSA + '2MASS/K', 'CDS/P/2MASS/J': IRSA + '2MASS/J', 'CDS/P/2MASS/H': IRSA + '2MASS/H',
           'CDS/P/DSS2/red': IRSA + 'DSS2/red', 'CDS/P/allWISE/W1': IRSA + 'AllWISE/W1'}  # HiPS tile copies
SURVEYS = (*MIRRORS, 'CDS/P/PanSTARRS/DR1/z', 'CDS/P/UKIDSS/LAS/DR9/K')
PIXELS = 500


def fetch_cutout(ra: float, dec: float, survey: str, fov: float, urls: tuple[str, ...] = HIPS2FITS):
    """North-up TAN cut-out ``fov`` arcmin across, centred on ra, dec (deg), from the first hips2fits that answers."""
    query = urlencode({'hips': survey, 'ra': ra, 'dec': dec, 'fov': fov / 60, 'width': PIXELS, 'height': PIXELS,
                       'projection': 'TAN', 'format': 'fits'})
    errors = []
    for url in urls:
        try:
            with urlopen(f"{url}?{query}", timeout=10) as response:
                return fits.open(BytesIO(response.read()))[0]
        except OSError as err:
            errors.append(f'{url}: {err}')
    raise OSError('No hips2fits service answered: ' + '; '.join(errors))


def plot_finder(hdu: fits.PrimaryHDU, name: str, survey: str, fov: float, path: Path) -> Path:
    """Inverted greyscale chart with a centre reticle and N/E arrows."""
    fig = Figure(figsize=(6, 6))
    ax = fig.add_subplot(projection=WCS(hdu.header))
    data = np.squeeze(hdu.data)
    vmin, vmax = np.nanpercentile(data, [1, 99.5])
    ax.imshow(data, origin='lower', cmap='gray_r', vmin=vmin, vmax=vmax)
    c, gap, arm = (PIXELS - 1) / 2, PIXELS / 30, PIXELS / 10  # 0-indexed image centre
    for sign in (-1, 1):  # reticle: four ticks leaving the target itself clear
        ax.plot([c + sign * gap, c + sign * (gap + arm)], [c, c], color='#DDAA33', lw=1.5)
        ax.plot([c, c], [c + sign * gap, c + sign * (gap + arm)], color='#DDAA33', lw=1.5)
    for dx, dy, label in ((0, 0.12, 'N'), (-0.12, 0, 'E')):  # hips2fits TAN is north up, east left
        ax.annotate(label, xy=(0.88, 0.08), xytext=(0.88 + dx, 0.08 + dy), xycoords='axes fraction',
                    color='#004488', ha='center', va='center',
                    arrowprops={'arrowstyle': '<-', 'color': '#004488'})
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_title(f'{name}   {survey}   {fov:g} arcmin')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    return path


def aladin_chart(ra: float, dec: float, survey: str, fov: float, path: Path, aladin: str) -> Path:
    """Chart drawn by Aladin (its own reticle, N/E arrows and scale) from a HiPS mirror; zoom snaps near ``fov``."""
    if survey not in MIRRORS:
        raise OSError(f'No HiPS mirror known for {survey}; use one of {list(MIRRORS)}')
    script = f"load {MIRRORS[survey]}\nsync\n{ra} {dec}\nzoom {fov}'\nsync\nsave {PIXELS}x{PIXELS} {path.name}\nquit\n"
    run = subprocess.run([aladin, '-nogui'], input=script, text=True, capture_output=True, timeout=300,
                         cwd=path.parent)  # relative file name: Aladin scripts do not quote paths
    if not path.is_file() or path.stat().st_size == 0:
        raise OSError(f'Aladin made no chart: {(run.stdout + run.stderr).strip()[-300:]}')
    return path


def make_finder(name: str, ra: float, dec: float, out_dir: Path, survey: str = SURVEYS[0], fov: float = 2.0,
                urls: tuple[str, ...] = HIPS2FITS) -> Path:
    """Fetch and plot one finder chart to ``out_dir/<name>_<fov>arcmin.png``; Aladin (if on PATH) when no
    hips2fits service answers."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f'{name}_{fov:g}arcmin.png'
    try:
        hdu = fetch_cutout(ra, dec, survey, fov, urls)
    except OSError as err:
        aladin = shutil.which('aladin')
        if aladin is None:
            raise OSError(f'{err}. Install Aladin (an "aladin" command on PATH) to use HiPS mirrors instead') from err
        return aladin_chart(ra, dec, survey, fov, path, aladin)
    return plot_finder(hdu, name, survey, fov, path)
