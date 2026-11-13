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


def chart_path(out_dir: Path, name: str, fov: float) -> Path:
    return out_dir / f'{name}_{fov:g}arcmin.png'


def aladin_charts(items: list[tuple[float, float, Path]], survey: str, fov: float, aladin: str) -> list[Path]:
    """Charts drawn by one Aladin session (its own reticle, N/E arrows and scale) from a HiPS mirror; zoom snaps
    near ``fov``. ``items`` are (ra, dec, path) in one folder; returns the paths actually written."""
    if survey not in MIRRORS:
        raise OSError(f'No HiPS mirror known for {survey}; use one of {list(MIRRORS)}')
    lines = [f'load {MIRRORS[survey]}', 'sync']
    for ra, dec, path in items:
        path.unlink(missing_ok=True)  # so a stale chart is never mistaken for a new one
        lines += [f'{ra} {dec}', f"zoom {fov}'", 'sync', f'save {PIXELS}x{PIXELS} {path.name}']
    try:  # relative file names: Aladin scripts do not quote paths
        subprocess.run([aladin, '-nogui'], input='\n'.join(lines + ['quit', '']), text=True, capture_output=True,
                       timeout=120 + 30 * len(items), cwd=items[0][2].parent)
    except subprocess.TimeoutExpired as err:
        raise OSError(f'Aladin timed out after {err.timeout:.0f} s') from err
    return [path for _, _, path in items if path.is_file() and path.stat().st_size > 0]


def make_finders(items: list[tuple[str, float, float]], out_dir: Path, survey: str = SURVEYS[0], fov: float = 2.0,
                 urls: tuple[str, ...] = HIPS2FITS, log=print) -> tuple[list[Path], list[str]]:
    """Charts for (name, ra, dec) items into ``out_dir``; returns (written paths, failed names).
    hips2fits first; from the first target no hips2fits service answers for, the rest go to one Aladin session
    (an "aladin" command on PATH) instead of retrying dead services per target."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for k, (name, ra, dec) in enumerate(items):
        try:
            hdu = fetch_cutout(ra, dec, survey, fov, urls)
        except OSError as err:
            rest = items[k:]
            log(f'hips2fits unavailable ({err})')
            aladin = shutil.which('aladin')
            if aladin is None:
                log('Install Aladin (an "aladin" command on PATH) to draw charts from HiPS mirrors instead')
                return paths, [n for n, _, _ in rest]
            log(f'Drawing {len(rest)} charts in one Aladin session from {MIRRORS.get(survey, survey)}')
            wanted = {n: chart_path(out_dir, n, fov) for n, _, _ in rest}
            try:
                made = aladin_charts([(r, d, wanted[n]) for n, r, d in rest], survey, fov, aladin)
            except OSError as aladin_err:
                log(f'Aladin failed: {aladin_err}')
                made = []
            return paths + made, [n for n, p in wanted.items() if p not in made]
        paths.append(plot_finder(hdu, name, survey, fov, chart_path(out_dir, name, fov)))
        log(f'Finder chart: {paths[-1].name}')
    return paths, []
