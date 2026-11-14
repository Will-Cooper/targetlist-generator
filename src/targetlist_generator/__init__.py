"""Plan an observing night from a target table: visibility, airmass-ordered blocks with nearby telluric standards,
visibility plots, finder charts and telescope target lists."""
from importlib.metadata import version

from .pipeline import make_finders, make_targetlist

__version__ = version('targetlist-generator')
__all__ = ['make_finders', 'make_targetlist', '__version__']
