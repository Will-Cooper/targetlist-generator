"""Cut target lists to a night's visibility, pair telluric standards, write IRTF/ING lists and finder charts."""
from importlib.metadata import version

from .pipeline import make_finders, make_targetlist

__version__ = version('targetlist-generator')
__all__ = ['make_finders', 'make_targetlist', '__version__']
