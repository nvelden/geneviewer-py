"""Gene cluster charts for Python.

The package builds chart specs for the geneviewer JavaScript library and
shows them in Shiny for Python. The import name is ``geneviewer``; the
distribution on PyPI is ``geneviewer-py``.
"""

from importlib.metadata import version

from .chart import GC
from .theme import geneviewer_palette, geneviewer_theme

__version__ = version("geneviewer-py")

__all__ = ["GC", "__version__", "geneviewer_palette", "geneviewer_theme"]
