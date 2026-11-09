"""Gene cluster charts for Python.

The package builds chart specs for the geneviewer JavaScript library and
shows them in Shiny for Python. The import name is ``geneviewer``; the
distribution on PyPI is ``geneviewer-py``.
"""

from importlib.metadata import version

__version__ = version("geneviewer-py")

__all__ = ["__version__"]
