"""Gene cluster charts for Python.

The package builds chart specs for the geneviewer JavaScript library and
shows them in Shiny for Python. The import name is ``geneviewer``; the
distribution on PyPI is ``geneviewer-py``.
"""

from importlib.metadata import version

from .align import mummer_alignment, protein_blast
from .chart import GC
from .processing import get_introns
from .readers import (
    gbk_features_to_df,
    genbank_to_fasta,
    read_bed,
    read_fasta,
    read_gbk,
    read_gff,
)
from .theme import geneviewer_palette, geneviewer_theme

__version__ = version("geneviewer-py")

__all__ = [
    "GC",
    "__version__",
    "gbk_features_to_df",
    "genbank_to_fasta",
    "geneviewer_palette",
    "geneviewer_theme",
    "get_introns",
    "mummer_alignment",
    "protein_blast",
    "read_bed",
    "read_fasta",
    "read_gbk",
    "read_gff",
]
