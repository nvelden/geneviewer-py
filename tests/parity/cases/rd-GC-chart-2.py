"""Parity case: the chart of fixture rd-GC-chart-2."""

from pathlib import Path

from geneviewer import GC

# The GenBank file of the R package's inst/extdata, as the reader tests use it
GBK = Path(__file__).resolve().parents[2] / "readers" / "inputs" / "BGC0000001.gbk"


def build(data):
    # R reads the GenBank file itself; so does GC()
    return GC(str(GBK), height="200px")
