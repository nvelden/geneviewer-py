"""Parity case: the chart of fixture rd-GC-chart-2."""

from geneviewer import GC


def build(data):
    # R reads a GenBank file and sets cluster and strand; the data is the result
    return GC(data, cluster="cluster", strand="strand", height="200px")
