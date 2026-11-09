"""Parity case: the chart of fixture rd-GC-coordinates-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="name", height="200px").coordinates()
