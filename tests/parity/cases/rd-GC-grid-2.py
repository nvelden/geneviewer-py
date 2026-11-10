"""Parity case: the chart of fixture rd-GC-grid-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").grid(
        height="120px", cluster=2
    )
