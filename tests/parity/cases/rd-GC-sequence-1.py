"""Parity case: the chart of fixture rd-GC-sequence-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="200px")
        .labels("name")
        .sequence(show=True, y=50, cluster=None)
    )
