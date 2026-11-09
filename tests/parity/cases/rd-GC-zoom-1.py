"""Parity case: the chart of fixture rd-GC-zoom-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="class", height="300px")
        .labels("name")
        .zoom(range=[0.2, 0.6])
    )
