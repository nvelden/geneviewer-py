"""Parity case: the chart of fixture cov-zoom-no-overview."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="class", height="300px").zoom(
        range=[0.3, 0.7], overview=False
    )
