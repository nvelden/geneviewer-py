"""Parity case: the chart of fixture cov-corner-radius."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="class", cluster="cluster", height="250px")
        .genes(marker="rbox", corner_radius=6)
        .legend(corner_radius=4)
    )
