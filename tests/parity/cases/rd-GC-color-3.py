"""Parity case: the chart of fixture rd-GC-color-3."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="100px").color(
        custom_colors={"A": "yellow", "B": "pink", "C": "purple"}
    )
