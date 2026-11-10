"""Parity case: the chart of fixture rd-GC-grid-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").grid(
        margin={"left": "50px", "right": "0px"}
    )
