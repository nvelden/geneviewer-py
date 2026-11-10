"""Parity case: the chart of fixture rd-GC-color-4."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="100px").color(
        color_scheme="schemeSet2", custom_colors={"C": "grey"}
    )
