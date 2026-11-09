"""Parity case: the chart of fixture rd-GC-color-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="100px").color(
        color_scheme="schemeCategory10"
    )
