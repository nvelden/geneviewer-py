"""Parity case: the chart of fixture rd-geneviewer-palette-1."""

from geneviewer import GC


def build(data):
    return GC(data, group="group", height="100px").genes(color_scheme="okabeIto")
