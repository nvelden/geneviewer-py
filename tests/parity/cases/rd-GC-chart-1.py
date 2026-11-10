"""Parity case: the chart of fixture rd-GC-chart-1."""

from geneviewer import GC


def build(data):
    return GC(data, group="group", cluster="cluster", height="200px").labels("name")
