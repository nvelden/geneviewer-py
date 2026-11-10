"""Parity case: the chart of fixture rd-GC-pattern-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="150px").pattern()
