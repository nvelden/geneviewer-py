"""Parity case: the chart of fixture rd-GC-trackMouse-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="220px").track_mouse()
