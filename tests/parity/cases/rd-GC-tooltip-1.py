"""Parity case: the chart of fixture rd-GC-tooltip-1."""

from geneviewer import GC


def build(data):
    return GC(data, group="group", height="200px").tooltip(
        formatter=" <b>Start:</b> {start}<br><b>end:</b> {end}"
    )
