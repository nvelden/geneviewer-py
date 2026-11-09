"""Parity case: the chart of fixture cov-tooltip-style."""

from geneviewer import GC


def build(data):
    return GC(data, group="class", cluster="cluster", height="250px").tooltip(
        formatter="<b>{name}</b>", background_color="#ffffcc"
    )
