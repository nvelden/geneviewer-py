"""Parity case: the chart of fixture rd-GC-scaleBar-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group")
        .scale_bar(cluster=1, title="1 kb", scale_bar_unit=1000)
        .scale_bar(cluster=2, title="2 kb", scale_bar_unit=2000)
    )
