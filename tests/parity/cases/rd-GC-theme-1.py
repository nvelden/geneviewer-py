"""Parity case: the chart of fixture rd-GC-theme-1."""

from geneviewer import GC


def build(data):
    chart = (
        GC(data, cluster="cluster", group="group", height="200px")
        .labels("name")
        .cluster_label()
    )
    return chart.theme("dark")
