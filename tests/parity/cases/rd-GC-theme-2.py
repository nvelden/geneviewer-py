"""Parity case: the chart of fixture rd-GC-theme-2."""

from geneviewer import GC


def build(data):
    chart = (
        GC(data, cluster="cluster", group="group", height="200px")
        .labels("name")
        .cluster_label()
    )
    return chart.theme(
        "classic", font_family="Georgia", palette=["#1b9e77", "#d95f02", "#7570b3"]
    )
