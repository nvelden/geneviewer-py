"""Parity case: the chart of fixture rd-GC-theme-6."""

from geneviewer import GC


def build(data):
    chart = (
        GC(data, cluster="cluster", group="group", height="200px")
        .labels("name")
        .cluster_label()
    )
    return chart.theme("classic", labels={"fontSize": "14px", "fill": "grey"}).labels(
        "name", font_size="10px"
    )
