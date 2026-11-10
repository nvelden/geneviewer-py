"""Parity case: the chart of fixture cov-chart-without-extras."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="class", cluster="cluster", height="250px", save_button=False)
        .labels("name")
        .cluster_label(show=False)
        .scale_bar(show=False)
        .tooltip(show=False)
    )
