"""Parity case: the chart of fixture rd-GC-title-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="400px")
        .labels("name")
        .title(title="Cluster 1 Data", subtitle="Detailed View", show=True)
    )
