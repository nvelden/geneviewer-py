"""Parity case: the chart of fixture rd-GC-clusterLabel-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").cluster_label(
        title=data["cluster"].unique().tolist()
    )
