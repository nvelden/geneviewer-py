"""Parity case: the chart of fixture rd-GC-clusterFooter-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group").cluster_footer(
        title="Cluster Footer", subtitle="Cluster subtitle"
    )
