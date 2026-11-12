"""Parity case: the chart of fixture cov-hide-shown."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", height="250px")
        .labels("name")
        .labels(show=False, cluster=2)
        .coordinates()
        .coordinates(show=False, cluster=1)
        .cluster_title("Title")
        .cluster_title(show=False, cluster=2)
    )
