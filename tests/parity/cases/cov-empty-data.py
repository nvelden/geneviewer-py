"""Parity case: the chart of fixture cov-empty-data."""

from geneviewer import GC


def build(data):
    # ophA_clusters without rows: the fixture's data is its columns
    return (
        GC(data, cluster="cluster", group="class", height="200px")
        .title("No genes")
        .labels("name")
        .legend()
        .coordinates()
        .zoom()
    )
