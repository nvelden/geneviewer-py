"""Parity case: the chart of fixture rd-coronavirus-genomes-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="category", height="700px")
        .labels("gene")
        .cluster_label()
        .zoom(range=[0.75, 1])
    )
