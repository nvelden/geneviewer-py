"""Parity case: the chart of fixture rd-GC-links-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", height="200px")
        .links(group="group", value1="A", value2="B", measure="identity", label=False)
        .labels(label="group")
    )
