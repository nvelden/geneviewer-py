"""Parity case: the chart of fixture rd-GC-labels-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").labels(
        group_style={"A": {"fontWeight": "bold"}}
    )
