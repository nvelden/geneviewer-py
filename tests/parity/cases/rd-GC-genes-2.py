"""Parity case: the chart of fixture rd-GC-genes-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").genes(
        cluster=2, item_style=[{"index": 2, "fill": "red", "stroke": "blue"}]
    )
