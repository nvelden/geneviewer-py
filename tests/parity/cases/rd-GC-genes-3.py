"""Parity case: the chart of fixture rd-GC-genes-3."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").genes(
        group_style={"C": {"marker": "box", "stroke": "red"}}
    )
