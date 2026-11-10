"""Parity case: the chart of fixture rd-GC-align-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="group", cluster="cluster", height="150px")
        .align("group", "A", align="left")
        .legend(False)
    )
