"""Parity case: the chart of fixture vig-geneviewer-25."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="class", height="150px")
        .align(id_column="class", id="NTF2-like", align="left")
        .legend(True)
    )
