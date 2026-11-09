"""Parity case: the chart of fixture rd-GC-links-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", height="200px")
        .labels(label="group")
        .links(
            group="group",
            value1=["A", "A"],
            value2=["B", "A"],
            label=False,
            cluster=[1, 2],
        )
    )
