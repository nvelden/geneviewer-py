"""Parity case: the chart of fixture cov-range-annotations."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="name", height="250px")
        .annotation(cluster=1, type="rectangle", position=[[[1100, 60], [1300, 40]]])
        .scale(axis_type="range")
        .annotation(cluster=1, type="text", text="x", x=1450)
        .annotation(
            cluster=2,
            type="rectangle",
            position=[[[250, 60], [450, 40]], [[650, 60], [850, 40]]],
        )
    )
