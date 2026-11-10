"""Parity case: the chart of fixture rd-GC-annotation-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="220px")
        .annotation(
            type="textMarker", cluster=1, position=24, text="Gene 1", arrow_size=8
        )
        .annotation(type="text", text="feature 1", x=91, y=71)
        .annotation(type="symbol", symbol="triangle", x=95, y=64, size=10, rotation=180)
        .annotation(type="terminator", x=81)
        .annotation(type="promoter", x=49)
        .track_mouse()
    )
