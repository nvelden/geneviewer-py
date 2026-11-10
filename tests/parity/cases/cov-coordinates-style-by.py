"""Parity case: the chart of fixture cov-coordinates-style-by."""

from geneviewer import GC


def build(data):
    return GC(data, group="class", cluster="cluster", height="250px").coordinates(
        style_by="name",
        group_style={"ophA": {"fill": "blue", "fontWeight": "bold"}},
    )
