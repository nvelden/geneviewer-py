"""Parity case: the chart of fixture cov-toolbar."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="class", cluster="cluster", height="250px")
        .events(drag=True)
        .toolbar(
            show="always",
            buttons=["undo", "redo", "save"],
            font_family="serif",
            tooltip_color="#333333",
            tooltip_border="1px solid #b7b9be",
            tooltip_radius="4px",
            tooltip_shadow="1px 2px 10px rgba(0, 0, 0, 0.2)",
        )
    )
