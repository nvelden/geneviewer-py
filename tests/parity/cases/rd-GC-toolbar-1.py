"""Parity case: the chart of fixture rd-GC-toolbar-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, height="150px")
        .labels("name")
        .events(select=True, drag=True)
        .toolbar(
            show="always",
            position="topleft",
            size=16,
            color="#9e9e9e",
            hover_color="black",
            icons={"reset": "arrows-to-circle", "save": "camera"},
        )
    )
