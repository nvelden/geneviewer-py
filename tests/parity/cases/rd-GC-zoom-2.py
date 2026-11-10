"""Parity case: the chart of fixture rd-GC-zoom-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="class", height="300px").zoom(
        mode="absolute",
        update="release",
        height=30,
        handle_style={"backgroundColor": "#d9d9d9"},
    )
