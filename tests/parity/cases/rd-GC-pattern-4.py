"""Parity case: the chart of fixture rd-GC-pattern-4."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="150px").pattern(
        patterns={"A": "checkerboard", "B": "rings", "C": "bricks"},
        color={"A": "black", "B": "black", "C": "white"},
        background={"A": "white", "B": "white", "C": "black"},
    )
