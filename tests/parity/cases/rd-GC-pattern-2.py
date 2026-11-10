"""Parity case: the chart of fixture rd-GC-pattern-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="150px").pattern(
        patterns={"A": "dots", "B": "crosshatch", "C": "none"},
        background=["white", "grey80"],
    )
