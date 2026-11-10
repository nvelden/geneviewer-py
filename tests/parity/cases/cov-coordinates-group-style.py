"""Parity case: the chart of fixture cov-coordinates-group-style."""

from geneviewer import GC


def build(data):
    return GC(data, group="class", cluster="cluster", height="250px").coordinates(
        group_style={
            "Monooxygenase": {"fill": "red", "tick": {"stroke": "red"}},
            "NTF2-like": {"show": False},
        }
    )
