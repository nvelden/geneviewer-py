"""Parity case: the chart of fixture cov-links-low-color."""

from geneviewer import GC


def build(data):
    return GC(data, group="class", cluster="cluster", height="250px").links(
        "class",
        measure="identity",
        low_color="#dddddd",
        normal_color="#1f77b4",
        opacity=0.6,
    )
