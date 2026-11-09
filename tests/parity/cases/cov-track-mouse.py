"""Parity case: the chart of fixture cov-track-mouse."""

from geneviewer import GC


def build(data):
    return GC(data, group="class", cluster="cluster", height="250px").track_mouse(
        show=True, cluster=1, background_color="black", color="white"
    )
