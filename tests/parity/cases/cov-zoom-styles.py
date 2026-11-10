"""Parity case: the chart of fixture cov-zoom-styles."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="class", height="300px").zoom(
        show=True,
        range=[0.1, 0.5],
        min=0.05,
        overview=True,
        color="grey",
        style={"backgroundColor": "#eeeeee"},
        window_style={"border": "1px solid red"},
    )
