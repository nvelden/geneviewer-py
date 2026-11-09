"""Parity case: the chart of fixture cov-events."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="class", cluster="cluster", height="250px")
        .legend(True)
        .events(
            click=False, dblclick=True, legend_click=False, zoom=True, highlight=True
        )
    )
