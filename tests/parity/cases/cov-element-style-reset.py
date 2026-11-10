"""Parity case: the chart of fixture cov-element-style-reset."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="class", cluster="cluster", height="250px")
        .labels("name")
        .element_style("labels", row_id=1, fill="red")
        .element_style("labels", row_id=2, fill="blue")
        .element_style("labels", row_id=2, reset=True)
    )
