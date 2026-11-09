"""Parity case: the chart of fixture rd-GC-labels-4."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="200px")
        .labels(label="group")
        .labels(
            cluster=1,
            item_style=[
                {"index": 0, "fill": "red", "fontSize": "14px", "fontWeight": "bold"}
            ],
        )
    )
