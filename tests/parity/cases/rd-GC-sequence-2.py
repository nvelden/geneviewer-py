"""Parity case: the chart of fixture rd-GC-sequence-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="200px")
        .scale(hidden=True, scale_breaks=True)
        .sequence(
            start=None,
            end=None,
            sequence_style={"stroke": "blue", "strokeWidth": 1},
            marker_style={
                "stroke": "blue",
                "strokeWidth": 1,
                "gap": 3,
                "tiltAmount": 5,
            },
        )
        .legend(False)
    )
