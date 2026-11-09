"""Parity case: the chart of fixture rd-GC-links-3."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="220px")
        .links(
            group="group",
            measure="identity",
            label=False,
            color_bar_options={"position": "legend"},
        )
        .legend(position="top")
    )
