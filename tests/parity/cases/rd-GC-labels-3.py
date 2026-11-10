"""Parity case: the chart of fixture rd-GC-labels-3."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").labels(
        label="group",
        show=True,
        x=0,
        y=50,
        dy="-1.2em",
        dx="0em",
        rotate=0,
        adjust_labels=True,
        font_size="12px",
        font_style="italic",
        fill="black",
        font_family="sans-serif",
        text_anchor="middle",
        cursor="default",
    )
