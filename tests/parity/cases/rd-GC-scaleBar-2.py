"""Parity case: the chart of fixture rd-GC-scaleBar-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="400px").scale_bar(
        title="1kb",
        scale_bar_unit=1000,
        cluster=None,
        x=0,
        y=50,
        label_style={
            "labelPosition": "left",
            "fontSize": "10px",
            "fontFamily": "sans-serif",
            "fill": "red",
            "cursor": "default",
        },
        text_padding=-2,
        scale_bar_line_style={"stroke": "black", "strokeWidth": 1},
        scale_bar_tick_style={"stroke": "black", "strokeWidth": 1},
    )
