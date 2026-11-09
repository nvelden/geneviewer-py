"""Parity case: the chart of fixture rd-GC-scale-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="200px")
        .scale(
            start=1,
            end=2600,
            padding=2,
            hidden=False,
            reverse=False,
            axis_position="bottom",
            scale_breaks=True,
            scale_break_threshold=20,
            scale_break_padding=1,
            ticks_count=10,
            ticks_format=",.0f",
            y=None,
            tick_style={"stroke": "grey", "strokeWidth": 1, "lineLength": 6},
            text_style={
                "fill": "black",
                "fontSize": "10px",
                "fontFamily": "Arial",
                "cursor": "default",
            },
            line_style={"stroke": "grey", "strokeWidth": 1},
        )
        .legend(False)
    )
