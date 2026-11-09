"""Parity case: the chart of fixture rd-GC-coordinates-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="name", height="200px")
        .coordinates()
        .coordinates(
            cluster=2,
            show=True,
            tick_values_top=[130, 170, 210, 240],
            tick_values_bottom=[160, 200],
            ticks_format=",.0f",
            rotate=-45,
            y_position_top=55,
            y_position_bottom=45,
            overlap_threshold=20,
            tick_style={"stroke": "black", "strokeWidth": 1, "lineLength": 6},
            text_style={
                "fill": "black",
                "fontSize": "12px",
                "fontFamily": "Arial",
                "cursor": "default",
            },
        )
    )
