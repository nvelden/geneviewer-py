"""Parity case: the chart of fixture rd-GC-legend-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").legend(
        position="top",
        orientation="horizontal",
        align="center",
        x=0,
        y=0,
        adjust_height=True,
        background_color="#0000",
        order=[],
        style={"backgroundColor": "#0000"},
        legend_options={"cursor": "pointer", "colorScheme": None, "customColors": None},
        legend_text_options={
            "cursor": "pointer",
            "textAnchor": "start",
            "dy": ".35em",
            "fontSize": "12px",
            "fontFamily": "sans-serif",
            "fill": "black",
        },
    )
