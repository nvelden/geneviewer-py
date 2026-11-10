"""Parity case: the chart of fixture rd-GC-clusterFooter-2."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group").cluster_footer(
        title="This is a footer",
        subtitle="Subtitle for the footer",
        show=True,
        cluster=1,
        x=6,
        y=-20,
        align="center",
        spacing=12,
        title_font={
            "fontSize": "12px",
            "fontWeight": "bold",
            "fontFamily": "sans-serif",
            "fill": "black",
            "cursor": "default",
        },
        subtitle_font={
            "fill": "grey",
            "fontSize": "10px",
            "fontStyle": "normal",
            "fontFamily": "sans-serif",
            "cursor": "default",
        },
    )
