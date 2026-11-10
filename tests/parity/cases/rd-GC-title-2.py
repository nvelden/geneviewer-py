"""Parity case: the chart of fixture rd-GC-title-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="500px")
        .labels("name")
        .title(
            title="Cluster 1 Data",
            subtitle="Detailed View",
            show=True,
            height="50px",
            cluster=1,
            x=0,
            y=25,
            align="center",
            spacing=20,
            style={"backgroundColor": "#0000"},
            title_font={
                "fontSize": "16px",
                "fontStyle": "normal",
                "fontWeight": "bold",
                "textDecoration": "normal",
                "fontFamily": "sans-serif",
                "cursor": "default",
                "fill": "black",
            },
            subtitle_font={
                "fontSize": "14px",
                "fontStyle": "normal",
                "fontWeight": "bold",
                "textDecoration": "normal",
                "fontFamily": "sans-serif",
                "cursor": "default",
                "fill": "black",
            },
        )
    )
