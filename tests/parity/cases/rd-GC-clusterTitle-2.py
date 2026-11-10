"""Parity case: the chart of fixture rd-GC-clusterTitle-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", group="group", height="400px")
        .labels("name")
        .cluster_title(
            title="Cluster 1 Data",
            subtitle="Detailed View",
            show=True,
            cluster=1,
            x=0,
            y=5,
            align="center",
            spacing=20,
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
