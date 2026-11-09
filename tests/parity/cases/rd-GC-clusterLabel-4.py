"""Parity case: the chart of fixture rd-GC-clusterLabel-4."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").cluster_label(
        title=["Cluster 1", "Cluster 2"],
        width="100px",
        x=0,
        y=0,
        position="left",
        align="left",
        wrap_label=True,
        wrap_options={
            "dyAdjust": 0,
            "lineHeightEms": 1.05,
            "lineHeightSquishFactor": 1,
            "splitOnHyphen": True,
            "centreVertically": True,
        },
        font_size="12px",
        font_style="normal",
        font_weight="bold",
        font_family="sans-serif",
        cursor="default",
    )
