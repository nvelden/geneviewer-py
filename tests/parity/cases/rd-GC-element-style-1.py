"""Parity case: the chart of fixture rd-GC-element-style-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data)
        .labels("name")
        .cluster_title("Cluster", subtitle="Subtitle")
        .element_style("genes", row_id=1, fill="red", marker="box")
        .element_style("labels", row_id=1, text="Edited", font_size="16px")
        .element_style(
            "clusterTitle", cluster=1, part="subtitle", font_style="italic", fill="grey"
        )
    )
