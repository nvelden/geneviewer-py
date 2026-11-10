"""Parity case: the chart of fixture rd-GC-genes-1."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").genes(
        group="group",
        show=True,
        marker="arrow",
        marker_size="medium",
        color_scheme=None,
        custom_colors=None,
        prevent_overlap=False,
        gene_overlap_spacing=40,
        cluster=1,
        x=1,
        y=50,
        stroke="black",
        stroke_width=1,
        arrowhead_width=None,
        arrowhead_height=None,
        arrow_height=None,
        marker_height=None,
    )
