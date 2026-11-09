"""Parity case: the chart of fixture rd-GC-cluster-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, height="150px")
        .cluster(prevent_gene_overlap="nested")
        .genes(
            style_by="geneTrack",
            group_style={"2": {"markerHeight": 8, "arrowheadHeight": 8}},
        )
        .labels("name")
    )
