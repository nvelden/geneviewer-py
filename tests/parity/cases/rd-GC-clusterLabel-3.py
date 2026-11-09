"""Parity case: the chart of fixture rd-GC-clusterLabel-3."""

from geneviewer import GC


def build(data):
    return GC(data, cluster="cluster", group="group", height="200px").cluster_label(
        title=["pQEB1", "pQEB1_inv1"],
        subtitle=["E. cloacae", "K. oxytoca"],
        subtitle_font={"fontStyle": "italic"},
        spacing=14,
    )
