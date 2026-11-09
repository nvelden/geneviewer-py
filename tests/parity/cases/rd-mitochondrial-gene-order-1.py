"""Parity case: the chart of fixture rd-mitochondrial-gene-order-1."""

from geneviewer import GC


def build(data):
    return (
        GC(
            data,
            start="plot_start",
            end="plot_end",
            cluster="cluster",
            group="block",
            strand="strand",
        )
        .labels("gene")
        .cluster_label()
    )
