"""Parity case: the chart of fixture cov-normalize."""

from _inputs import before_normalize

from geneviewer import GC


def build(data):
    # ophA_clusters, rebuilt from the result until the datasets ship (M7.8)
    genes = before_normalize(data)
    return GC(genes, group="class", cluster="cluster", height="250px").normalize(
        normalize=True, preserve_gene_length=False
    )
