"""Parity case: the chart of fixture rd-GC-normalize-1."""

import pandas as pd

from geneviewer import GC


def build(data):
    genes = pd.DataFrame(
        {
            "start": [10, 90, 130, 10, 90, 130],
            "end": [40, 120, 160, 40, 120, 160],
            "name": ["Gene 1", "Gene 2", "Gene 3"] * 2,
            "group": ["A", "B", "C"] * 2,
            "cluster": [1, 1, 1, 2, 2, 2],
        }
    )
    return GC(genes, group="group", cluster="cluster", height="150px").normalize(
        cluster=2, preserve_gene_length=True, gap=None
    )
