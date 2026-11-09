"""Parity case: the chart of fixture rd-GC-cluster-1."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="group", height="150px")
        .cluster(separate_strands=True, strand_spacing=0)
        .legend(False)
    )
