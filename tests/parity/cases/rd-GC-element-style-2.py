"""Parity case: the chart of fixture rd-GC-element-style-2."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="family")
        .element_style("legend", group="A", part="label", text="Family A")
        .element_style("legend", group="A", part="marker", fill="red")
    )
