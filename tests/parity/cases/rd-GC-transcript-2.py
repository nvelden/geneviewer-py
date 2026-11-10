"""Parity case: the chart of fixture rd-GC-transcript-2."""

from _inputs import TRANSCRIPTS

from geneviewer import GC


def build(data):
    return GC(TRANSCRIPTS, start="start", end="end", height="200px").transcript(
        transcript="transcript",
        type="type",
        selection=2,
        item_style_exons=[{"index": 0, "fill": "red"}],
    )
