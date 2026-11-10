"""Parity case: the chart of fixture rd-GC-transcript-1."""

from _inputs import TRANSCRIPTS

from geneviewer import GC

_SIZE = {
    "arrowheadWidth": None,
    "arrowheadHeight": None,
    "markerHeight": None,
    "cornerRadius": None,
}


def build(data):
    return GC(TRANSCRIPTS, start="start", end="end", height="200px").transcript(
        transcript="transcript",
        strand="strand",
        type="type",
        group=None,
        show=True,
        selection=None,
        color_scheme=None,
        custom_colors=None,
        style_exons={
            "show": True,
            "strokeWidth": 0,
            "cursor": "default",
            "marker": "box",
            "markerSize": "medium",
            **_SIZE,
        },
        style_introns={
            "show": True,
            "strokeWidth": 1,
            "fill": "none",
            "cursor": "default",
            "marker": "intron",
            "markerSize": "medium",
            **_SIZE,
        },
        style_utrs={
            "show": True,
            "fontSize": "10px",
            "fontStyle": "normal",
            "fontFamily": "sans-serif",
            "cursor": "default",
            "color": "black",
            "fill": "#FFF",
            "strokeWidth": 1,
            "marker": "box",
            "markerSize": "medium",
            **_SIZE,
        },
        label_options={
            "show": True,
            "xOffset": 2,
            "yOffset": 0,
            "fontSize": "12px",
            "fontStyle": "normal",
            "fontWeight": "normal",
            "fontFamily": "sans-serif",
            "cursor": "default",
            "color": "black",
        },
        item_style_exons=[],
        item_style_introns=[],
        item_style_utrs=[],
    )
