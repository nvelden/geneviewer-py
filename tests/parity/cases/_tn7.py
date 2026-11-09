"""The Tn7-CRISPR vignette chart, step by step, as in its R code.

Fixtures vig-Tn7-CRISPR-1 to -6 are the chart after each step; the
vig-Examples-15 fixture makes the same chart in one call.
"""

from geneviewer import GC

COLORS = {
    "insertion site": "#F3C79B",
    "tnsA": "#D62728",
    "tnsB": "#AEB6F0",
    "tnsC": "#D7F2F2",
    "tniQ / tnsD": "#2CA02C",
    "tnsE": "#C7E02C",
    "cas8f / cas5f": "#8C1515",
    "cas7f": "#F2A5C4",
    "cas6f": "#7B3FD4",
    "cas5": "#B5651D",
    "cas, other": "#C9A0DC",
    "CRISPR array": "#555555",
    "right end": "#9FD4DE",
    "left end": "#F28E2B",
    "other": "#E8E8E8",
}

_END_BAR = {"marker": "box", "markerHeight": 38, "arrowheadHeight": 38}
_GENE_BOX = {"marker": "box", "markerHeight": 20, "arrowheadHeight": 20}
SHAPES = {
    "R": _END_BAR,
    "L": _END_BAR,
    "CRISPR": {**_GENE_BOX, "stroke": "black", "strokeWidth": 1},
    "tnsA (fused)": _GENE_BOX,
}
PATTERN = {
    "patterns": {"right end": "diagonal", "left end": "diagonal"},
    "background": None,
    "size": 4,
}

_WHITE = {"fill": "white"}
LABEL_STYLES = {
    "tnsA": _WHITE,
    "cas8f": _WHITE,
    "cas6f": _WHITE,
    "cas5": _WHITE,
    "SRP-RNA": {"y": 55, "dy": "0px", "textAnchor": "end", "x": -4},
}
_GENE_NAMES = [
    "tnsA",
    "tnsB",
    "tnsC",
    "tniQ / tnsD",
    "tnsE",
    "cas8f / cas5f",
    "cas7f",
    "cas6f",
    "cas5",
]
LEGEND_STYLES = {name: {"fontStyle": "italic"} for name in _GENE_NAMES}

TITLE_ARGS = ("Tn7-like transposons that carry a minimal CRISPR-Cas system",)
TITLE = {
    "subtitle": "Redrawn from Fig. 3 of Peters et al. (2017), "
    "Proc Natl Acad Sci USA 114:E7358-E7366",
    "height": "70px",
    "title_font": {"fontSize": "19px", "fontWeight": "bold"},
    "subtitle_font": {"fontSize": "14px", "fontStyle": "italic"},
}
TOOLTIP = "<b>{gene}</b> {product}<br>{start} - {end}"


def cluster_label(elements):
    """Return the arguments of the cluster labels."""
    return {
        "title": elements["cluster"].tolist(),
        "subtitle": elements["strain_short"].tolist(),
        "width": "190px",
        "font_size": "12px",
        "subtitle_font": {"fontSize": "10px", "fontStyle": "italic"},
        "align": "right",
        "x": -10,
    }


def step(data, n):
    """Return the chart of fixture vig-Tn7-CRISPR-<n>."""
    chart = (
        GC(data, cluster="cluster", group="category", strand="strand", height="900px")
        .color(custom_colors=COLORS)
        .scale(axis_type="range", hidden=True)
        .legend(position="bottom", align="center", order=list(COLORS))
    )
    if n >= 2:
        chart = chart.scale(
            axis_type="range",
            hidden=True,
            scale_breaks=True,
            scale_break_padding="500bp",
        )
    if n >= 3:
        chart = chart.align(id_column="gene", id="R", align="left")
    if n >= 4:
        chart = chart.genes(
            marker_size="small", marker_height=20, style_by="gene", group_style=SHAPES
        ).pattern(**PATTERN)
    if n >= 5:
        chart = chart.labels(
            label="label",
            font_size="9px",
            y=50,
            dy="14.5px",
            adjust_labels=False,
            style_by="label",
            group_style=LABEL_STYLES,
        )
    if n >= 6:
        elements = data[["cluster", "strain_short"]].drop_duplicates()
        chart = (
            chart.cluster_label(**cluster_label(elements))
            .title(*TITLE_ARGS, **TITLE)
            .tooltip(formatter=TOOLTIP)
            .legend(legend_text_options={"fontSize": "12px"}, group_style=LEGEND_STYLES)
            .theme("publication")
        )
    return chart
