"""Parity case: the chart of fixture vig-Examples-15."""

import _tn7 as tn7

from geneviewer import GC


def build(data):
    elements = data[["cluster", "strain_short"]].drop_duplicates()
    return (
        GC(data, cluster="cluster", group="category", strand="strand", height="900px")
        .genes(
            marker_size="small",
            marker_height=20,
            style_by="gene",
            group_style=tn7.SHAPES,
        )
        .labels(
            label="label",
            font_size="9px",
            y=50,
            dy="14.5px",
            adjust_labels=False,
            style_by="label",
            group_style=tn7.LABEL_STYLES,
        )
        .scale(
            axis_type="range",
            hidden=True,
            scale_breaks=True,
            scale_break_padding="500bp",
        )
        .align("gene", "R", align="left")
        .cluster_label(**tn7.cluster_label(elements))
        .title(*tn7.TITLE_ARGS, **tn7.TITLE)
        .tooltip(formatter=tn7.TOOLTIP)
        .color(custom_colors=tn7.COLORS)
        .pattern(**tn7.PATTERN)
        .legend(
            position="bottom",
            align="center",
            order=list(tn7.COLORS),
            legend_text_options={"fontSize": "12px"},
            group_style=tn7.LEGEND_STYLES,
        )
        .theme("publication")
    )
