"""Parity case: the chart of fixture rd-GC-links-4."""

from geneviewer import GC


def build(data):
    return (
        GC(data, cluster="cluster", height="200px")
        .links(
            group="group",
            data=None,
            curve=True,
            measure="identity",
            show_links=True,
            label=True,
            normal_color="#1f77b4",
            inverted_color="#d62728",
            use_group_colors=False,
            color_bar=True,
            color_bar_options={
                "x": 0,
                "y": 24,
                "width": 10,
                "height": 60,
                "labelOptions": {"fontSize": 8, "xOffset": 2, "yOffset": 0},
                "titleOptions": {"fontSize": 10, "xOffset": 2, "yOffset": 0},
                "barOptions": {"stroke": "#000", "strokeWidth": 0.5, "opacity": 1},
            },
            link_width=1,
            link_style={"stroke": "black", "strokeWidth": 0.5, "fillOpacity": 0.4},
            label_style={"fontSize": "8px"},
        )
        .labels(label="group", cluster=1)
        .cluster_label()
    )
