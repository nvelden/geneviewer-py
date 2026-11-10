"""Example Shiny for Python app: three charts, events and ``zoom_to``.

Run it with ``shiny run examples/shiny/app.py``. It needs the ``shiny``
extra: ``pip install "geneviewer-py[shiny]"``.

Click, double-click or drag parts of the first chart, or move its zoom
slider, and the panel on the right shows each event input. The buttons move
the slider from the server, and the theme menu sends a new chart.

The second chart links the genes of two clusters by a column they share, and
the third draws transcripts; both use the package's example datasets.
"""

import json

import pandas as pd
from shiny import App, reactive, render, ui

from geneviewer import GC, datasets
from geneviewer.shiny import EVENTS, output_geneviewer, render_geneviewer, zoom_to

# The ophA gene cluster of Omphalotus olearius, from the docs site's
# getting started
genes = pd.DataFrame(
    {
        "name": ["ophB1", "ophC", "ophA", "ophD", "ophB2", "ophP", "ophE"],
        "start": [2522, 5286, 9536, 12616, 13183, 19346, 20170],
        "end": [4276, 4718, 10904, 11859, 15046, 16016, 21484],
        "class": [
            "Monooxygenase",
            "NTF2-like",
            "Methyltransferase",
            "O-acyltransferase",
            "Monooxygenase",
            "Prolyloligopeptidase",
            "F-box/RNHI-like",
        ],
    }
)

# Two ophA clusters, and the transcripts of BRCA1
clusters = datasets.load("ophA_clusters")
transcripts = datasets.load("BRCA1_splice_variants")

THEMES = ["default", "classic", "dark", "vintage", "infographic", "publication", "bw"]

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_select("theme", "Theme", THEMES),
        ui.input_action_button("first_half", "Zoom to the first half"),
        ui.input_action_button("everything", "Show everything"),
        width=260,
    ),
    ui.layout_columns(
        ui.card(
            ui.card_header("Genes and labels"),
            output_geneviewer("chart", height="260px"),
        ),
        ui.card(
            ui.card_header("Events"),
            ui.output_code("events"),
        ),
        col_widths=[8, 4],
    ),
    ui.layout_columns(
        ui.card(
            ui.card_header("Links between clusters"),
            output_geneviewer("links", height="300px"),
        ),
        ui.card(
            ui.card_header("Transcripts"),
            output_geneviewer("transcript", height="300px"),
        ),
        col_widths=[6, 6],
    ),
    title="geneviewer in Shiny for Python",
)


def server(input, output, session):
    """Draw the charts, move the zoom slider and show the events."""

    @render_geneviewer
    def chart():
        return (
            GC(genes, group="class", height="260px")
            .labels("name")
            .zoom()
            .events(
                highlight=True,
                dblclick=True,
                select=True,
                drag=True,
                zoom=True,
                rendered=True,
            )
            .theme(input.theme())
        )

    @render_geneviewer
    def links():
        return (
            GC(clusters, cluster="cluster", group="class", height="300px")
            .links(group="class")
            .labels("name")
            .theme(input.theme())
        )

    @render_geneviewer
    def transcript():
        return (
            GC(transcripts, height="300px")
            .transcript(transcript="ensembl_transcript_id", strand="strand")
            .theme(input.theme())
        )

    @reactive.effect
    @reactive.event(input.first_half)
    def _():
        zoom_to("chart", [0, 0.5])

    @reactive.effect
    @reactive.event(input.everything)
    def _():
        zoom_to("chart", None)

    # The last value of each event input, input.chart_<event>()
    last = reactive.value({})

    def keep(name):
        event = input[f"chart_{name}"]

        @reactive.effect
        @reactive.event(event)
        def _():
            last.set({**last.get(), name: event()})

    for name in EVENTS:
        keep(name)

    @render.code
    def events():
        if not last():
            return "No events yet"
        return "\n\n".join(
            f"input.chart_{name}():\n{json.dumps(value, indent=1)}"
            for name, value in last().items()
        )


app = App(app_ui, server)
