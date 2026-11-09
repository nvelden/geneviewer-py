import asyncio
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

pytest.importorskip("shiny")

from htmltools import TagList  # noqa: E402
from shiny import module  # noqa: E402
from shiny.module import ResolvedId  # noqa: E402

from geneviewer import GC  # noqa: E402
from geneviewer.shiny import (  # noqa: E402
    output_geneviewer,
    render_geneviewer,
    zoom_to,
)

GENES = pd.DataFrame(
    {"start": [10, 90], "end": [40, 120], "name": ["a", "b"], "group": ["A", "B"]}
)


def fake_session(ns=""):
    sent = []
    return SimpleNamespace(ns=ResolvedId(ns), _send_message_sync=sent.append), sent


# output_geneviewer ----------------------------------------------------------


def test_output_is_a_div_with_the_binding_class():
    tag = output_geneviewer("chart")
    assert tag.name == "div"
    assert tag.attrs["id"] == "chart"
    assert tag.attrs["class"] == "geneviewer-output"
    assert "width:100%;" in tag.attrs["style"]
    assert "height:400px;" in tag.attrs["style"]


def test_output_size_in_pixels_or_css():
    style = output_geneviewer("chart", width=500, height="50vh").attrs["style"]
    assert "width:500px;" in style
    assert "height:50vh;" in style


def test_output_brings_the_library_and_the_binding():
    (dep,) = TagList(output_geneviewer("chart")).get_dependencies()
    assert dep.name == "geneviewer"
    scripts = [script["src"] for script in dep.script]
    assert scripts == ["geneviewer.min.js", "geneviewer-shiny.js"]
    assert [sheet["href"] for sheet in dep.stylesheet] == ["geneviewer.css"]
    folder = dep.source["subdir"]
    for file in (*scripts, "geneviewer.css"):
        assert (Path(folder) / file).is_file()


def test_output_id_in_a_module():
    @module.ui
    def chart_ui():
        return output_geneviewer("chart")

    assert chart_ui("genes").attrs["id"] == "genes-chart"


# render_geneviewer ----------------------------------------------------------


def test_render_sends_the_spec():
    chart = GC(GENES, group="group").labels("name")
    renderer = render_geneviewer(lambda: chart)
    value = asyncio.run(renderer.transform(chart))
    assert value == {"spec": chart.to_spec()}


def test_render_needs_a_chart():
    renderer = render_geneviewer(lambda: GENES)
    with pytest.raises(TypeError, match="chart made with GC"):
        asyncio.run(renderer.transform(GENES))


# zoom_to --------------------------------------------------------------------


def test_zoom_to_sends_a_sorted_range():
    session, sent = fake_session()
    zoom_to("chart", [0.5, 0.25], session=session)
    assert sent == [
        {"custom": {"geneviewer-zoom": {"id": "chart", "range": [0.25, 0.5]}}}
    ]


def test_zoom_to_everything_and_ranges_from_inputs():
    session, sent = fake_session()
    zoom_to("chart", None, session=session)
    zoom_to("chart", (12000, 18000), session=session)
    message = [m["custom"]["geneviewer-zoom"]["range"] for m in sent]
    assert message == [None, [12000.0, 18000.0]]


def test_zoom_to_in_a_module():
    session, sent = fake_session("genes")
    zoom_to("chart", [0, 1], session=session)
    assert sent[0]["custom"]["geneviewer-zoom"]["id"] == "genes-chart"


@pytest.mark.parametrize(
    "window",
    [[1], [1, 2, 3], ["0", "1"], [True, False], [0, float("nan")], [0, float("inf")]],
)
def test_zoom_to_checks_the_range(window):
    session, sent = fake_session()
    with pytest.raises(ValueError, match="'range' must be None or two numbers"):
        zoom_to("chart", window, session=session)
    assert sent == []


def test_zoom_to_checks_the_id():
    session, _ = fake_session()
    with pytest.raises(ValueError, match="output id"):
        zoom_to(1, [0, 1], session=session)


def test_zoom_to_needs_a_session():
    with pytest.raises(RuntimeError, match="works in a Shiny session"):
        zoom_to("chart", [0, 1])
