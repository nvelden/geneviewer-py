"""The Jupyter widget: its spec, events and zoom_to, without a browser."""

import subprocess
import sys

import pandas as pd
import pytest

pytest.importorskip("anywidget")

from geneviewer import GC, datasets  # noqa: E402
from geneviewer.jupyter import EVENTS, GCWidget, _module  # noqa: E402
from geneviewer.shiny import EVENTS as SHINY_EVENTS  # noqa: E402


@pytest.fixture
def chart():
    genes = datasets.load("ophA_clusters")
    return GC(genes, cluster="cluster", group="class").zoom()


def test_widget_holds_the_spec(chart):
    widget = chart.widget()
    assert isinstance(widget, GCWidget)
    assert widget.spec == chart.to_spec()
    assert widget.chart is chart
    assert widget.trait_metadata("spec", "sync") is True


def test_setting_the_chart_changes_the_spec(chart):
    widget = chart.widget()
    other = chart.title("Another")
    widget.chart = other
    assert widget.chart is other
    assert widget.spec["title"] == {"title": "Another"}
    with pytest.raises(TypeError, match="made with GC"):
        widget.chart = {"data": []}
    with pytest.raises(TypeError, match="made with GC"):
        GCWidget(pd.DataFrame())


def test_events_are_those_of_shiny():
    assert EVENTS == SHINY_EVENTS


def test_an_event_sets_its_attribute_and_calls_back(chart):
    widget = chart.widget()
    assert all(getattr(widget, name) is None for name in EVENTS)
    seen = []
    widget.on("select", seen.append)
    gene = {"component": "genes", "data": {"name": "ophA"}}
    widget._message(widget, {"type": "event", "name": "select", "detail": gene}, [])
    assert widget.select == gene
    # The same event again still calls back
    widget._message(widget, {"type": "event", "name": "select", "detail": gene}, [])
    assert seen == [gene, gene]
    widget.on("select", seen.append, remove=True)
    widget._message(widget, {"type": "event", "name": "select", "detail": None}, [])
    assert widget.select is None
    assert len(seen) == 2
    # Other messages are left alone
    widget._message(widget, {"type": "event", "name": "nope", "detail": 1}, [])
    widget._message(widget, "text", [])
    with pytest.raises(ValueError, match="There is no event 'nope'"):
        widget.on("nope", print)


def test_zoom_to_sends_the_window(chart, monkeypatch):
    widget = chart.widget()
    sent = []
    monkeypatch.setattr(widget, "send", sent.append)
    widget.zoom_to([0.5, 0])
    widget.zoom_to(None)
    assert sent == [
        {"type": "zoom", "range": [0.0, 0.5]},
        {"type": "zoom", "range": None},
    ]
    with pytest.raises(ValueError, match="two numbers"):
        widget.zoom_to([1])


def test_a_chart_shows_as_a_widget_in_a_notebook(chart):
    bundle = chart._repr_mimebundle_()
    data = bundle[0] if isinstance(bundle, tuple) else bundle
    assert "application/vnd.jupyter.widget-view+json" in data


def test_module_is_the_library_and_the_widget_code():
    module = _module()
    assert "__GENEVIEWER_NEW_CHART__" not in module
    assert module.count("export default") == 1
    assert module.rstrip().endswith("export default { render: __geneviewerRender };")
    assert "sourceMappingURL" not in module


def test_jupyter_is_optional():
    # geneviewer.jupyter loads anywidget; the package itself does not, and a
    # chart without it shows as text
    code = (
        "import sys, geneviewer; assert 'anywidget' not in sys.modules; "
        "sys.modules['anywidget'] = None; "
        "import pandas as pd; "
        "chart = geneviewer.GC(pd.DataFrame({'start': [1], 'end': [2]})); "
        "assert chart._repr_mimebundle_() is None"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
