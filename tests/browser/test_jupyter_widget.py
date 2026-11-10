"""The widget's JavaScript in a browser, with a stand-in for the notebook's model."""

import pytest

pytest.importorskip("anywidget")

from playwright.sync_api import Page  # noqa: E402

from geneviewer import GC, datasets  # noqa: E402
from geneviewer.jupyter import _module  # noqa: E402

pytestmark = pytest.mark.browser

# What anywidget gives render(): a model with get, on, off and send. The page
# keeps what the widget sends, and can change the spec and send a message
# ruff: noqa: E501
SETUP = """async ([source, spec]) => {
  const module = await import(URL.createObjectURL(new Blob([source], { type: "text/javascript" })));
  const listeners = {};
  window.sent = [];
  window.state = { spec };
  window.model = {
    get: (name) => window.state[name],
    on: (name, callback) => { (listeners[name] = listeners[name] || []).push(callback); },
    off: (name, callback) => { listeners[name] = (listeners[name] || []).filter((c) => c !== callback); },
    send: (content) => window.sent.push(content),
    emit: (name, ...args) => (listeners[name] || []).forEach((callback) => callback(...args)),
    count: (name) => (listeners[name] || []).length,
  };
  const el = document.createElement("div");
  el.id = "output";
  document.body.appendChild(el);
  window.cleanup = module.default.render({ model: window.model, el });
}"""


@pytest.fixture
def widget_page(page: Page) -> Page:
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    genes = datasets.load("ophA_clusters")
    chart = (
        GC(genes, cluster="cluster", group="class", height="300px")
        .labels("name")
        .zoom()
        .events(select=True, zoom=True, rendered=True)
    )
    page.set_viewport_size({"width": 900, "height": 600})
    page.set_content("<body></body>")
    page.evaluate(SETUP, [_module(), chart.to_spec()])
    page.locator("#output .gene").first.wait_for()
    yield page
    assert errors == []


def events(page: Page, name: str) -> list:
    return page.evaluate(
        "name => window.sent.filter((m) => m.type === 'event' && m.name === name)"
        ".map((m) => m.detail)",
        name,
    )


def test_widget_draws_at_the_size_of_the_spec(widget_page: Page):
    assert widget_page.locator("#output .gene").count() == 17
    box = widget_page.locator("#output .geneviewer-jupyter").bounding_box()
    assert box["height"] == 300
    assert box["width"] > 800
    assert events(widget_page, "rendered")[0]["height"] == 300


def test_widget_sends_events(widget_page: Page):
    widget_page.locator("#output .gene").nth(2).click()
    widget_page.wait_for_function("window.sent.some((m) => m.name === 'select')")
    assert events(widget_page, "select")[-1]["data"]["name"] == "ophA"
    assert [g["data"]["name"] for g in events(widget_page, "selection")[-1]] == ["ophA"]


def test_zoom_message_moves_the_window(widget_page: Page):
    gene = widget_page.locator("#output .gene").first
    width = gene.bounding_box()["width"]
    widget_page.evaluate(
        "window.model.emit('msg:custom', { type: 'zoom', range: [0, 0.5] })"
    )
    widget_page.wait_for_function(
        "window.sent.some((m) => m.name === 'zoom' && m.detail.action === 'set')"
    )
    assert events(widget_page, "zoom")[-1]["range"] == [0, 0.5]
    assert gene.bounding_box()["width"] == pytest.approx(2 * width, rel=0.1)


def test_a_new_spec_redraws_in_place(widget_page: Page):
    genes = datasets.load("ophA_clusters")
    one = genes[genes["cluster"] == "ophA"]
    spec = GC(one, group="class", height="200px").to_spec()
    widget_page.evaluate(
        "spec => { window.state.spec = spec; window.model.emit('change:spec'); }", spec
    )
    widget_page.wait_for_function(
        "document.querySelectorAll('#output .gene').length === 7"
    )
    assert widget_page.locator("#output .geneviewer-jupyter").count() == 1
    box = widget_page.locator("#output .geneviewer-jupyter").bounding_box()
    assert box["height"] == 200


def test_cleanup_removes_the_chart(widget_page: Page):
    widget_page.evaluate("window.cleanup()")
    assert widget_page.locator("#output .gene").count() == 0
    assert widget_page.evaluate("window.model.count('change:spec')") == 0
    assert widget_page.evaluate("window.model.count('msg:custom')") == 0


def test_two_widgets_each_draw_their_own_chart(widget_page: Page):
    # As in a notebook: each widget loads the module anew
    genes = datasets.load("ophA_clusters")
    one = genes[genes["cluster"] == "ophA"]
    spec = GC(one, group="class", height="200px").to_spec()
    widget_page.evaluate(
        """async ([source, spec]) => {
          const url = URL.createObjectURL(new Blob([source + "\\n// again"], { type: "text/javascript" }));
          const module = await import(url);
          const el = document.createElement("div");
          el.id = "second";
          document.body.appendChild(el);
          const model = { get: () => spec, on() {}, off() {}, send() {} };
          module.default.render({ model, el });
        }""",
        [_module(), spec],
    )
    widget_page.locator("#second .gene").first.wait_for()
    assert widget_page.locator("#second .gene").count() == 7
    assert widget_page.locator("#output .gene").count() == 17
