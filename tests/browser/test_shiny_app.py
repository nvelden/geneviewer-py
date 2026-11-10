"""The example Shiny app in a browser: drawing, events, zoom_to and redraws."""

import json
from pathlib import Path

import pytest

pytest.importorskip("shiny")

from playwright.sync_api import Page, expect  # noqa: E402
from shiny.pytest import create_app_fixture  # noqa: E402

pytestmark = pytest.mark.browser

APP = Path(__file__).resolve().parents[2] / "examples" / "shiny" / "app.py"
app = create_app_fixture(APP)


def events(page: Page) -> dict:
    """Read the events panel: the last value of each input, by event name."""
    text = page.locator("#events").text_content()
    values = {}
    for block in text.split("\n\n"):
        header, _, value = block.partition("\n")
        if header.startswith("input.chart_"):
            values[header[len("input.chart_") : -len("():")]] = json.loads(value)
    return values


def wait_for_event(page: Page, name: str, check=lambda value: True) -> dict:
    page.wait_for_function(
        """([name]) => document.querySelector("#events").textContent
             .includes("input.chart_" + name + "()")""",
        arg=[name],
    )
    for _ in range(50):
        value = events(page).get(name)
        if value is not None and check(value):
            return value
        page.wait_for_timeout(100)
    raise AssertionError(f"input.chart_{name}() is {events(page).get(name)}")


@pytest.fixture
def chart_page(page: Page, app) -> Page:
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(app.url)
    page.locator("#chart svg").first.wait_for()
    yield page
    assert errors == []


def test_app_draws_every_chart(chart_page: Page):
    expect(chart_page.locator(".geneviewer-output")).to_have_count(3)
    expect(chart_page.locator("#chart .gene")).to_have_count(7)
    # Two clusters of ophA_clusters, with a link for each pair of genes of
    # the same class
    expect(chart_page.locator("#links .gene")).to_have_count(17)
    expect(chart_page.locator("#links path.GeneLink")).to_have_count(14)
    # The five transcripts of BRCA1_splice_variants, with their introns
    expect(chart_page.locator("#transcript .exon").first).to_be_visible()
    expect(chart_page.locator("#transcript .intron").first).to_be_visible()
    expect(chart_page.locator("#transcript .five-prime-label")).to_have_count(5)
    # The first draw reports rendered and the window of the zoom slider
    assert wait_for_event(chart_page, "rendered")["chartId"] == "chart"
    zoom = wait_for_event(chart_page, "zoom")
    assert zoom["action"] == "render"
    assert zoom["range"] is None


def test_clicking_a_gene_sets_the_select_input(chart_page: Page):
    chart_page.locator("#chart .gene").nth(2).click()
    select = wait_for_event(chart_page, "select")
    assert select["component"] == "genes"
    assert select["data"]["name"] == "ophA"
    selection = wait_for_event(chart_page, "selection")
    assert [item["data"]["name"] for item in selection] == ["ophA"]
    # The same gene again still reaches the server: priority "event"
    chart_page.locator("#chart .gene").nth(2).click()
    wait_for_event(chart_page, "selection", lambda value: value == [])


def test_zoom_to_moves_the_window(chart_page: Page):
    wait_for_event(chart_page, "zoom")
    gene = chart_page.locator("#chart .gene").first
    width = gene.bounding_box()["width"]
    chart_page.click("#first_half")
    zoom = wait_for_event(chart_page, "zoom", lambda value: value["action"] == "set")
    assert zoom["range"] == [0, 0.5]
    assert zoom["zoomed"] is True
    # Half the range in the same space: each gene is drawn about twice as wide
    assert gene.bounding_box()["width"] == pytest.approx(2 * width, rel=0.1)
    chart_page.click("#everything")
    zoom = wait_for_event(
        chart_page,
        "zoom",
        lambda value: value["action"] == "set" and not value["zoomed"],
    )
    assert zoom["range"] is None
    assert gene.bounding_box()["width"] == pytest.approx(width, rel=0.01)


def test_a_new_chart_redraws_in_place(chart_page: Page):
    first = wait_for_event(chart_page, "rendered")
    chart = chart_page.locator("#chart")
    gene = chart.locator(".gene").first
    parts = "el => el.querySelectorAll('svg, div').length"
    style = "el => getComputedStyle(el).backgroundColor"
    fill = "el => getComputedStyle(el).fill"
    background = chart.evaluate(style)
    color = gene.evaluate(fill)
    width = gene.bounding_box()["width"]
    chart_page.click("#first_half")
    wait_for_event(chart_page, "zoom", lambda value: value["action"] == "set")
    count = chart.evaluate(parts)
    chart_page.select_option("#theme", "dark")
    wait_for_event(chart_page, "rendered", lambda value: value["id"] != first["id"])
    # One chart in the output, in the new theme
    assert chart.evaluate(parts) == count
    expect(chart.locator(".gene")).to_have_count(7)
    assert chart.evaluate(style) != background
    assert gene.evaluate(fill) != color
    # setSpec() on the same chart, not a new chart: with the same data, the
    # window picked before is kept (the theme's gene shapes differ a little)
    assert gene.bounding_box()["width"] > 1.5 * width
    assert events(chart_page)["zoom"]["range"] == [0, 0.5]
