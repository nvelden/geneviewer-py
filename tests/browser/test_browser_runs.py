import pytest

pytestmark = pytest.mark.browser


def test_browser_runs(page):
    page.set_content("<p id='ok'>ok</p>")
    assert page.text_content("#ok") == "ok"
