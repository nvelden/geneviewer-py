import pytest

pytestmark = pytest.mark.parity


def test_fixtures_found(js_checkout):
    specs = list((js_checkout / "test" / "fixtures").glob("*/spec.json"))
    assert specs, f"no fixture specs under {js_checkout}"
