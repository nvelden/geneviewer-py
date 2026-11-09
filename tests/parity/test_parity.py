"""Specs built in Python against R's fixtures in geneviewer-js.

Each file ``cases/<fixture>.py`` defines ``build(data)`` (or ``build(data,
links)``), which makes the fixture's chart with the Python API. ``data`` is
the fixture's ``spec.data`` as a data frame, ``links`` the data of its
first links entry. The spec the chart makes must equal ``spec.json`` under
the rules of ``compare.py``. Every fixture spec must also be valid against
the vendored schema.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from pathlib import Path

import pandas as pd
import pytest

from geneviewer.chart import validate_spec

from .compare import differences

pytestmark = pytest.mark.parity

CASES = Path(__file__).parent / "cases"

# Fixtures that are not cases yet, by the step that adds them
LATER = {}

# Known differences, waiting for a decision: fixture, reason
DIFFERENCES = {
    "rd-GC-genes-1": (
        "R keeps the names of options passed through ..., so prevent_overlap "
        "and gene_overlap_spacing stay snake_case; Python writes options in "
        "camelCase. Neither is an option of the library."
    ),
}


def _frame(columns: dict) -> pd.DataFrame:
    # R writes a column of one row as a scalar
    return pd.DataFrame(
        {name: v if isinstance(v, list) else [v] for name, v in columns.items()}
    )


def _load(path: Path):
    if str(CASES) not in sys.path:
        sys.path.insert(0, str(CASES))
    spec = importlib.util.spec_from_file_location(f"case_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fixtures(js_checkout: Path) -> Path:
    return js_checkout / "test" / "fixtures"


def _links_data(spec: dict) -> pd.DataFrame | None:
    links = spec.get("links")
    if links is None:
        return None
    first = links if isinstance(links, dict) else links[0]
    return _frame(first["data"]) if "data" in first else None


@pytest.mark.parametrize(
    "case", sorted(p.stem for p in CASES.glob("*.py") if not p.stem.startswith("_"))
)
def test_case(case, js_checkout):
    if case in DIFFERENCES:
        pytest.xfail(DIFFERENCES[case])
    expected = json.loads(
        (_fixtures(js_checkout) / case / "spec.json").read_text(encoding="utf-8")
    )
    build = _load(CASES / f"{case}.py").build
    data = _frame(expected["data"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if build.__code__.co_argcount > 1:
            chart = build(data, _links_data(expected))
        else:
            chart = build(data)
    actual = chart.to_spec()
    diffs = differences(expected, actual)
    assert not diffs, "\n".join(
        f"{where}: R {json.dumps(e)[:120]}, Python {json.dumps(a)[:120]}"
        for where, e, a in diffs[:20]
    )


def test_every_fixture_spec_is_valid(js_checkout):
    invalid = []
    for path in sorted(_fixtures(js_checkout).glob("*/spec.json")):
        try:
            validate_spec(json.loads(path.read_text(encoding="utf-8")))
        except ValueError as error:
            invalid.append(f"{path.parent.name}: {error}")
    assert not invalid, "\n".join(invalid)


def test_every_rd_and_cov_fixture_has_a_case(js_checkout):
    cases = {p.stem for p in CASES.glob("*.py")}
    missing = [
        path.parent.name
        for path in sorted(_fixtures(js_checkout).glob("*/spec.json"))
        if path.parent.name.startswith(("rd-", "cov-"))
        and path.parent.name not in cases
        and path.parent.name not in LATER
    ]
    assert not missing, f"no case for: {', '.join(missing)}"
