"""The readers give R's results on the fixtures of dev/export_reader_fixtures.R.

Each fixture in expected/ names the calls, their input (a file in inputs/,
or a table) and R's result, warnings and error. The test runs the same
calls from this folder, so paths in messages are the same, and compares:
the same columns in the same order, the same types (character as
``string``, integer as ``Int64``, numeric as ``float64``, logical as
``boolean``, list as lists), the same values (NA as ``<NA>`` or ``NaN``),
the same warnings and the same error message.
"""

import json
import math
import warnings
from pathlib import Path

import pandas as pd
import pytest

import geneviewer
from geneviewer.readers._common import make_unique

HERE = Path(__file__).resolve().parent
FIXTURES = sorted((HERE / "expected").glob("*.json"))

# Fixtures where R is wrong and Python does not copy it (none: B65 and B66
# were fixed in R)
KNOWN: dict[str, str] = {}

DTYPES = {
    "character": "string",
    "integer": "Int64",
    "numeric": "float64",
    "logical": "boolean",
    "list": "object",
}


def to_frame(table: dict) -> pd.DataFrame:
    """Build a DataFrame with R's types from a fixture table."""
    columns = {}
    for column in table["columns"]:
        name, kind = column["name"], column["type"]
        values = table["data"][name]
        if kind == "numeric":
            values = [math.nan if v is None else v for v in values]
        if kind == "list":
            columns[name] = pd.Series(values, dtype=object)
        else:
            columns[name] = pd.array(values, dtype=DTYPES[kind])
    return pd.DataFrame(columns)


def check_table(actual: pd.DataFrame, expected: dict) -> None:
    names = [c["name"] for c in expected["columns"]]
    assert list(actual.columns) == names
    nrows = len(expected["data"][names[0]]) if names else 0
    assert len(actual) == nrows
    for column in expected["columns"]:
        name, kind = column["name"], column["type"]
        series = actual[name]
        assert str(series.dtype) == DTYPES[kind], f"{name}: {series.dtype}, R {kind}"
        for i, (a, e) in enumerate(
            zip(series.tolist(), expected["data"][name], strict=True)
        ):
            assert same(a, e), f"{name}[{i}]: {a!r}, R {e!r}"


def same(actual, expected) -> bool:
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            return False
        return all(same(a, e) for a, e in zip(actual, expected, strict=True))
    if expected is None:
        return (
            actual is None
            or actual is pd.NA
            or (isinstance(actual, float) and math.isnan(actual))
        )
    if isinstance(expected, bool):
        return actual is not pd.NA and bool(actual) is expected
    if isinstance(expected, (int, float)):
        if actual is None or actual is pd.NA:
            return False
        return math.isclose(float(actual), float(expected), rel_tol=1e-12)
    return actual == expected


def gbk_from_r(value: dict) -> dict:
    """R's read_gbk() with features as dicts, repeated names as in Python."""
    result = {}
    for name, record in value.items():
        sections = {}
        for section, text in record.items():
            if section != "FEATURES":
                sections[section] = text
                continue
            features = {}
            for key, items in (text or {}).items():
                if items is None:
                    features[key] = None
                    continue
                features[key] = [
                    dict(
                        zip(
                            make_unique(p[0] for p in pairs),
                            [p[1] for p in pairs],
                            strict=True,
                        )
                    )
                    for pairs in items
                ]
            sections[section] = features
        result[name] = sections
    return result


def run(fixture: dict):
    value = fixture["input"]
    if fixture["data"] is not None:
        value = to_frame(fixture["data"])
    for step in fixture["calls"]:
        function = getattr(geneviewer, step["fun"])
        args = dict(step["args"])
        if step["fun"] == "genbank_to_fasta":
            args["output_dir"] = None
        value = function(value, **args)
    return value


def message(text: str) -> str:
    return text.replace("\\", "/")


@pytest.mark.parametrize("path", FIXTURES, ids=[p.stem for p in FIXTURES])
def test_reader_matches_r(path, monkeypatch, tmp_path):
    if path.stem in KNOWN:
        pytest.xfail(KNOWN[path.stem])
    fixture = json.loads(path.read_text(encoding="utf-8"))
    monkeypatch.chdir(HERE)
    if fixture["calls"][-1]["fun"] == "genbank_to_fasta":
        # Write next to a copy of the input, not into inputs/
        copy = tmp_path / Path(fixture["input"]).name
        copy.write_bytes((HERE / fixture["input"]).read_bytes())
        fixture["input"] = str(copy)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if fixture["error"] is not None:
            with pytest.raises(Exception) as error:
                run(fixture)
            assert message(str(error.value)) == fixture["error"]
            return
        result = run(fixture)
    found = [message(str(w.message)) for w in caught if w.category is UserWarning]
    assert found == fixture["warnings"]

    last = fixture["calls"][-1]["fun"]
    if last == "read_gbk":
        assert result == gbk_from_r(fixture["result"])
    elif last == "genbank_to_fasta":
        lines = Path(result).read_text(encoding="utf-8").split("\n")[:-1]
        assert lines == fixture["result"]
    else:
        check_table(result, fixture["result"])


def test_every_reader_has_fixtures():
    funs = {
        json.loads(p.read_text(encoding="utf-8"))["calls"][-1]["fun"] for p in FIXTURES
    }
    assert funs == {
        "read_gbk",
        "gbk_features_to_df",
        "genbank_to_fasta",
        "read_gff",
        "read_fasta",
        "read_bed",
        "get_introns",
    }
