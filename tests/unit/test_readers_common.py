"""R's rules the readers copy, checked against R's results (R 4.5)."""

import math

import pandas as pd
import pytest

from geneviewer.readers._common import (
    bind_rows,
    location_span,
    make_names,
    read_tsv,
    type_convert,
    widen_attributes,
)


@pytest.mark.parametrize(
    ("values", "dtype"),
    [
        (["T", "F"], "boolean"),
        (["TRUE", ""], "boolean"),
        (["true"], "string"),
        (["T", "false"], "string"),
        ([" TRUE"], "string"),
        ([" 1", "2"], "Int64"),
        (["1 "], "float64"),
        (["1", "", "NA"], "Int64"),
        (["+5", "-0"], "Int64"),
        (["2147483648"], "float64"),
        (["0x10", "3"], "float64"),
        (["1e", ".5", "5.", "Infinity", "nan", "1e400"], "float64"),
        (["1L"], "string"),
        (["1,5"], "string"),
        (["1d5"], "string"),
        (["-"], "string"),
        (["."], "string"),
        (["", "NA", None], "boolean"),
        (["  "], "boolean"),
        (["", "a"], "string"),
    ],
)
def test_type_convert_as_r(values, dtype):
    assert str(type_convert(values).dtype) == dtype


def test_type_convert_values():
    assert type_convert(["0x10", "3"]).tolist() == [16.0, 3.0]
    assert type_convert(["0x1p3"]).tolist() == [8.0]
    assert type_convert(["1e400", "-inf"]).tolist() == [math.inf, -math.inf]
    assert type_convert(["", "a", "NA"]).tolist() == ["", "a", pd.NA]


def test_read_tsv_fills_short_rows_and_keeps_rows_of_spaces():
    data = read_tsv(["a\t1\tx", "b\t2", "", "  ", "c\t\tz"], ["V1", "V2", "V3"])
    assert data["V1"].tolist() == ["a", "b", "  ", "c"]
    assert data["V2"].tolist() == [1, 2, pd.NA, pd.NA]
    assert data["V3"].tolist() == ["x", "", "", "z"]
    with pytest.raises(ValueError, match="more columns than column names"):
        read_tsv(["a\t1\tx\textra"], ["V1", "V2", "V3"])


def test_make_names_as_r():
    names = ["a", "b", "b", "function", "db-xref", "1x", "NA", "TRUE", "_y", ".z"]
    names += [".2a", "a.1", None, None]
    assert make_names(names) == [
        "a",
        "b",
        "b.1",
        "function.",
        "db.xref",
        "X1x",
        "NA.",
        "TRUE.",
        "X_y",
        ".z",
        "X.2a",
        "a.1",
        "NA..1",
        "NA..2",
    ]


def test_widen_attributes_joins_repeated_keys():
    ids, data = widen_attributes(
        ["2", "1", "2", "2", "3", "3", "3"],
        ["Note", "gene", "Note", "Note", "gene", "Note", "Note"],
        ["first", "abc", "second", "first", None, None, None],
    )
    assert ids == ["2", "1", "3"]
    assert data["Note"].tolist() == ["first,second,first", pd.NA, pd.NA]
    assert data["gene"].tolist() == [pd.NA, "abc", pd.NA]
    with pytest.raises(ValueError, match="empty"):
        widen_attributes(["1"], [""], ["x"])
    with pytest.raises(ValueError, match="unique"):
        widen_attributes(["1", "2"], [None, "NA"], ["x", "y"])


@pytest.mark.parametrize(
    ("location", "span"),
    [
        (
            "join(complement(5441..5600),complement(5095..5373),complement(5024..5040))",
            (5024, 5600),
        ),
        ("complement(join(6915..7069,7157..7232))", (6915, 7232)),
        ("join(<10..120,200..>300)", (10, 300)),
        ("join(30883..31317,1..108)", (30883, 31317)),
        ("complement(123)", (123, 123)),
    ],
)
def test_location_span(location, span):
    assert location_span(location) == span
    assert all(math.isnan(v) for v in location_span(None))


def test_bind_rows_combines_types():
    a = pd.DataFrame(
        {"x": pd.array([1], dtype="Int64"), "y": pd.array([pd.NA], dtype="boolean")}
    )
    b = pd.DataFrame({"x": [2.5], "y": pd.array(["s"], dtype="string"), "z": [1.0]})
    result = bind_rows([a, b])
    assert list(result.columns) == ["x", "y", "z"]
    assert str(result["x"].dtype) == "float64"
    assert str(result["y"].dtype) == "string"
    assert result["z"].isna().tolist() == [True, False]
    with pytest.raises(TypeError, match="Can't combine"):
        bind_rows(
            [
                pd.DataFrame({"x": [1.0]}),
                pd.DataFrame({"x": pd.array(["a"], dtype="string")}),
            ]
        )
