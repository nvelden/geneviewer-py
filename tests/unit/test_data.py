import datetime as dt
import math

import numpy as np
import pandas as pd
import pytest

from geneviewer import data as d

GENES = pd.DataFrame(
    {
        "start": [10, 90, 130, 170, 210],
        "end": [40, 120, 160, 200, 240],
        "name": ["Gene 1", "Gene 3", "Gene 4", "Gene 5", "Gene 6"],
        "group": ["A", "B", "B", "A", "C"],
        "cluster": [1, 1, 2, 2, 2],
    }
)


# prepare_data -------------------------------------------------------------


def test_prepare_data_adds_columns_and_keeps_the_input():
    before = GENES.copy()
    out = d.prepare_data(GENES, cluster="cluster", group="group")
    pd.testing.assert_frame_equal(GENES, before)
    assert list(out.columns) == [*GENES.columns, "rowID", "strand"]
    assert out["rowID"].tolist() == [1, 2, 3, 4, 5]
    assert out["cluster"].tolist() == ["1", "1", "2", "2", "2"]
    assert out["strand"].tolist() == ["forward"] * 5


def test_prepare_data_other_start_and_end_columns():
    data = pd.DataFrame({"from": [5, 1], "to": [1, 5]})
    out = d.prepare_data(data, start="from", end="to")
    assert list(out.columns) == ["from", "to", "rowID", "start", "end", "strand"]
    assert out["strand"].tolist() == ["reverse", "forward"]


def test_prepare_data_leaves_out_rows_without_start_or_end():
    data = pd.DataFrame({"start": [1, None, 3], "end": [2, 3, None]})
    with pytest.warns(UserWarning, match="filtered out"):
        out = d.prepare_data(data)
    assert out["rowID"].tolist() == [1]
    # No rows left: the columns stay
    with pytest.warns(UserWarning):
        empty = d.prepare_data(pd.DataFrame({"start": [None], "end": [1]}))
    assert len(empty) == 0
    assert list(empty.columns) == ["start", "end", "rowID", "strand"]


@pytest.mark.parametrize(
    ("args", "message"),
    [
        ({"start": "x"}, "start column not found in data"),
        ({"end": "x"}, "end column not found in data"),
        ({"cluster": "x"}, "cluster column not found in data"),
        ({"group": "x"}, "group column not found in data"),
        ({"strand": "x"}, "strand column not found in data"),
    ],
)
def test_prepare_data_checks_columns(args, message):
    with pytest.raises(ValueError, match=message):
        d.prepare_data(GENES, **args)


def test_prepare_data_needs_a_data_frame():
    with pytest.raises(TypeError):
        d.prepare_data({"start": [1], "end": [2]})


# add_strand ---------------------------------------------------------------


def test_add_strand_from_a_column_swaps_start_and_end():
    data = pd.DataFrame(
        {"start": [3, 2, 5, 1], "end": [2, 4, 6, 9], "strand": [1, -1, "+", 0]}
    )
    out = d.add_strand(data)
    assert out["strand"].tolist() == ["forward", "reverse", "forward", "reverse"]
    assert out["start"].tolist() == [2, 4, 5, 9]
    assert out["end"].tolist() == [3, 2, 6, 1]


def test_add_strand_values_as_r_reads_them():
    data = pd.DataFrame({"start": [1, 2], "end": [2, 3], "dir": [1.0, -1.0]})
    assert d.add_strand(data, "dir")["strand"].tolist() == ["forward", "reverse"]
    data["dir"] = ["complement", "x"]
    with pytest.raises(ValueError, match="Invalid strand values"):
        d.add_strand(data, "dir")


def test_add_strand_no_rows():
    assert d.add_strand(None) is None
    empty = d.add_strand(GENES.iloc[0:0])
    assert len(empty) == 0
    assert list(empty.columns) == [*GENES.columns, "strand"]


# chart_fields, clusters ---------------------------------------------------


def test_chart_fields():
    spec = {"data": GENES, "fields": {"cluster": "cluster", "group": "group"}}
    state = d.chart_fields(spec)
    assert state["start"] == "start" and state["group"] == "group"
    assert state["clusters"] == ["1", "2"]
    assert state["columns"] == [*GENES.columns, "rowID", "strand"]
    assert d.chart_fields({"data": GENES})["clusters"] == ["cluster"]


def test_chart_fields_after_transcript():
    spec = {
        "data": GENES.assign(transcript=["t1", "t1", "t2", "t2", "t3"]),
        "transcript": [{"transcript": "transcript", "strand": "s"}, {}],
    }
    state = d.chart_fields(spec)
    assert state["cluster"] == "transcript" and state["strand"] == "s"
    assert state["clusters"] == ["t1", "t2", "t3"]


def test_chart_fields_no_rows():
    data = pd.DataFrame({"start": [None], "end": [1.0], "cluster": ["a"]})
    state = d.chart_fields({"data": data, "fields": {"cluster": "cluster"}})
    assert state["columns"] == ["start", "end", "cluster", "rowID", "strand"]
    assert state["clusters"] == []
    assert d.chart_fields({"data": data})["clusters"] == []


def test_get_updated_clusters():
    clusters = ["a", "b", "c"]
    assert d.get_updated_clusters(clusters, None) == clusters
    assert d.get_updated_clusters(clusters, 2) == ["b"]
    assert d.get_updated_clusters(clusters, [3, 1]) == ["c", "a"]
    assert d.get_updated_clusters(clusters, "b") == ["b"]
    with pytest.warns(UserWarning, match="out of range"):
        assert d.get_updated_clusters(clusters, [1, 4]) is None
    with pytest.warns(UserWarning, match="not valid"):
        assert d.get_updated_clusters(clusters, ["a", "x"]) is None


def test_cluster_rows():
    spec = {"data": GENES, "fields": {"cluster": "cluster"}}
    rows = d.cluster_rows(d.chart_state(spec), "2")
    assert rows["rowID"].tolist() == [3, 4, 5]


# jsonable -----------------------------------------------------------------


def test_jsonable_data_frame_as_columns():
    data = pd.DataFrame(
        {
            "i": pd.array([1, None], dtype="Int64"),
            "f": [1.5, math.nan],
            "inf": [math.inf, 1.0],
            "s": ["a", None],
            "b": pd.array([True, None], dtype="boolean"),
            "c": pd.Categorical(["y", None], categories=["y", "x"]),
            "day": pd.to_datetime(["2020-01-02", None]),
            "time": pd.to_datetime(["2020-01-02 03:04:05", "2020-01-02 00:00:00"]),
            "date": [dt.date(2020, 1, 2), None],
        }
    )
    assert d.jsonable(data) == {
        "i": [1, None],
        "f": [1.5, None],
        "inf": [None, 1.0],
        "s": ["a", None],
        "b": [True, None],
        "c": ["y", None],
        "day": ["2020-01-02", None],
        "time": ["2020-01-02T03:04:05Z", "2020-01-02T00:00:00Z"],
        "date": ["2020-01-02", None],
    }


def test_jsonable_values():
    value = {
        "a": (1, np.int64(2)),
        "b": np.array([0.5, np.nan]),
        "c": pd.Series(["x"]),
        "d": np.bool_(True),
        "e": [{"x": np.float64(1.0)}],
    }
    assert d.jsonable(value) == {
        "a": [1, 2],
        "b": [0.5, None],
        "c": ["x"],
        "d": True,
        "e": [{"x": 1.0}],
    }
    out = d.jsonable(value)
    assert type(out["a"][1]) is int and type(out["d"]) is bool
