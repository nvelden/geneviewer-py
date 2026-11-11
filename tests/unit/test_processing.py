"""Links and normalize, with the values R gives (geneviewer m6-js-library)."""

import warnings

import pandas as pd
import pytest

from geneviewer import GC
from geneviewer.data import chart_state
from geneviewer.processing import get_links


def links(data, group="grp", **kwargs):
    return get_links(chart_state(GC(data, cluster="cl")._spec), group, **kwargs)


def genes(groups, clusters):
    n = len(groups)
    return pd.DataFrame(
        {
            "start": [10 * (i + 1) for i in range(n)],
            "end": [10 * (i + 1) + 5 for i in range(n)],
            "name": [f"g{i}" for i in range(n)],
            "grp": groups,
            "cl": clusters,
        }
    )


# get_links ------------------------------------------------------------------


def test_links_by_group_columns():
    data = genes(["A", "B", "A", "C"], [1, 1, 2, 2])
    result = links(data)
    assert list(result["name1"]) == ["g0"]
    assert list(result["name2"]) == ["g2"]
    # As R's merge(): the group column once, the cluster column renamed
    assert "grp1" in result.columns and "grp2" not in result.columns
    assert list(result["cluster1"]) == ["1"] and list(result["cluster2"]) == ["2"]
    assert {"rowID1", "rowID2", "strand1", "strand2"} <= set(result.columns)


@pytest.mark.parametrize("position", [0, 5])
def test_links_with_another_column_named_cluster(position):
    # The clusters are "cl"; a column "cluster" before or after it is not
    # read in their place, and is left out of the links
    data = genes(["A", "B", "A", "B"], [1, 1, 2, 2])
    data.insert(position, "cluster", ["x", "x", "x", "x"])
    result = links(data)
    assert list(result["name1"]) == ["g0", "g1"]
    assert list(result["name2"]) == ["g2", "g3"]
    assert list(result["cluster1"]) == ["1", "1"]
    assert list(result["cluster2"]) == ["2", "2"]
    assert list(result.columns).count("cluster1") == 1
    with pytest.raises(ValueError, match="cannot be the group of the links"):
        links(data, group="cluster")


def test_links_rows_in_rs_merge_order():
    # R's merge() shell-sorts the rows of a group: with 9 rows they do not
    # keep their order (R: g0 g5 g2 g8 g3 g1 g6 g7 g4)
    data = genes(["M", "N", *"MNMNNMNNM"], [1, 1] + [2] * 9)
    data["name"] = ["a", "b"] + [f"g{i}" for i in range(9)]
    result = links(data)
    assert list(result["name2"]) == [
        "g0",
        "g5",
        "g2",
        "g8",
        "g3",
        "g1",
        "g6",
        "g7",
        "g4",
    ]


def test_links_leave_out_missing_and_no_hit():
    data = genes(["A", None, "No Hit", "A", None, "No Hit"], [1, 1, 1, 2, 2, 2])
    assert list(links(data)["name2"]) == ["g3"]


def test_links_value_pairs_swap_for_the_next_clusters():
    # B is not in cluster 1, so value1 and value2 swap, and stay swapped for
    # clusters 2 and 3 (R: g0-g1, g2-g3)
    data = genes(["A", "B", "A", "B"], [1, 2, 2, 3])
    result = links(data, value1="B", value2="A")
    assert list(result["name1"]) == ["g0", "g2"]
    assert list(result["name2"]) == ["g1", "g3"]
    assert {"grp1", "grp2"} <= set(result.columns)


def test_links_numeric_groups():
    data = genes([1.0, 2.0, 2.0, 1.0], [1, 1, 2, 2])
    assert list(links(data)["name2"]) == ["g3", "g2"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"group": "x"}, "Column x not found in chart data."),
        ({"value1": "Z"}, "Some values in 'value1' are not present"),
        ({"value2": "Z"}, "Some values in 'value2' are not present"),
        ({"cluster": 1}, "At least two clusters need to be provided."),
        # value1 and value2 of different lengths link nothing, as in R
        ({"value1": ["A", "B"], "value2": "A"}, "No links found"),
    ],
)
def test_links_errors(kwargs, message):
    data = genes(["A", "B", "A", "C"], [1, 1, 2, 2])
    with pytest.raises(ValueError, match=message):
        links(data, **kwargs)


def test_links_need_clusters():
    state = chart_state(GC(genes(["A", "A"], [1, 1]))._spec)
    with pytest.raises(ValueError, match="Please define cluster"):
        get_links(state, "grp")


# normalize ------------------------------------------------------------------

POSITIONS = pd.DataFrame(
    {
        "start": [100, 300, 700, 1000],
        "end": [200, 250, 900, 1100],
        "grp": ["A", "B", "A", "C"],
    }
)


def normalized(**kwargs):
    return GC(POSITIONS, group="grp").normalize(**kwargs).to_spec()["data"]


def test_normalize_keeps_gene_lengths_and_strand():
    data = normalized()
    assert data["start"] == [100, 433, 617, 1000]
    assert data["end"] == [200, 383, 817, 1100]
    assert data["original_start"] == [100, 300, 700, 1000]
    assert data["original_end"] == [200, 250, 900, 1100]


def test_normalize_mean_length_and_gap():
    data = normalized(preserve_gene_length=False, gap=0.1)
    assert data["start"] == [100, 550, 650, 925]
    assert data["end"] == [275, 375, 825, 1100]


def test_normalize_keeps_the_range_of_the_cluster():
    # Evenly spaced genes stay where they are; a given gap is the gap drawn
    even = pd.DataFrame({"start": [10, 50, 90], "end": [40, 80, 120]})
    data = GC(even).normalize().to_spec()["data"]
    assert data["start"] == [10, 50, 90]
    assert data["end"] == [40, 80, 120]
    data = GC(even).normalize(gap=0.1).to_spec()["data"]
    assert data["start"] == [10, 51, 92]
    assert data["end"] == [40, 81, 122]
    single = pd.DataFrame({"start": [300], "end": [400]})
    data = GC(single).normalize(preserve_gene_length=False).to_spec()["data"]
    assert (data["start"], data["end"]) == ([300], [400])


def test_normalize_custom_gaps_only():
    data = normalized(normalize=False, custom_gaps={"A": 2})
    assert data["start"] == [100, 350, 750, 1150]
    assert data["end"] == [200, 300, 950, 1250]
    assert data["original_start"] == [100, 300, 700, 1000]


def test_normalize_one_cluster_and_other_start_column():
    data = pd.DataFrame(
        {
            "from": [10, 90, 130, 10, 90, 130],
            "to": [40, 120, 160, 40, 120, 160],
            "cl": [1, 1, 1, 2, 2, 2],
        }
    )
    spec = GC(data, start="from", end="to", cluster="cl").normalize(cluster=2)
    result = spec.to_spec()["data"]
    assert result["start"] == [10, 90, 130, 10, 70, 130]
    assert result["from"] == result["start"]
    assert result["original_start"] == [None, None, None, 10, 90, 130]


def test_links_made_before_normalize_move_with_their_genes():
    data = pd.DataFrame(
        {
            "s": [10, 100, 150, 10, 60, 110],
            "e": [40, 130, 180, 40, 90, 140],
            "name": ["A", "B", "C"] * 2,
            "cl": [1, 1, 1, 2, 2, 2],
        }
    )
    chart = GC(data, start="s", end="e", cluster="cl", group="name")
    before = (
        chart.links("name")
        .links("name", value1="A", value2="A")
        .normalize(cluster=1)
        .to_spec()["links"]
    )
    after = (
        chart.normalize(cluster=1)
        .links("name")
        .links("name", value1="A", value2="A")
        .to_spec()["links"]
    )
    columns = ["rowID1", "rowID2", "start1", "end1", "start2", "end2", "s1", "e1"]
    for old, new in zip(before, after, strict=True):
        assert {c: old["data"][c] for c in columns} == {
            c: new["data"][c] for c in columns
        }
    assert before[0]["data"]["start1"] == [10, 80, 150]

    one = chart.links("name").normalize(cluster=1).to_spec()["links"]
    assert one["data"]["end1"] == [40, 110, 180]


def test_links_given_as_data_keep_their_positions_after_normalize():
    data = pd.DataFrame(
        {"start": [10, 100, 150, 10], "end": [40, 130, 180, 40], "cl": [1, 1, 1, 2]}
    )
    alignment = pd.DataFrame(
        {
            "cluster1": ["1"],
            "cluster2": ["2"],
            "start1": [100],
            "end1": [120],
            "start2": [10],
            "end2": [30],
            "identity": [90],
        }
    )
    with pytest.warns(UserWarning, match="keep their positions"):
        spec = GC(data, cluster="cl").links(data=alignment).normalize().to_spec()
    assert spec["links"]["data"]["start1"] == [100]
    assert "links" not in GC(data, cluster="cl").normalize().to_spec()


def test_normalize_warnings():
    with pytest.warns(UserWarning, match="Please define cluster"):
        GC(POSITIONS).normalize(cluster=1)
    with pytest.warns(UserWarning, match="Custom gaps will not be applied"):
        GC(POSITIONS).normalize(custom_gaps={"A": 2})
    with pytest.warns(UserWarning, match="Z not found in the specified group"):
        GC(POSITIONS, group="grp").normalize(custom_gaps={"Z": 2})
    with pytest.raises(ValueError, match="group column not found"):
        GC(POSITIONS).normalize(group="x")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        GC(POSITIONS, group="grp").normalize()


def test_a_chart_without_rows_takes_transcript_and_normalize():
    data = pd.DataFrame(
        {
            "start": pd.Series([], dtype=float),
            "end": pd.Series([], dtype=float),
            "name": pd.Series([], dtype=str),
            "cl": pd.Series([], dtype=str),
        }
    )
    spec = GC(data, cluster="cl").normalize().labels("name").to_spec()
    assert spec["labels"] == {"label": "name"}
    spec = GC(data).transcript(transcript="cl").to_spec()
    assert spec["data"]["start"] == []
