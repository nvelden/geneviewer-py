"""The example datasets, compared with the files R exported."""

import csv
import json
import math
from pathlib import Path

import pandas as pd
import pytest

import geneviewer
from geneviewer import datasets

FOLDER = Path(datasets.__file__).parent
INDEX = json.loads((FOLDER / "datasets.json").read_text(encoding="utf-8"))
DTYPES = {"character": "string", "integer": "Int64", "numeric": "float64"}


def test_names():
    assert datasets.names() == sorted(INDEX)
    assert len(datasets.names()) == 11
    assert sorted(p.stem for p in FOLDER.glob("*.csv")) == datasets.names()
    assert geneviewer.datasets is datasets


def test_titles():
    titles = datasets.titles()
    assert list(titles) == datasets.names()
    assert all(titles.values())


@pytest.mark.parametrize("name", sorted(INDEX))
def test_dataset_matches_export(name):
    data = datasets.load(name)
    entry = INDEX[name]
    columns = [c["name"] for c in entry["columns"]]
    assert list(data.columns) == columns
    assert len(data) == entry["rows"]
    with open(FOLDER / f"{name}.csv", encoding="utf-8", newline="") as file:
        header, *rows = csv.reader(file)
    assert header == columns
    for j, column in enumerate(entry["columns"]):
        series = data[column["name"]]
        assert str(series.dtype) == DTYPES[column["type"]], column["name"]
        for value, row in zip(series.tolist(), rows, strict=True):
            text = row[j]
            if text == "NA":
                assert value is pd.NA or math.isnan(value)
            elif column["type"] == "character":
                assert value == text
            elif column["type"] == "integer":
                assert value == int(text)
            else:
                assert value == float(text)


def test_empty_text_is_not_missing():
    # R has both "" and NA in this column
    gene = datasets.load("erythromycin_BlastP")["Gene"]
    assert int(gene.isna().sum()) == 48
    assert int((gene == "").sum()) == 7


def test_load_gives_a_copy():
    first = datasets.load("ophA_clusters")
    first.loc[0, "name"] = "changed"
    assert datasets.load("ophA_clusters").loc[0, "name"] == "ophB1"


def test_unknown_dataset():
    with pytest.raises(ValueError, match="There is no dataset 'nope'"):
        datasets.load("nope")


def test_chart_from_dataset():
    data = datasets.load("ophA_clusters")
    chart = geneviewer.GC(data, cluster="cluster", group="class")
    assert chart is not None
