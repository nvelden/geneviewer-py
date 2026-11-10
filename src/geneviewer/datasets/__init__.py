"""Example datasets, exported from the R package.

`names` lists them and `load` returns one as a DataFrame with R's column
types. The files are written by ``dev/export_datasets.R`` of the R package:
a CSV file per dataset, and ``datasets.json`` with the titles and the column
types.
"""

from __future__ import annotations

import json
import math
from functools import cache
from importlib import resources

import pandas as pd

__all__ = ["load", "names", "titles"]

# R's column types as pandas types
_DTYPES = {"character": "string", "integer": "Int64", "numeric": "float64"}


@cache
def _index() -> dict[str, dict]:
    file = resources.files(__name__) / "datasets.json"
    return json.loads(file.read_text(encoding="utf-8"))


def names() -> list[str]:
    """Return the names of the example datasets.

    Returns
    -------
    list of str
        The names, in alphabetical order.

    Examples
    --------
    >>> from geneviewer import datasets
    >>> datasets.names()[:2]
    ['BRCA1_splice_variants', 'coronavirus_genomes']
    """
    return list(_index())


def titles() -> dict[str, str]:
    """Return what each example dataset holds.

    Returns
    -------
    dict of str to str
        The title of each dataset, by name.

    Examples
    --------
    >>> from geneviewer import datasets
    >>> datasets.titles()["ophA_clusters"]
    'ophA Gene Cluster from Omphalotus olearius'
    """
    return {name: entry["title"] for name, entry in _index().items()}


def load(name: str) -> pd.DataFrame:
    """Return an example dataset.

    The columns have the types of the R package's dataset: text is pandas
    ``string``, whole numbers ``Int64`` and other numbers ``float64``. A
    missing value is ``<NA>``, or ``NaN`` in a ``float64`` column.

    Parameters
    ----------
    name : str
        The dataset, one of `names`.

    Returns
    -------
    pandas.DataFrame
        The dataset; a new copy on each call.

    Raises
    ------
    ValueError
        When there is no dataset of that name.

    Examples
    --------
    >>> from geneviewer import datasets
    >>> datasets.load("ophA_clusters").head(3)
        name   start      end              class cluster
    0  ophB1  2522.0   4276.0      Monooxygenase    ophA
    1   ophC  5286.0   4718.0          NTF2-like    ophA
    2   ophA  9536.0  10904.0  Methyltransferase    ophA
    """
    index = _index()
    if name not in index:
        raise ValueError(
            f"There is no dataset {name!r}. The datasets are: " + ", ".join(index) + "."
        )
    types = {c["name"]: c["type"] for c in index[name]["columns"]}
    with (resources.files(__name__) / f"{name}.csv").open(encoding="utf-8") as file:
        text = pd.read_csv(file, dtype="string", keep_default_na=False, na_values="NA")
    columns = {}
    for column, kind in types.items():
        values = text[column].tolist()
        if kind == "numeric":
            # float() gives the nearest number in every pandas version
            values = [math.nan if v is pd.NA else float(v) for v in values]
        elif kind == "integer":
            values = [v if v is pd.NA else int(v) for v in values]
        columns[column] = pd.array(values, dtype=_DTYPES[kind])
    return pd.DataFrame(columns)
