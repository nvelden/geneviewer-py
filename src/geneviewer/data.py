"""The data of a chart: preparation, its columns and clusters, and JSON.

Ported from R: ``prepare_data``, ``chart_fields``, ``chart_state``,
``state_clusters`` and ``cluster_rows`` (``R/spec.R``), ``add_strand`` and
``getUpdatedClusters`` (``R/utils.R``). `jsonable` converts values for the
spec as R's htmlwidgets (jsonlite) writes them.
"""

from __future__ import annotations

import datetime as dt
import math
import warnings
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd

from . import _r


def prepare_data(
    data: pd.DataFrame,
    start: str = "start",
    end: str = "end",
    cluster: str | None = None,
    group: str | None = None,
    strand: str | None = None,
    warn: bool = True,
) -> pd.DataFrame | None:
    """Return the data of a chart as ``GC()`` prepares it.

    Checks the columns, leaves out rows without a start or end, and adds the
    columns ``rowID``, ``start``, ``end`` and ``strand``; cluster names
    become strings. The library does the same with the spec's data; Python
    needs it to check arguments against the data and to process it.

    Parameters
    ----------
    data : pandas.DataFrame
        The data.
    start, end, cluster, group, strand : str, optional
        Column names, as for ``GC()``.
    warn : bool
        Whether to warn about rows left out.

    Returns
    -------
    pandas.DataFrame or None
        A new data frame; ``None`` when no rows are left.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    columns = list(data.columns)
    if start not in columns:
        raise ValueError("start column not found in data")
    if end not in columns:
        raise ValueError("end column not found in data")
    if cluster is not None and cluster not in columns:
        raise ValueError("cluster column not found in data")
    if group is not None and group not in columns:
        raise ValueError("group column not found in data")
    if strand is not None and strand not in columns:
        raise ValueError("strand column not found in data")

    # Rows with NA in the start or end column are left out
    kept = data[start].notna() & data[end].notna()
    data = data[kept].reset_index(drop=True)
    if warn and not kept.all():
        warnings.warn(
            "Rows with NA values in 'start' or 'end' columns have been filtered out.",
            UserWarning,
            stacklevel=3,
        )

    data["rowID"] = np.arange(1, len(data) + 1)
    data_tmp = data.copy()
    data["start"] = data_tmp[start]
    data["end"] = data_tmp[end]
    if cluster is not None:
        data[cluster] = [_r.as_character(value) for value in data[cluster]]
    return add_strand(data, strand)


_FORWARD = ("1", "forward", "+", "sense")
_REVERSE = ("-1", "0", "reverse", "-", "antisense", "complement")


def add_strand(
    data: pd.DataFrame | None, strand: str | None = None
) -> pd.DataFrame | None:
    """Set the strand of each gene, and its start and end to match it.

    Without a strand column (``strand`` or a column named "strand"), a gene
    whose start is before its end is on the forward strand, else on the
    reverse one. With one, its values become "forward" and "reverse", and
    start and end are swapped where they do not match the strand.

    Parameters
    ----------
    data : pandas.DataFrame or None
        Data with ``start`` and ``end`` columns.
    strand : str, optional
        The strand column: 1, "forward", "+" or "sense" for the forward
        strand, -1, 0, "reverse", "-", "antisense" or "complement" for the
        reverse strand.

    Returns
    -------
    pandas.DataFrame or None
        A new data frame with a ``strand`` column; ``None`` for ``None``.
    """
    if data is None:
        return None
    data = data.copy()
    # Data without rows keeps its columns, with an empty strand
    if len(data) == 0:
        data["strand"] = pd.Series([], dtype=object)
        return data
    columns = list(data.columns)
    if strand is None and "strand" in columns:
        strand = "strand"
    if "start" not in columns:
        raise ValueError("start column not found in data")
    if "end" not in columns:
        raise ValueError("end column not found in data")

    if strand is None:
        data["strand"] = np.where(data["start"] < data["end"], "forward", "reverse")
        return data
    if strand not in columns:
        raise ValueError("strand column not found in data")
    values = [_r.as_character(value) for value in data[strand]]
    if not all(value in _FORWARD or value in _REVERSE for value in values):
        raise ValueError(
            "Invalid strand values found. Valid values are 1/forward/+/sense or "
            "-1/0/reverse/-/antisense/complement."
        )
    data["strand"] = ["forward" if value in _FORWARD else "reverse" for value in values]
    swap = ((data["strand"] == "forward") & (data["start"] > data["end"])) | (
        (data["strand"] == "reverse") & (data["start"] < data["end"])
    )
    if swap.any():
        starts = data.loc[swap, "start"].copy()
        data.loc[swap, "start"] = data.loc[swap, "end"]
        data.loc[swap, "end"] = starts
    return data


def chart_fields(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Return the columns and clusters of a chart, from its spec.

    What checking an argument needs, without preparing the data: ``columns``
    and ``clusters`` are those of `prepare_data`. After ``transcript()``,
    the clusters are the transcripts, as in the chart it draws.

    Parameters
    ----------
    spec : Mapping
        The spec of the chart.

    Returns
    -------
    dict
        The column names ``start``, ``end``, ``cluster``, ``group`` and
        ``strand`` (``None`` when not given); ``columns``, the names of the
        columns of the prepared data; and ``clusters``, the names of the
        clusters in order (none when the data has no rows).
    """
    fields = spec.get("fields") or {}
    state: dict[str, Any] = {
        "start": fields.get("start", "start"),
        "end": fields.get("end", "end"),
        "cluster": fields.get("cluster"),
        "group": fields.get("group"),
        "strand": fields.get("strand"),
    }
    transcript = spec.get("transcript")
    if transcript is not None:
        first = transcript if isinstance(transcript, Mapping) else transcript[0]
        state["cluster"] = first.get("transcript", "transcript")
        state["strand"] = first.get("strand")
    data = spec["data"]
    if state["start"] in data.columns and state["end"] in data.columns:
        kept = data[state["start"]].notna() & data[state["end"]].notna()
    else:
        kept = pd.Series(False, index=data.index)
    # Also for data without rows: its columns can still be named
    columns = list(data.columns)
    columns += [c for c in ("rowID", "start", "end", "strand") if c not in columns]
    state["columns"] = columns
    if not kept.any():
        state["clusters"] = []
    elif state["cluster"] is None:
        state["clusters"] = ["cluster"]
    elif state["cluster"] not in data.columns:
        state["clusters"] = []
    else:
        values = [_r.as_character(value) for value in data.loc[kept, state["cluster"]]]
        state["clusters"] = list(dict.fromkeys(values))
    return state


def chart_state(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Return `chart_fields` and the prepared rows of a chart, as ``data``.

    For what works on the rows: links, introns, normalization, alignment
    and editing single rows.

    Parameters
    ----------
    spec : Mapping
        The spec of the chart.

    Returns
    -------
    dict
        The state.
    """
    state = chart_fields(spec)
    state["data"] = prepare_data(
        spec["data"],
        state["start"],
        state["end"],
        state["cluster"],
        state["group"],
        state["strand"],
        warn=False,
    )
    return state


def state_clusters(data: pd.DataFrame, cluster: str | None) -> list[str]:
    """Return the names of the clusters of prepared data.

    Parameters
    ----------
    data : pandas.DataFrame
        Prepared data.
    cluster : str or None
        The cluster column.

    Returns
    -------
    list of str
        The names, in order.
    """
    if len(data) == 0:
        return []
    if cluster is None:
        return ["cluster"]
    return list(dict.fromkeys(_r.as_character(value) for value in data[cluster]))


def cluster_rows(state: Mapping[str, Any], name: str) -> pd.DataFrame:
    """Return the genes of one cluster, as the library splits the data.

    Parameters
    ----------
    state : Mapping
        From `chart_state`.
    name : str
        The cluster.

    Returns
    -------
    pandas.DataFrame
        Its rows, with a ``cluster`` column.
    """
    data = state["data"]
    if state["cluster"] is not None:
        data = data[data[state["cluster"]] == name]
    data = data.copy()
    data["cluster"] = name
    return data


def get_updated_clusters(clusters: list[str], cluster: Any) -> list[str] | None:
    """Return the clusters an argument selects (R's ``getUpdatedClusters``).

    Parameters
    ----------
    clusters : list of str
        The names of the chart's clusters, in order.
    cluster : int, str, list or None
        Cluster numbers (from 1) or names; ``None`` for all.

    Returns
    -------
    list of str or None
        The names; ``None``, with a warning, when ``cluster`` names clusters
        the chart does not have.
    """
    if cluster is None:
        return list(clusters)
    selected = _r.as_vector(cluster)
    if selected and all(_r.is_number(value) for value in selected):
        if any(value > len(clusters) or value < 1 for value in selected):
            warnings.warn(
                "Some cluster numbers provided are out of range. "
                "Please check the cluster numbers.",
                UserWarning,
                stacklevel=3,
            )
            return None
        return [clusters[int(value) - 1] for value in selected]
    names = [_r.as_character(value) for value in selected]
    if not all(name in clusters for name in names):
        warnings.warn(
            "Some cluster names provided are not valid. "
            "Please check the cluster names.",
            UserWarning,
            stacklevel=3,
        )
        return None
    return names


def check_clusters(spec: Mapping[str, Any], cluster: Any) -> None:
    """Warn when ``cluster`` names clusters the chart does not have.

    The library leaves such an entry out.

    Parameters
    ----------
    spec : Mapping
        The spec of the chart.
    cluster : int, str, list or None
        The clusters an argument selects.
    """
    if cluster is not None:
        get_updated_clusters(chart_fields(spec)["clusters"], cluster)


def _scalar(value: Any) -> Any:
    # One value, as jsonlite writes it: NA, NaN and Inf as null, dates as
    # "2020-01-02", times as "2020-01-02T03:04:05Z" (UTC)
    if _r.is_na(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return None if math.isinf(value) else value
    if isinstance(value, (pd.Timestamp, dt.datetime)):
        stamp = pd.Timestamp(value)
        if stamp.tzinfo is not None:
            stamp = stamp.tz_convert("UTC")
        return stamp.strftime("%Y-%m-%dT%H:%M:%SZ")
    if isinstance(value, dt.date):
        return value.strftime("%Y-%m-%d")
    if isinstance(value, np.datetime64):
        return _scalar(pd.Timestamp(value))
    if isinstance(value, np.str_):
        return str(value)
    return value


def _column(series: pd.Series) -> list[Any]:
    if isinstance(series.dtype, pd.CategoricalDtype):
        # A factor: its labels
        return [None if _r.is_na(v) else str(v) for v in series.astype(object)]
    if pd.api.types.is_datetime64_any_dtype(series):
        if all(
            _r.is_na(v) or (v.tzinfo is None and v == v.normalize()) for v in series
        ):
            # Dates without a time, as R's Date
            return [None if _r.is_na(v) else v.strftime("%Y-%m-%d") for v in series]
    return [jsonable(value) for value in series.astype(object)]


def jsonable(value: Any) -> Any:
    """Return a value of the spec as JSON values, as R writes the spec.

    Parameters
    ----------
    value : Any
        A value of the spec: data frames, dicts, lists, numpy and pandas
        values.

    Returns
    -------
    Any
        Dicts, lists, strings, numbers, booleans and ``None``. A data frame
        becomes a dict of columns.
    """
    if isinstance(value, pd.DataFrame):
        return {str(name): _column(value[name]) for name in value.columns}
    if isinstance(value, pd.Series):
        return _column(value)
    if isinstance(value, Mapping):
        return {str(name): jsonable(item) for name, item in value.items()}
    if isinstance(value, (list, tuple, pd.Index)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return [jsonable(item) for item in value.tolist()]
    return _scalar(value)
