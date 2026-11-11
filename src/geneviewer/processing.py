"""Chart data processing: links, introns and normalize.

Ported from R with the same results: ``get_links`` from
``R/chart_links.R``, ``get_introns`` from ``R/transcript.R``, and
``normalize`` from ``GC_normalize()`` (``R/chart_align.R``) with the helpers
of ``R/gene_normalization.R`` and ``update_cluster_data()``
(``R/utils.R``).
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd

from . import _r
from . import data as _data


def _warn(message: str) -> None:
    warnings.warn(message, UserWarning, stacklevel=4)


def _key(value: Any) -> str | None:
    # Values compared as R compares a column with a value: as text
    return None if _r.is_na(value) else _r.as_character(value)


def _keys(series: pd.Series) -> list[str | None]:
    return [_key(v) for v in series]


def _sort_keys(series: pd.Series) -> list[Any]:
    # The order merge() sorts its key by: numbers as numbers, other values as
    # text in the C locale (the order of the fixtures' R)
    if pd.api.types.is_bool_dtype(series) or pd.api.types.is_numeric_dtype(series):
        return [float(v) for v in series]
    return _keys(series)


def _rows(data: pd.DataFrame, keep: list[bool]) -> pd.DataFrame:
    # Rows by a mask; data[[]] would select no columns
    return data.loc[np.asarray(keep, dtype=bool)]


def _suffixed(data: pd.DataFrame, suffix: str) -> pd.DataFrame:
    return data.rename(columns=lambda name: f"{name}{suffix}").reset_index(drop=True)


def _join(
    x: pd.DataFrame, y: pd.DataFrame, pairs: list[tuple[int, int]]
) -> pd.DataFrame:
    xi = [i for i, _ in pairs]
    yi = [j for _, j in pairs]
    return pd.concat(
        [x.iloc[xi].reset_index(drop=True), y.iloc[yi].reset_index(drop=True)],
        axis=1,
    )


def _isort_with_index(values: list[int]) -> list[int]:
    # R's isort_with_index() (src/main/sort.c), a shell sort that merge()
    # uses: tied values do not keep their order
    x = list(values)
    index = list(range(len(x)))
    n = len(x)
    h = 1
    while h <= n // 9:
        h = 3 * h + 1
    while h > 0:
        for i in range(h, n):
            v, iv = x[i], index[i]
            j = i
            while j >= h and x[j - h] > v:
                x[j], index[j] = x[j - h], index[j - h]
                j -= h
            x[j], index[j] = v, iv
        h //= 3
    return index


def _merge_pairs(x_keys: list[Any], y_keys: list[Any]) -> list[tuple[int, int]]:
    # The rows of R's merge(x, y) before it sorts them: the row pairs
    # .Internal(merge()) gives (src/main/unique.c), from 0
    y_set = set(y_keys)
    common = [k for k in x_keys if k in y_set]
    first = {}
    for position, k in enumerate(common, start=1):
        first.setdefault(k, position)
    x_inds = [first.get(k, 0) for k in x_keys]
    y_inds = [first.get(k, 0) for k in y_keys]
    ix, iy = _isort_with_index(x_inds), _isort_with_index(y_inds)
    xs = [x_inds[i] for i in ix]
    ys = [y_inds[i] for i in iy]
    pairs = []
    i = next((n for n, v in enumerate(xs) if v > 0), len(xs))
    j = next((n for n, v in enumerate(ys) if v > 0), len(ys))
    while i < len(xs):
        value = xs[i]
        next_i = i
        while next_i < len(xs) and xs[next_i] == value:
            next_i += 1
        while j < len(ys) and ys[j] < value:
            j += 1
        next_j = j
        while next_j < len(ys) and ys[next_j] == value:
            next_j += 1
        pairs += [(ix[a], iy[b]) for a in range(i, next_i) for b in range(j, next_j)]
        i, j = next_i, next_j
    return pairs


def _merge_on_group(x: pd.DataFrame, y: pd.DataFrame, key: str) -> pd.DataFrame:
    # R's merge(x, y, by.x = key1, by.y = key2): each row of x with each row
    # of y with the same key, sorted by the key (a stable sort, after the
    # unstable one of _merge_pairs); y's key column is left out
    x = _suffixed(x, "1")
    y = _suffixed(y, "2")
    pairs = _merge_pairs(_keys(x[f"{key}1"]), _keys(y[f"{key}2"]))
    sort_by = _sort_keys(x[f"{key}1"])
    pairs.sort(key=lambda pair: sort_by[pair[0]])
    return _join(x, y.drop(columns=f"{key}2"), pairs)


def get_links(
    state: Mapping[str, Any],
    group: Any,
    value1: Any = None,
    value2: Any = None,
    cluster: Any = None,
) -> pd.DataFrame:
    """Return links between genes of neighbouring clusters with the same group.

    Genes link when their value in the ``group`` column is the same, in each
    pair of neighbouring clusters. With ``value1`` and ``value2``, genes of
    ``value1[i]`` link to genes of ``value2[i]`` instead. Rows without a group
    and the group ``"No Hit"`` are left out.

    Parameters
    ----------
    state : Mapping
        From ``data.chart_state()``.
    group : str
        Column of groups.
    value1, value2 : str or list, optional
        Groups to link, pairwise.
    cluster : list, optional
        At least two clusters, by number (from 1) or name. Default: all.

    Returns
    -------
    pandas.DataFrame
        One row per link: the columns of the two genes, with ``1`` and ``2``
        added to their names. Without ``value1`` and ``value2``, the group
        column is there once, as ``<group>1``. The cluster column is named
        ``cluster1`` and ``cluster2``.

    Raises
    ------
    ValueError
        When the group column or values are not in the data, fewer than two
        clusters are given, the chart has no clusters, or no genes link.
    """
    data = state["data"]
    if group is None:
        raise ValueError("Either 'group' or 'data' must be provided.")
    if group not in data.columns:
        raise ValueError(f"Column {group} not found in chart data.")

    keys = _keys(data[group])
    data = _rows(data, [k is not None and k != "No Hit" for k in keys])
    present = set(_keys(data[group]))
    if value1 is not None and not all(_key(v) in present for v in _r.as_vector(value1)):
        raise ValueError("Some values in 'value1' are not present in the group column.")
    if value2 is not None and not all(_key(v) in present for v in _r.as_vector(value2)):
        raise ValueError("Some values in 'value2' are not present in the group column.")

    if cluster is not None and len(_r.as_vector(cluster)) < 2:
        raise ValueError("At least two clusters need to be provided.")
    clusters = _data.get_updated_clusters(state["clusters"], cluster) or []

    cluster_column = state["cluster"]
    if cluster_column is None:
        raise ValueError("Please define cluster in GC().")
    # The links name their clusters in cluster1 and cluster2. Another column
    # named cluster is left out: it would be read in place of the clusters
    # of the chart
    if cluster_column != "cluster" and "cluster" in data.columns:
        if group == "cluster":
            raise ValueError(
                'The column "cluster" cannot be the group of the links when the '
                f'clusters of the chart are another column ("{cluster_column}"). '
                "Rename the column."
            )
        data = data.drop(columns="cluster")
    data = data.rename(columns={cluster_column: "cluster"})
    data = data[data["cluster"].isin(clusters)]

    # R changes value1 and value2 in place, for the cluster pairs that follow
    values1 = None if value1 is None else [_key(v) for v in _r.as_vector(value1)]
    values2 = None if value2 is None else [_key(v) for v in _r.as_vector(value2)]

    links = []
    for name1, name2 in zip(clusters[:-1], clusters[1:], strict=True):
        cluster1 = data[data["cluster"] == name1]
        cluster2 = data[data["cluster"] == name2]
        keys1, keys2 = _keys(cluster1[group]), _keys(cluster2[group])

        if values1 is None or values2 is None:
            common = set(keys1) & set(keys2)
            cluster1 = _rows(cluster1, [k in common for k in keys1])
            cluster2 = _rows(cluster2, [k in common for k in keys2])
            links.append(_merge_on_group(cluster1, cluster2, group))
        elif len(values1) == len(values2):
            for i in range(len(values1)):
                if values1[i] not in keys1:
                    values1[i], values2[i] = values2[i], values1[i]
                rows1 = _rows(cluster1, [k == values1[i] for k in keys1])
                rows2 = _rows(cluster2, [k == values2[i] for k in keys2])
                if len(rows1) == 0 or len(rows2) == 0:
                    continue
                pairs = [(a, b) for a in range(len(rows1)) for b in range(len(rows2))]
                links.append(_join(_suffixed(rows1, "1"), _suffixed(rows2, "2"), pairs))

    links = [frame for frame in links if len(frame) > 0]
    if not links:
        raise ValueError("No links found with the specified parameters.")
    return pd.concat(links, ignore_index=True)


def get_introns(data: pd.DataFrame) -> pd.DataFrame:
    """Add the introns between the exons of each transcript.

    Exons and UTRs that overlap, nest or touch are one covered region; each
    gap of at least one base between covered regions becomes an intron. An
    intron row copies the other columns of the exon or UTR before it. Exons
    whose start is after their end (reverse strand) are read as their
    span, and keep their positions in the result.

    Parameters
    ----------
    data : pandas.DataFrame
        Exons and UTRs, with ``start`` and ``end``. Optional columns:
        ``type`` (rows whose type contains ``exon`` or ``utr``, in any case,
        are used; default ``"exon"``) and ``transcript`` (default
        ``"transcript1"``).

    Returns
    -------
    pandas.DataFrame
        The rows of `data` and the introns (``type`` ``"intron"``), by
        transcript in the order they first appear, then by ``start``.
        ``start`` and ``end`` are floats, as R's numbers.

    Raises
    ------
    ValueError
        Without ``start`` or ``end``.

    Examples
    --------
    >>> exons = pd.DataFrame({"transcript": "tx1", "type": "exon",
    ...                       "start": [100, 400, 900], "end": [250, 600, 1100]})
    >>> get_introns(exons)[["start", "end", "type"]]
        start     end    type
    0   100.0   250.0    exon
    1   251.0   399.0  intron
    2   400.0   600.0    exon
    3   601.0   899.0  intron
    4   900.0  1100.0    exon
    """
    data = data.copy().reset_index(drop=True)
    if "type" not in data.columns:
        data["type"] = pd.array(["exon"] * len(data), dtype="string")
    if "transcript" not in data.columns:
        data["transcript"] = pd.array(["transcript1"] * len(data), dtype="string")
    if "start" not in data.columns or "end" not in data.columns:
        raise ValueError("The 'start' and 'end' columns are required in the data.")

    order = {t: i for i, t in enumerate(pd.unique(data["transcript"]))}
    is_part = [
        not pd.isna(t) and re.search("exon|utr", str(t), re.IGNORECASE) is not None
        for t in data["type"]
    ]
    # A boolean Series: an empty list would select columns, not rows
    parts = data[pd.Series(is_part, index=data.index, dtype=bool)].copy()
    start = parts["start"].astype("float64")
    end = parts["end"].astype("float64")
    swap = start > end
    parts["start"], parts["end"] = start.where(~swap, end), end.where(~swap, start)

    # The running end merges overlapping and nested parts; touching parts have
    # no intron between them
    parts = parts.sort_values(["transcript", "start"], kind="stable")
    by = parts.groupby("transcript", dropna=False, sort=False)
    covered = by["end"].cummax()
    following = by["start"].shift(-1)
    gap = following.notna() & (following > covered + 1)
    introns = parts[gap].copy()
    introns["start"] = covered[gap] + 1
    introns["end"] = following[gap] - 1
    introns["type"] = "intron"

    result = pd.concat([data, introns], ignore_index=True)
    for name in ("start", "end"):
        result[name] = result[name].astype("float64")
    rank = result["transcript"].map(order)
    result = (
        result.assign(_order=rank)
        .sort_values(["_order", "start"], kind="stable", na_position="last")
        .drop(columns="_order")
        .reset_index(drop=True)
    )
    return result


def _recycle(value: Any, i: int) -> Any:
    # R's recycle_at(): the i-th value (from 0), recycled; None stays None
    if value is None:
        return None
    values = _r.as_vector(value)
    return values[i % len(values)]


def _order_genes(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy().reset_index(drop=True)
    data["original_order"] = np.arange(1, len(data) + 1)
    start = data["start"].astype("float64")
    end = data["end"].astype("float64")
    data["is_reverse"] = end < start
    data["actual_start"] = np.minimum(start, end)
    data["actual_end"] = np.maximum(start, end)
    return data.sort_values(
        "actual_start", kind="stable", na_position="last"
    ).reset_index(drop=True)


def _restore_gene_order(data: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    data = data.sort_values("original_order", kind="stable")
    return data.drop(columns=columns).reset_index(drop=True)


def _new_positions(data: pd.DataFrame, new_start: Any, new_end: Any) -> pd.DataFrame:
    reverse = data["is_reverse"].to_numpy(dtype=bool)
    new_start = np.asarray(new_start, dtype="float64")
    new_end = np.asarray(new_end, dtype="float64")
    data["start"] = np.where(reverse, new_end, new_start)
    data["end"] = np.where(reverse, new_start, new_end)
    return data


def _custom_gaps(
    data: pd.DataFrame, group: str, gene_name: str, factor: float = 2
) -> pd.DataFrame:
    """Scale the gap after each gene of a group.

    The genes after it move by ``(factor - 1)`` times the gap. The last gene
    has no gap after it. ``original_start`` and ``original_end`` are added
    when the data has none yet.

    Parameters
    ----------
    data : pandas.DataFrame
        The genes of one cluster, with ``start`` and ``end``.
    group : str
        Column of groups.
    gene_name : str
        The group whose gaps change.
    factor : float
        The factor. Default 2.

    Returns
    -------
    pandas.DataFrame
        The genes, in the order given.
    """
    if gene_name not in _keys(data[group]):
        _warn(f"{gene_name} not found in the specified group column.")
        return data
    data = _order_genes(data)
    n = len(data)
    index = [i for i, k in enumerate(_keys(data[group])) if k == gene_name]
    index = [i for i in index if i < n - 1]

    actual_start = data["actual_start"].to_numpy()
    actual_end = data["actual_end"].to_numpy()
    shift = np.zeros(n)
    for i in index:
        shift[i + 1 :] += (factor - 1) * (actual_start[i + 1] - actual_end[i])

    if "original_start" not in data.columns:
        data["original_start"] = data["start"]
    if "original_end" not in data.columns:
        data["original_end"] = data["end"]
    data = _new_positions(data, actual_start + shift, actual_end + shift)
    return _restore_gene_order(
        data, ["actual_start", "actual_end", "original_order", "is_reverse"]
    )


def _normalize_positions(
    data: pd.DataFrame, gap: float | None = None, preserve_gene_length: bool = True
) -> pd.DataFrame:
    """Space the genes of one cluster evenly.

    The genes keep their order, strand and the span of the cluster. The total
    of the gaps is kept when ``gap`` is ``None``; otherwise each gap is
    ``gap`` times the span. Positions are rounded to whole numbers (halves to
    even, as R rounds). The old positions go to ``original_start`` and
    ``original_end``.

    Parameters
    ----------
    data : pandas.DataFrame
        The genes, with ``start`` and ``end``.
    gap : float, optional
        The gap, as a share of the span.
    preserve_gene_length : bool
        Whether genes keep their length; otherwise all get the mean length.

    Returns
    -------
    pandas.DataFrame
        The genes, in the order given.
    """
    if "start" not in data.columns or "end" not in data.columns:
        _warn("The data does not contain 'start' or 'end' columns.")
        return data
    data = _order_genes(data)
    data["original_start"] = data["start"]
    data["original_end"] = data["end"]

    actual_start = data["actual_start"].to_numpy()
    actual_end = data["actual_end"].to_numpy()
    total_length = np.nanmax(actual_end) - np.nanmin(actual_start)
    n = len(data)
    if gap is None:
        total_gap = np.nansum(actual_start[1:] - actual_end[:-1])
    else:
        total_gap = total_length * gap * (n - 1)
    mean_gene_length = (total_length - total_gap) / n
    # The gaps lie between the genes: one fewer than there are genes
    step = total_gap / (n - 1) if n > 1 else 0.0

    if preserve_gene_length:
        gene_length = np.abs(actual_end - actual_start)
        offsets = np.concatenate(([0.0], np.cumsum(gene_length + step)[:-1]))
    else:
        gene_length = np.full(n, mean_gene_length)
        offsets = np.arange(n) * (mean_gene_length + step)
    new_start = np.nanmin(actual_start) + offsets
    new_end = new_start + gene_length
    data = _new_positions(data, np.round(new_start), np.round(new_end))
    return _restore_gene_order(
        data, ["actual_start", "actual_end", "original_order", "is_reverse"]
    )


def _update_cluster_data(
    data: pd.DataFrame, updated: pd.DataFrame, by: str = "rowID"
) -> pd.DataFrame:
    """Put the rows of one cluster back into the data, by ``rowID``.

    R's ``update_cluster_data()``: values of ``updated`` replace those of
    ``data`` where they are not missing; its other columns are added, missing
    for the other rows. The rows are sorted by ``rowID``.

    Parameters
    ----------
    data : pandas.DataFrame
        All rows.
    updated : pandas.DataFrame
        The changed rows.
    by : str
        The column that identifies a row.

    Returns
    -------
    pandas.DataFrame
        The data.
    """
    if by not in data.columns:
        raise ValueError(f"The column {by} does not exist in the 'data' data frame.")
    if by not in updated.columns:
        raise ValueError(
            f"The column {by} does not exist in the 'updated_data' data frame."
        )
    result = data.merge(updated, on=by, how="left", suffixes=("", "_new"))
    for column in [c for c in updated.columns if c != by and c in data.columns]:
        new = result.pop(f"{column}_new")
        old = result[column]
        if new.dtype != old.dtype and not (
            pd.api.types.is_numeric_dtype(new) and pd.api.types.is_numeric_dtype(old)
        ):
            old = old.astype(object)
            new = new.astype(object)
        result[column] = new.where(new.notna(), old)
    return result.sort_values(by, kind="stable").reset_index(drop=True)


def normalize(
    spec: Mapping[str, Any],
    group: str | None = None,
    cluster: Any = None,
    normalize: Any = True,
    preserve_gene_length: Any = True,
    gap: Any = None,
    custom_gaps: Mapping[str, float] | None = None,
) -> Any:
    """Return the spec with the gene positions normalised (``GC_normalize``).

    Parameters
    ----------
    spec : Mapping
        The spec of the chart.
    group, cluster, normalize, preserve_gene_length, gap, custom_gaps
        As for ``GC.normalize()``. ``preserve_gene_length`` and ``gap`` may
        have one value per cluster, recycled.

    Returns
    -------
    dict
        The new spec. Its data is the chart's prepared data (with ``rowID``
        and ``strand``) and the new positions, also in the columns given as
        start and end.
    """
    state = _data.chart_state(spec)
    if state["cluster"] is None and cluster is not None:
        _warn("Could not normalize selected genes. Please define cluster in GC().")
        return spec
    if group is None and state["group"] is not None:
        group = state["group"]
    if group is not None and len(_r.as_vector(group)) != 1:
        raise ValueError("group must be a single column name")
    if group is not None and group not in state["data"].columns:
        raise ValueError("group column not found in data")

    clusters = _data.get_updated_clusters(state["clusters"], cluster) or []
    gaps = dict(custom_gaps or {})
    data = state["data"]
    for i, name in enumerate(clusters):
        cluster_data = _data.cluster_rows(state, name)
        if normalize:
            cluster_data = _normalize_positions(
                cluster_data,
                gap=_recycle(gap, i),
                preserve_gene_length=bool(_recycle(preserve_gene_length, i)),
            )
        if gaps:
            if group is None:
                _warn(
                    "group is None. Custom gaps will not be applied because the "
                    "group column is not specified."
                )
            else:
                for gene_name, factor in gaps.items():
                    cluster_data = _custom_gaps(cluster_data, group, gene_name, factor)
        data = _update_cluster_data(data, cluster_data)

    # The spec holds the shifted positions, also in the columns given as start
    # and end
    for field in ("start", "end"):
        data[state[field]] = data[field]
    result = {**spec, "data": data}
    if spec.get("links") is not None:
        fields = list(dict.fromkeys(["start", "end", state["start"], state["end"]]))
        result["links"] = _normalize_links(spec["links"], data, fields)
    return result


def _normalize_links(links: Any, data: pd.DataFrame, fields: list[str]) -> Any:
    """Put links made before ``normalize()`` at the new positions of their genes.

    R's ``normalize_links()``: links made by ``links(group=)`` name their
    genes in ``rowID1`` and ``rowID2``, and their positions are those of the
    genes. Links given as data, e.g. alignments, have positions inside the
    genes; they keep them, with a warning.

    Parameters
    ----------
    links : dict or list of dict
        The links of the spec: an entry or a list of entries.
    data : pandas.DataFrame
        The chart's data with the new positions and ``rowID``.
    fields : list of str
        The position columns: ``start``, ``end`` and the columns given as
        start and end.

    Returns
    -------
    dict or list of dict
        The links.
    """
    single = isinstance(links, Mapping)
    entries = [links] if single else list(links)
    kept = False
    position = {row: i for i, row in enumerate(data["rowID"].tolist())}
    for i, entry in enumerate(entries):
        links_data = entry.get("data")
        if not isinstance(links_data, pd.DataFrame) or not {"rowID1", "rowID2"} <= set(
            links_data.columns
        ):
            kept = True
            continue
        links_data = links_data.copy()
        for side in ("1", "2"):
            rows = [position.get(row) for row in links_data["rowID" + side].tolist()]
            for field in fields:
                column = field + side
                if column not in links_data.columns:
                    continue
                values = data[field].tolist()
                old = links_data[column].tolist()
                new = [old[j] if r is None else values[r] for j, r in enumerate(rows)]
                links_data[column] = pd.Series(new, index=links_data.index)
        entries[i] = {**entry, "data": links_data}
    if kept:
        _warn(
            "Links given as data keep their positions after normalize(). "
            "Add them after normalize() to place them on the new positions."
        )
    return entries[0] if single else entries
