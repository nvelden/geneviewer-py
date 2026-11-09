"""Chart data processing: links, introns and normalize.

``get_introns`` is ported from ``R/transcript.R``. ``get_links``
(``R/chart_links.R``) and the gene normalization of
``R/gene_normalization.R`` follow in M7.6; until then, the chart methods
that need them raise ``NotImplementedError``.
"""

from __future__ import annotations

import re
from typing import Any

import pandas as pd


def _not_yet(what: str) -> NotImplementedError:
    return NotImplementedError(f"{what} is not available yet in this version")


def get_links(state: Any, group: Any, value1: Any, value2: Any, cluster: Any) -> Any:
    """Return links between genes of neighbouring clusters with the same group.

    Parameters
    ----------
    state : Any
        From ``data.chart_state()``.
    group, value1, value2, cluster : Any
        As for ``GC.links()``.

    Returns
    -------
    pandas.DataFrame
        The links.
    """
    raise _not_yet("links(group=...) without data=")


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
    parts = data[is_part].copy()
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


def normalize(spec: Any, **kwargs: Any) -> Any:
    """Return the spec with gene positions normalised.

    Parameters
    ----------
    spec : Any
        The spec of the chart.
    **kwargs
        As for ``GC.normalize()``.

    Returns
    -------
    dict
        The new spec.
    """
    raise _not_yet("normalize()")
