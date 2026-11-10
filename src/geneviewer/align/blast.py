"""Protein alignment between clusters: ``protein_blast`` (``R/cluster_blast.R``).

R aligns with ``pwalign::pairwiseAlignment()``; Python aligns with Biopython's
``PairwiseAligner`` with the same scores (see `aligner`). Where several
alignments have the best score, the two packages can pick different ones, so
identity and similarity can differ a little from R's (``docs/differences.md``).
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from .. import _r
from ..processing import _merge_pairs
from ..readers.genbank import read_gbk_cds

# pwalign's defaults for two character vectors
GAP_OPENING = 10
GAP_EXTENSION = 4
QUALITY = 22

_SIMILAR_ACIDS = ("GAVLI", "FYW", "CM", "ST", "KRH", "DENQ", "P")
_ACID_GROUP = {acid: k for k, group in enumerate(_SIMILAR_ACIDS) for acid in group}


def _key(value: Any) -> str | None:
    return None if _r.is_na(value) else _r.as_character(value)


def substitution_scores(alphabet_length: int) -> tuple[float, float]:
    """Return the match and mismatch scores pwalign uses for character input.

    ``pwalign::pairwiseAlignment()`` on two character vectors, with its
    defaults (no ``substitutionMatrix``), scores by base quality: every
    letter has Phred quality 22 (error probability ``e = 10^-2.2``), and the
    alphabet is the letters of all patterns and subjects together (``n``).
    With ``E = 2e - n/(n-1) e^2``, a match scores ``log2((1 - E) n)`` and a
    mismatch ``log2(E n / (n-1))`` (``qualitySubstitutionMatrices()``).

    Parameters
    ----------
    alphabet_length : int
        The number of different letters.

    Returns
    -------
    tuple of float
        The match and mismatch scores.
    """
    n = max(alphabet_length, 2)
    e = 10 ** (-QUALITY / 10)
    combined = 2 * e - n / (n - 1) * e * e
    return math.log2((1 - combined) * n), math.log2(combined * n / (n - 1))


def aligner(alphabet_length: int) -> Any:
    """Return a Biopython aligner scoring as pwalign's defaults.

    A global alignment (end gaps count), the quality-based scores of
    `substitution_scores`, and gaps as pwalign scores them: a gap of ``k``
    letters costs ``10 + 4k`` (``gapOpening = 10``, ``gapExtension = 4``).

    Parameters
    ----------
    alphabet_length : int
        The number of different letters of all sequences aligned.

    Returns
    -------
    Bio.Align.PairwiseAligner
        The aligner.
    """
    from Bio.Align import PairwiseAligner

    match, mismatch = substitution_scores(alphabet_length)
    return PairwiseAligner(
        mode="global",
        match_score=match,
        mismatch_score=mismatch,
        open_gap_score=-(GAP_OPENING + GAP_EXTENSION),
        extend_gap_score=-GAP_EXTENSION,
    )


def _trim(pattern: str, subject: str) -> tuple[str, str]:
    # pwalign's pattern() and subject() leave out the end gaps: the columns
    # before the first and after the last column with a letter in both
    both = [
        i
        for i, (a, b) in enumerate(zip(pattern, subject, strict=True))
        if a != "-" and b != "-"
    ]
    if not both:
        return "", ""
    return pattern[both[0] : both[-1] + 1], subject[both[0] : both[-1] + 1]


def align_pairs(
    patterns: Sequence[str], subjects: Sequence[str]
) -> list[tuple[str, str]]:
    """Align each pattern with its subject as R's ``protein_blast()`` does.

    The sequences are aligned from their ends, as pwalign fills its matrices,
    which makes Biopython pick pwalign's alignment more often where several
    have the best score.

    Parameters
    ----------
    patterns, subjects : sequence of str
        The sequences, pairwise.

    Returns
    -------
    list of tuple of str
        The aligned pattern and subject of each pair, with gaps ``-``, without
        the end gaps (as ``pwalign::pattern()`` and ``subject()`` give them).
    """
    if len(patterns) != len(subjects):
        raise ValueError(
            "The lists of patterns and subjects must be of the same length"
        )
    if not patterns:
        return []
    letters = set("".join(patterns)) | set("".join(subjects))
    al = aligner(len(letters))
    done: dict[tuple[str, str], tuple[str, str]] = {}
    result = []
    for pattern, subject in zip(patterns, subjects, strict=True):
        pair = (pattern, subject)
        if pair not in done:
            if not pattern or not subject:
                done[pair] = ("", "")
            else:
                best = al.align(pattern[::-1], subject[::-1])[0]
                done[pair] = _trim(best[0][::-1], best[1][::-1])
        result.append(done[pair])
    return result


def compute_identity(patterns: Sequence[str], subjects: Sequence[str]) -> pd.DataFrame:
    """Return the identity and similarity of aligned pairs (R's ``compute_identity``).

    Positions where both have a gap do not count. Similar amino acids are in
    the same group of G A V L I, F Y W, C M, S T, K R H, D E N Q and P.

    Parameters
    ----------
    patterns, subjects : sequence of str
        Aligned sequences of equal length, pairwise.

    Returns
    -------
    pandas.DataFrame
        ``identity`` and ``similarity``, in percent, rounded to one decimal
        (``round(x, 3) * 100``, as R).
    """
    if len(patterns) != len(subjects):
        raise ValueError(
            "The lists of patterns and subjects must be of the same length"
        )
    if any(len(p) != len(s) for p, s in zip(patterns, subjects, strict=True)):
        raise ValueError("The two sequences must be of the same length")
    identity, similarity = [], []
    for pattern, subject in zip(patterns, subjects, strict=True):
        length = matches = similar = 0
        for a, b in zip(pattern, subject, strict=True):
            if a == b and a in "-.":
                continue
            length += 1
            if a == b:
                matches += 1
            elif (
                a in _ACID_GROUP
                and b in _ACID_GROUP
                and _ACID_GROUP[a] == _ACID_GROUP[b]
            ):
                similar += 1
        if length == 0:
            identity.append(math.nan)
            similarity.append(math.nan)
        else:
            identity.append(round(matches / length, 3) * 100)
            similarity.append(round((matches + similar) / length, 3) * 100)
    return pd.DataFrame({"identity": identity, "similarity": similarity})


def synteny_score(
    order1: Sequence[float],
    order2: Sequence[float],
    identity: Sequence[float],
    i: float = 0.5,
) -> float:
    """Return the synteny score of a cluster (R's ``synteny_score``).

    One point for each pair of neighbouring genes whose hits are neighbours
    too (pairs do not overlap), plus ``i`` times the summed identity / 100.

    Parameters
    ----------
    order1 : sequence of float
        The places of the query genes, sorted.
    order2 : sequence of float
        The places of their hits in the cluster.
    identity : sequence of float
        The identities.
    i : float
        The weight of the identity.

    Returns
    -------
    float
        The score.
    """
    score = 0.0
    j = 0
    n = len(order1)
    while j < n - 1:
        if abs(order1[j + 1] - order1[j]) < 2 and abs(order2[j + 1] - order2[j]) < 2:
            score += 1
            j += 2
        else:
            j += 1
    return score + (sum(identity) / 100) * i


def _merge_by_row(
    x: pd.DataFrame, x_key: str, y: pd.DataFrame, suffixes: tuple[str, str]
) -> pd.DataFrame:
    # R's merge(x, y, by.x = x_key, by.y = "rowID", all.x = TRUE) when every
    # row of x has its row in y: R's row order (_merge_pairs, then a stable
    # sort by the key); the key first, then x's columns, then y's
    pairs = _merge_pairs([_key(v) for v in x[x_key]], [_key(v) for v in y["rowID"]])
    keys = x[x_key].to_numpy()
    pairs.sort(key=lambda pair: keys[pair[0]])
    xi = [a for a, _ in pairs]
    yi = [b for _, b in pairs]
    x_rows = x.iloc[xi].reset_index(drop=True)
    y_rows = y.drop(columns="rowID").iloc[yi].reset_index(drop=True)
    common = set(x_rows.columns) & set(y_rows.columns)
    x_rows = x_rows.rename(columns={c: c + suffixes[0] for c in common})
    y_rows = y_rows.rename(columns={c: c + suffixes[1] for c in common})
    columns = [x_key] + [c for c in x_rows.columns if c != x_key]
    return pd.concat([x_rows[columns], y_rows], axis=1)


def _protein_combinations(
    data: pd.DataFrame, pair: tuple[str, str], cluster: str
) -> pd.DataFrame:
    # R's get_protein_combinations(): every gene of the first cluster with
    # every gene of the second, with their sequences and places
    sides = []
    for name in pair:
        rows = data[data[cluster] == name]
        first = np.fmin(
            rows["start"].to_numpy(dtype="float64"),
            rows["end"].to_numpy(dtype="float64"),
        )
        rows = rows.iloc[np.argsort(first, kind="stable")].reset_index(drop=True)
        rows = rows[["rowID", "translation"]].assign(order=np.arange(1, len(rows) + 1))
        sides.append(rows)
    first, second = sides
    grid = pd.DataFrame(
        {
            "rowID.x": np.tile(first["rowID"].to_numpy(), len(second)),
            "rowID.y": np.repeat(second["rowID"].to_numpy(), len(first)),
        }
    )
    combinations = _merge_by_row(grid, "rowID.x", first, ("1", "2"))
    combinations = _merge_by_row(combinations, "rowID.y", second, ("1", "2"))
    combinations["cluster1"] = pair[0]
    combinations["cluster2"] = pair[1]
    return combinations


def _left_join(
    x: pd.DataFrame, y: pd.DataFrame, x_key: str, y_key: str
) -> pd.DataFrame:
    # dplyr's left_join(): x's rows in order, each with its rows of y in
    # order (or missing values); names both have get ".x" and ".y"
    positions: dict[Any, list[int]] = {}
    for j, value in enumerate(y[y_key]):
        positions.setdefault(_key(value), []).append(j)
    xi, yi = [], []
    for i, value in enumerate(x[x_key]):
        for j in positions.get(_key(value), [None]):
            xi.append(i)
            yi.append(j)
    x_rows = x.iloc[xi].reset_index(drop=True)
    y_other = y.drop(columns=y_key)
    matched = [j for j in yi if j is not None]
    y_rows = y_other.iloc[matched].reset_index(drop=True)
    y_rows.index = [k for k, j in enumerate(yi) if j is not None]
    y_rows = y_rows.reindex(range(len(yi)))
    common = set(x_rows.columns) & set(y_rows.columns)
    x_rows = x_rows.rename(columns={c: f"{c}.x" for c in common})
    y_rows = y_rows.rename(columns={c: f"{c}.y" for c in common})
    return pd.concat([x_rows, y_rows], axis=1)


def protein_blast(
    data: pd.DataFrame | str | os.PathLike[str],
    query: str,
    id: str = "protein_id",  # noqa: A002
    start: str = "start",
    end: str = "end",
    cluster: str = "cluster",
    genes: Sequence[str] | None = None,
    identity: float = 30,
) -> pd.DataFrame:
    """Find the best hit of each gene in the query cluster's proteins.

    Every protein is aligned globally with every protein of the query
    cluster (Needleman-Wunsch, as R's ``protein_blast()`` with
    ``pwalign::pairwiseAlignment()``; see `aligner` for the scores). A gene
    keeps the query protein with the highest identity as its hit, when the
    identity is at least ``identity``; equal best hits give a row each. The
    clusters are then ordered by a synteny score, the query first.

    Needs Biopython: ``pip install "geneviewer-py[align]"``.

    Parameters
    ----------
    data : pandas.DataFrame or path
        Genes with their protein sequence (column ``translation`` or
        ``sequence``), or a GenBank file or folder of them (one cluster per
        file).
    query : str
        The query cluster.
    id : str
        Column of unique protein ids.
    start, end : str
        Columns of the start and end positions.
    cluster : str
        Column of cluster names.
    genes : list of str, optional
        Only these genes of the query cluster (by id). Default: all.
    identity : float
        The lowest identity (percent) of a hit. Default 30.

    Returns
    -------
    pandas.DataFrame
        `data` (rows without an id, sequence or cluster left out) with
        ``rowID``, ``start``, ``end``, and ``identity``, ``similarity``,
        ``BlastP`` (the id of the hit, ``"No Hit"`` without one; a hit only
        one gene has is ``"No Hit"`` too) and ``score`` (the cluster's
        synteny score). Ordered by cluster, the query first, then by score.

    Raises
    ------
    ImportError
        Without Biopython.
    ValueError
        When a column or the query is not in the data.

    Examples
    --------
    >>> genes = pd.DataFrame({
    ...     "protein_id": ["A1", "A2", "B1", "B2"],
    ...     "cluster": ["A", "A", "B", "B"],
    ...     "start": [1, 200, 1, 200], "end": [150, 350, 150, 350],
    ...     "translation": ["MKTAYIAKQRQISFVKSHFSRQ", "MSDNELLKRLGEAVLAHPTTG",
    ...                     "MKTAYIAKQRQISFVKSHFSRQ", "MLPGGRTSWVQEHCDKAYNIF"]})
    >>> protein_blast(genes, query="A")[["protein_id", "identity", "BlastP"]]
      protein_id  identity  BlastP
    0         A1     100.0      A1
    1         A2     100.0  No Hit
    2         B1     100.0      A1
    3         B2       NaN  No Hit
    """
    try:
        import Bio.Align  # noqa: F401
    except ImportError:
        raise ImportError(
            "Biopython is not installed. Install with "
            'pip install "geneviewer-py[align]".'
        ) from None

    if isinstance(data, (str, os.PathLike)):
        data = read_gbk_cds(data, start, end)
        cluster = "cluster"
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame, or a GenBank file or folder")
    columns = list(data.columns)
    if id not in columns:
        raise ValueError(
            "Please use the id variable to specificy the column name containing "
            "unique protein identifiers."
        )
    if start not in columns:
        raise ValueError(
            "Please use the start variable to specificy the column name specifying "
            "the start positions."
        )
    if end not in columns:
        raise ValueError(
            "Please use the end variable to specificy the column name specifying "
            "the end positions."
        )
    if cluster not in columns:
        raise ValueError(
            "Please use the cluster variable to specificy the column name specifying "
            "the cluster names."
        )
    query = _key(query)
    if query not in {_key(v) for v in data[cluster]}:
        raise ValueError("Query not found in cluster columnm.")
    if "sequence" not in columns and "translation" not in columns:
        raise ValueError("No 'sequence' or 'translation' column found in data.")

    data = data.copy()
    if "sequence" in columns and "translation" not in columns:
        data = data.rename(columns={"sequence": "translation"})
    keep = data[id].notna() & data["translation"].notna() & data[cluster].notna()
    data = data[keep.to_numpy(dtype=bool)].reset_index(drop=True)
    data["rowID"] = pd.array(np.arange(1, len(data) + 1), dtype="Int64")
    data_tmp = data.copy()
    data["start"] = data_tmp[start]
    data["end"] = data_tmp[end]
    data[cluster] = pd.array([_key(v) for v in data[cluster]], dtype="string")

    clusters = list(pd.unique(data[cluster]))
    combination_data = data
    if genes is not None:
        wanted = {_key(g) for g in _r.as_vector(genes)}
        ids = [_key(v) for v in data[id]]
        query_genes = {
            k
            for k, c in zip(ids, data[cluster], strict=True)
            if c == query and k not in wanted
        }
        combination_data = data[[k not in query_genes for k in ids]]

    combinations = pd.concat(
        [
            _protein_combinations(combination_data, (query, c), cluster)
            for c in clusters
        ],
        ignore_index=True,
    )
    to_align = combinations[combinations["cluster1"] != combinations["cluster2"]]
    to_align = to_align.reset_index(drop=True)
    query_rows = combinations[combinations["rowID.x"] == combinations["rowID.y"]]
    query_rows = query_rows.reset_index(drop=True)

    # With only the query cluster there is nothing to align (R's B67)
    aligned = align_pairs(
        list(to_align["translation1"]), list(to_align["translation2"])
    )
    scores = compute_identity([p for p, _ in aligned], [s for _, s in aligned])
    to_align["identity"] = scores["identity"].to_numpy()
    to_align["similarity"] = scores["similarity"].to_numpy()
    to_align = to_align[(to_align["identity"] >= identity).to_numpy()]
    query_rows = query_rows.assign(identity=100.0, similarity=100.0)

    hits = pd.concat([to_align, query_rows], ignore_index=True)[
        ["rowID.y", "rowID.x", "order1", "order2", "identity", "similarity"]
    ]
    # The best hits of each gene, by rowID.y (dplyr's slice_max keeps ties)
    hits = hits.iloc[np.argsort(hits["rowID.y"].to_numpy(), kind="stable")]
    best = hits.groupby("rowID.y", sort=False)["identity"].transform("max")
    hits = hits[(hits["identity"] == best).to_numpy()].reset_index(drop=True)

    ids_by_row = dict(zip(data["rowID"], data[id], strict=True))
    blastp = [ids_by_row[r] for r in hits["rowID.x"]]
    counts = pd.Series([_key(v) for v in blastp]).value_counts()
    # A hit only one gene has is not shown as a group
    hits["BlastP"] = pd.Series(
        [None if counts[_key(v)] == 1 else v for v in blastp], dtype=object
    )

    data = _left_join(data, hits, "rowID", "rowID.y")

    # The synteny score of each cluster
    scored = data[data["order1"].notna().to_numpy()]
    scored = scored.iloc[
        np.argsort(scored["order1"].to_numpy(dtype="float64"), kind="stable")
    ]
    synteny = {
        name: synteny_score(
            list(rows["order1"].astype("float64")),
            list(rows["order2"].astype("float64")),
            list(rows["identity"].astype("float64")),
        )
        for name, rows in scored.groupby(cluster, sort=False)
    }
    data["score"] = [synteny.get(name, math.nan) for name in data[cluster]]
    score = data["score"].to_numpy(dtype="float64")
    # dplyr's arrange(desc(score)), then order(cluster != query, -score):
    # both stable, missing scores last
    by_score = np.argsort(np.where(np.isnan(score), np.inf, -score), kind="stable")
    data = data.iloc[by_score].reset_index(drop=True)
    score = data["score"].to_numpy(dtype="float64")
    not_query = (data[cluster] != query).to_numpy(dtype=bool)
    order = np.lexsort((np.where(np.isnan(score), np.inf, -score), not_query))
    data = data.iloc[order].reset_index(drop=True)

    data["rowID"] = pd.array(np.arange(1, len(data) + 1), dtype="Int64")
    # The ids of the hits; with "No Hit", numeric ids become text, as in R
    missing = data["BlastP"].isna()
    if pd.api.types.is_numeric_dtype(data[id]) and not missing.any():
        data["BlastP"] = data["BlastP"].astype(data[id].dtype)
    else:
        data["BlastP"] = pd.array(
            [
                "No Hit" if m else _key(v)
                for v, m in zip(data["BlastP"], missing, strict=True)
            ],
            dtype="string",
        )
    return data.drop(columns=["rowID.x", "order1", "order2"])
