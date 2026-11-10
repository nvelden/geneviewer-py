"""BED reader, ported from R's ``read_bed()`` (``R/read_bed.R``)."""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

from ._common import read_lines, read_tables, read_tsv

FIELDS = (
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
)
PATTERN = r"\.bed[0-9]*$"
# UCSC directives and comments, wherever they appear; words inside
# chromosome or feature names are data, not declarations
_HEADER = re.compile(r"^\s*(track|browser)(\s|$)|^\s*(#|$)")


def read_bed(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read BED files into one table.

    BED files count from 0; the table counts from 1. ``chromStart`` and
    ``thickStart`` are increased by 1; end positions and ``blockSizes`` are
    kept as they are. ``blockStarts`` become 1-based positions relative to
    ``chromStart`` (the first block starts at 1), so the start of a block is
    ``chromStart + blockStart - 1``. ``track`` and ``browser`` lines and
    comments are left out.

    Parameters
    ----------
    path : str or path-like
        A BED file, or a folder of them (``.bed``, ``.bed6``, ``.bed12``, ...).

    Returns
    -------
    pandas.DataFrame
        One row per feature, with as many of the twelve BED columns as the
        first line has (``chrom``, ``chromStart``, ``chromEnd``, ``name``,
        ...). ``blockSizes`` and ``blockStarts`` are lists of numbers. For a
        folder, ``filename``: the file's name without its extension.

    Raises
    ------
    FileNotFoundError
        When the path does not exist, or a folder has no BED files.
    """
    return read_tables(path, PATTERN, "bed", process_bed)


def _blocks(value: object) -> list[float]:
    if value is None or value is pd.NA or (isinstance(value, float) and value != value):
        return [float("nan")]
    text = value if isinstance(value, str) else str(value)
    return [float(v) if v not in ("", "NA") else float("nan") for v in _split(text)]


def _split(text: str) -> list[str]:
    # R's strsplit(x, ","): no empty string after a final comma
    parts = text.split(",")
    return parts[:-1] if parts and parts[-1] == "" else parts


def process_bed(path: Path) -> pd.DataFrame:
    """Read one BED file.

    Parameters
    ----------
    path : pathlib.Path
        The file.

    Returns
    -------
    pandas.DataFrame
        The table of `read_bed`, without ``filename``.
    """
    lines = [line for line in read_lines(path) if not _HEADER.search(line)]
    if lines:
        parts = _split_tabs(lines[0])
        used = FIELDS[: len(parts)]
    else:
        used = FIELDS
    data = read_tsv(lines, used)
    raw = [line.split("\t") for line in lines if line != ""]

    # 0-based starts to 1-based. BED ends are exclusive, so they are already
    # the 1-based end positions.
    for name in ("chromStart", "thickStart"):
        if name in data.columns:
            data[name] = data[name].astype("float64") + 1
    for name, shift in (("blockSizes", 0), ("blockStarts", 1)):
        if name in data.columns:
            j = used.index(name)
            data[name] = pd.Series(
                [
                    [v + shift for v in _blocks(row[j] if j < len(row) else "")]
                    for row in raw
                ],
                dtype=object,
            )
    return data


def _split_tabs(line: str) -> list[str]:
    # R's strsplit(line, "\t"): no empty string after a final tab
    parts = line.split("\t")
    return parts[:-1] if parts and parts[-1] == "" else parts
