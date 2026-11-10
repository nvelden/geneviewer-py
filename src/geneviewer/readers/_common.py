"""What the readers share, ported from R's base functions where they must agree.

The readers return the tables R's readers return: the same columns, in the
same order, with the same types. R's types map to pandas as follows:
character to ``string``, integer to ``Int64``, double to ``float64``,
logical to ``boolean``; ``NA`` is ``<NA>`` (``NaN`` in ``float64``).
"""

from __future__ import annotations

import gzip
import math
import os
import re
import warnings
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

STRING = pd.StringDtype()

# R's type.convert() ----------------------------------------------------------

_LOGICAL = {"T": True, "F": False, "TRUE": True, "FALSE": False}
_INTEGER = re.compile(r"^\s*[+-]?\d+$")
_DOUBLE = re.compile(
    r"^\s*[+-]?("
    r"(\d+\.?\d*|\.\d+)([eE][+-]?\d*)?"
    r"|0[xX]([0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)([pP][+-]?\d+)?"
    r"|[iI][nN][fF]([iI][nN][iI][tT][yY])?"
    r"|[nN][aA][nN]"
    r")\s*$"
)
_INT_MAX = 2147483647


def _is_blank(value: str | None) -> bool:
    return value is None or value.strip() == ""


def _r_double(value: str) -> float:
    # R_strtod: what _DOUBLE accepts
    text = value.strip()
    sign = -1.0 if text.startswith("-") else 1.0
    body = text.lstrip("+-")
    lower = body.lower()
    if lower.startswith("inf"):
        return sign * math.inf
    if lower == "nan":
        return math.nan
    if lower.startswith("0x"):
        return sign * float.fromhex(body if "p" in lower else body + "p0")
    mantissa, _, exponent = lower.partition("e")
    if exponent in ("", "+", "-"):
        return sign * float(mantissa)
    try:
        return sign * float(f"{mantissa}e{exponent}")
    except OverflowError:  # pragma: no cover (float() gives inf instead)
        return sign * math.inf


def type_convert(values: Sequence[str | None]) -> pd.Series:
    """Convert a column of text as R's ``read.table()`` does (``type.convert``).

    ``NA`` is ``None`` or ``"NA"``; in a column that is not text, a field of
    spaces is ``NA`` too. The column is logical when every other value is
    ``T``, ``F``, ``TRUE`` or ``FALSE``; integer when every value is a whole
    number that fits in 32 bits (spaces before it are allowed, after it
    not); double when every value is a number as R reads it (``1e5``,
    ``0x1F``, ``Inf``, ``NaN``, spaces around it); and text otherwise.

    Parameters
    ----------
    values : sequence of str or None
        The fields of one column.

    Returns
    -------
    pandas.Series
        A ``boolean``, ``Int64``, ``float64`` or ``string`` column.
    """
    values = list(values)
    missing = [v is None or v == "NA" for v in values]
    present = [v for v, na in zip(values, missing, strict=True) if not na]
    filled = [v for v in present if not _is_blank(v)]
    if all(v in _LOGICAL for v in filled):
        return pd.Series(
            [
                pd.NA if na or _is_blank(v) else _LOGICAL[v]
                for v, na in zip(values, missing, strict=True)
            ],
            dtype="boolean",
        )
    if all(_INTEGER.match(v) and abs(int(v)) <= _INT_MAX for v in filled):
        return pd.Series(
            [
                pd.NA if na or _is_blank(v) else int(v)
                for v, na in zip(values, missing, strict=True)
            ],
            dtype="Int64",
        )
    if all(_DOUBLE.match(v) for v in filled):
        return pd.Series(
            [
                math.nan if na or _is_blank(v) else _r_double(v)
                for v, na in zip(values, missing, strict=True)
            ],
            dtype="float64",
        )
    return pd.Series(
        [pd.NA if na else v for v, na in zip(values, missing, strict=True)],
        dtype=STRING,
    )


def read_tsv(lines: Sequence[str], col_names: Sequence[str]) -> pd.DataFrame:
    r"""Read tab-separated lines as R's ``read.table()`` with the readers' options.

    The options are ``sep = "\\t"``, ``quote = ""``, ``comment.char = ""``,
    ``fill = TRUE`` and ``header = FALSE``: empty lines are left out, short
    rows are filled with empty fields, and every column goes through
    `type_convert`.

    Parameters
    ----------
    lines : sequence of str
        The lines.
    col_names : sequence of str
        The names of the columns.

    Returns
    -------
    pandas.DataFrame
        The table.

    Raises
    ------
    ValueError
        Without lines, or for a row with more fields than names.
    """
    rows = [line.split("\t") for line in lines if line != ""]
    if not rows:
        raise ValueError("no lines available in input")
    n = len(col_names)
    if any(len(row) > n for row in rows):
        raise ValueError("more columns than column names")
    columns = {
        name: type_convert([row[j] if j < len(row) else "" for row in rows])
        for j, name in enumerate(col_names)
    }
    return pd.DataFrame(columns)


# R's names -------------------------------------------------------------------

_R_RESERVED = {
    "if",
    "else",
    "repeat",
    "while",
    "function",
    "for",
    "next",
    "break",
    "TRUE",
    "FALSE",
    "NULL",
    "Inf",
    "NaN",
    "NA",
    "NA_integer_",
    "NA_real_",
    "NA_character_",
    "NA_complex_",
    "in",
}


def make_unique(names: Iterable[str]) -> list[str]:
    """Make names unique as R's ``make.unique()``: ``a``, ``a.1``, ``a.2``.

    Parameters
    ----------
    names : iterable of str
        The names.

    Returns
    -------
    list of str
        The names, repeated ones with a suffix.
    """
    names = list(names)
    seen = set(names)
    counts: dict[str, int] = {}
    result = []
    used: set[str] = set()
    for name in names:
        if name not in used:
            used.add(name)
            result.append(name)
            continue
        count = counts.get(name, 0)
        while True:
            count += 1
            candidate = f"{name}.{count}"
            if candidate not in seen and candidate not in used:
                break
        counts[name] = count
        used.add(candidate)
        result.append(candidate)
    return result


def make_names(names: Iterable[str | None]) -> list[str]:
    """Make syntactic, unique column names as R's ``data.frame()`` does.

    ``make.names(unique = TRUE)``: characters other than letters, digits,
    ``.`` and ``_`` become ``.``; a name that does not start with a letter,
    or with a dot not followed by a digit, gets an ``X`` in front; reserved
    words get a ``.`` after them (``function.``, ``NA.``).

    Parameters
    ----------
    names : iterable of str or None
        The names; ``None`` is R's ``NA``.

    Returns
    -------
    list of str
        The names.
    """
    result = []
    for name in names:
        if name is None:
            name = "NA"
        name = re.sub(r"[^A-Za-z0-9._]", ".", name)
        if not re.match(r"^([A-Za-z]|\.(?![0-9]))", name):
            name = "X" + name
        if name in _R_RESERVED:
            name += "."
        result.append(name)
    return make_unique(result)


# Metadata --------------------------------------------------------------------


def widen_attributes(
    ids: Sequence[str], keys: Sequence[str | None], values: Sequence[str | None]
) -> tuple[list[str], pd.DataFrame]:
    """Spread key-value metadata into one row per record (R's ``widen_attributes``).

    Records and keys keep the order in which they first appear. A key
    repeated in a record keeps its values joined with ``","`` in input
    order, without the missing ones; when all are missing, the value is
    ``NA``.

    Parameters
    ----------
    ids : sequence of str
        The record of each value.
    keys : sequence of str or None
        The key of each value; ``None`` becomes the column ``"NA"``.
    values : sequence of str or None
        The values.

    Returns
    -------
    tuple
        The record ids, in order, and a table of ``string`` columns, one row
        per record.

    Raises
    ------
    ValueError
        For an empty key, or two keys that give the same column name.
    """
    record_ids = list(dict.fromkeys(ids))
    key_order = list(dict.fromkeys(keys))
    row = {record: i for i, record in enumerate(record_ids)}
    col = {key: j for j, key in enumerate(key_order)}
    cells: dict[tuple[int, int], list[str | None]] = {}
    for record, key, value in zip(ids, keys, values, strict=True):
        cells.setdefault((row[record], col[key]), []).append(value)
    table = np.full((len(record_ids), len(key_order)), None, dtype=object)
    for (i, j), cell in cells.items():
        if len(cell) == 1:
            table[i, j] = cell[0]
        else:
            kept = [v for v in cell if v is not None]
            table[i, j] = ",".join(kept) if kept else None
    names = ["NA" if key is None else key for key in key_order]
    if any(name == "" for name in names):
        raise ValueError("Names can't be empty.")
    if len(set(names)) < len(names):
        raise ValueError("Names must be unique.")
    data = pd.DataFrame(
        {name: pd.array(table[:, j], dtype=STRING) for j, name in enumerate(names)}
    )
    return record_ids, data


# Locations -------------------------------------------------------------------

_SPAN = re.compile(r"[<>]?\d+\.\.[<>]?\d+")


def location_span(location: str | None) -> tuple[float, float]:
    """Return the first and last base of a GenBank location (R's ``location_span``).

    All parts of ``join()`` and ``order()`` count. Parts listed from right to
    left, as in ``join(complement(5441..5600),complement(5095..5373))``, are
    read in order; a location whose parts are not in order (one crossing the
    origin of a circular genome) keeps its first part.

    Parameters
    ----------
    location : str or None
        The location.

    Returns
    -------
    tuple of float
        Start and end; ``NaN`` for a missing or unreadable location.
    """
    if location is None or location is pd.NA:
        return (math.nan, math.nan)
    parts = _SPAN.findall(location)
    if not parts:
        bases = re.findall(r"\d+", location)
        if len(bases) == 1:
            return (float(bases[0]), float(bases[0]))
        return (math.nan, math.nan)
    pos = [tuple(float(p.strip("<>")) for p in part.split("..", 1)) for part in parts]
    starts = [p[0] for p in pos]
    if re.match(r"^\s*(join|order)\(\s*complement", location):
        starts = starts[::-1]
    if any(b < a for a, b in zip(starts, starts[1:], strict=False)):
        return pos[0]
    flat = [v for p in pos for v in p]
    return (min(flat), max(flat))


# Files -----------------------------------------------------------------------


def read_lines(path: str | os.PathLike[str]) -> list[str]:
    r"""Read the lines of a text file as R's ``readLines()``.

    Lines may end in ``\\n``, ``\\r\\n`` or ``\\r``; gzip files are read
    uncompressed.

    Parameters
    ----------
    path : str or path-like
        The file.

    Returns
    -------
    list of str
        The lines, without their line endings.
    """
    with open(path, "rb") as file:
        raw = file.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    text = raw.decode("utf-8", errors="replace")
    if not text:
        # An empty file has no lines (not one empty line)
        return []
    return (
        text.replace("\r\n", "\n")
        .replace("\r", "\n")
        .split("\n")[: -1 if text.endswith(("\n", "\r")) else None]
    )


def file_path_sans_ext(path: str | os.PathLike[str]) -> str:
    """Return a file name without its extension, as R's ``tools::file_path_sans_ext``.

    Parameters
    ----------
    path : str or path-like
        The file.

    Returns
    -------
    str
        The base name without its last extension (``a.fasta.gz`` gives
        ``a.fasta``).
    """
    return re.sub(r"\.[A-Za-z0-9]+$", "", Path(path).name)


def list_files(path: str | os.PathLike[str], pattern: str) -> list[Path]:
    """List the files of a folder whose names match a pattern, sorted by name.

    Parameters
    ----------
    path : str or path-like
        The folder.
    pattern : str
        A regular expression searched in each file name.

    Returns
    -------
    list of pathlib.Path
        The files.
    """
    regex = re.compile(pattern)
    return sorted(
        (p for p in Path(path).iterdir() if regex.search(p.name)),
        key=lambda p: p.name,
    )


def dir_files(path: str | os.PathLike[str], pattern: str, type: str) -> list[Path]:  # noqa: A002
    """List the files of a type in a folder; an error when there are none.

    Parameters
    ----------
    path : str or path-like
        The folder.
    pattern : str
        A regular expression of the file extension.
    type : str
        The file type named in the error.

    Returns
    -------
    list of pathlib.Path
        The files.
    """
    files = list_files(path, pattern)
    if not files:
        raise FileNotFoundError(f"No .{type} files found in the specified directory.")
    return files


def read_tables(
    path: str | os.PathLike[str],
    pattern: str,
    type: str,  # noqa: A002
    read_file: Callable[[Path], pd.DataFrame],
) -> pd.DataFrame:
    """Read a file, or every file of a folder into one table (R's ``read_tables``).

    Each file of a folder gets the column ``filename``: its name without
    the extension.

    Parameters
    ----------
    path : str or path-like
        A file or a folder.
    pattern : str
        A regular expression of the file extension.
    type : str
        The file type named in the error.
    read_file : callable
        Reads one file into a table.

    Returns
    -------
    pandas.DataFrame
        The table.
    """
    path = Path(path)
    if path.is_dir():
        tables = []
        for file in dir_files(path, pattern, type):
            data = read_file(file)
            data["filename"] = pd.array(
                [re.sub(pattern, "", file.name)] * len(data), dtype=STRING
            )
            tables.append(data)
        return bind_rows(tables)
    if path.exists():
        return read_file(path)
    raise FileNotFoundError("The specified path does not exist.")


def _combine(name: str, columns: list[pd.Series]) -> str | np.dtype | Any:
    # The type of a column of bind_rows(): as vctrs combines them
    dtypes = [c.dtype for c in columns if not (c.dtype == "boolean" and c.isna().all())]
    if not dtypes:
        return "boolean"
    kinds = {str(d) for d in dtypes}
    if len(kinds) == 1:
        return dtypes[0]
    if kinds <= {"Int64", "float64"}:
        return "float64"
    raise TypeError(
        f"Can't combine column `{name}` of types {', '.join(sorted(kinds))}."
    )


def bind_rows(tables: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Stack tables as dplyr's ``bind_rows()``.

    Columns are matched by name, in the order they first appear; a column a
    table does not have is ``NA`` there. Integer and double columns combine
    to double; a logical column of only ``NA`` takes the type of the others.

    Parameters
    ----------
    tables : sequence of pandas.DataFrame
        The tables.

    Returns
    -------
    pandas.DataFrame
        One table; empty without tables.
    """
    tables = list(tables)
    if not tables:
        return pd.DataFrame()
    names = list(dict.fromkeys(name for table in tables for name in table.columns))
    result = {}
    for name in names:
        present = [t[name] for t in tables if name in t.columns]
        dtype = _combine(name, present)
        parts = []
        for table in tables:
            if name in table.columns:
                column = table[name]
                if column.dtype == "boolean" and column.isna().all():
                    column = pd.Series([pd.NA] * len(table))
                parts.append(column.astype(object).tolist())
            else:
                parts.append([pd.NA] * len(table))
        values = [v for part in parts for v in part]
        if str(dtype) == "float64":
            values = [math.nan if v is pd.NA or v is None else v for v in values]
        if str(dtype) == "object":
            values = [None if v is pd.NA else v for v in values]
        result[name] = pd.array(values, dtype=dtype)
    return pd.DataFrame(result)


def warn(message: str) -> None:
    """Warn as R's ``warning()`` in a reader.

    Parameters
    ----------
    message : str
        The message.
    """
    warnings.warn(message, UserWarning, stacklevel=3)
