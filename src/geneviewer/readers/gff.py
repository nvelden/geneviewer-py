"""GFF reader, ported from R's ``read_gff()`` (``R/read_gff.R``)."""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from pathlib import Path
from urllib.parse import unquote_to_bytes

import pandas as pd

from ._common import read_lines, read_tables, read_tsv, widen_attributes

FIELDS = (
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
)
PATTERN = r"\.(gff[23]?|gtf)$"


def read_gff(
    path: str | os.PathLike[str], fields: Sequence[str] | None = None
) -> pd.DataFrame:
    """Read GFF files into one table.

    Reads one GFF file, or every ``.gff``, ``.gff2``, ``.gff3`` and ``.gtf``
    file of a folder. Each attribute becomes a column. GFF3 attributes are
    read as ``key=value`` pairs, URL-decoded; GTF and GFF2 attributes, rows
    without any ``=``, as ``key "value"`` pairs, without the quotes. A key
    repeated in one row, such as ``tag`` in GTF files, keeps all its values
    joined with ``","``, as GFF3 writes several values of one key. Comment
    lines are left out, and reading stops at a ``##FASTA`` section.

    Parameters
    ----------
    path : str or path-like
        A GFF file, or a folder of them.
    fields : sequence of str, optional
        The columns to keep, e.g. ``["seqid", "type", "start", "end",
        "strand", "Name"]``. A folder's table keeps its ``filename`` column.

    Returns
    -------
    pandas.DataFrame
        One row per feature: ``seqid``, ``source``, ``type``, ``start``,
        ``end``, ``score``, ``strand``, ``phase``, then the attributes; for a
        folder, ``filename``, the file's name without its extension.

    Raises
    ------
    FileNotFoundError
        When the path does not exist, or a folder has no GFF files.
    ValueError
        When a field is not a column of the table.
    """
    data = read_tables(path, PATTERN, "gff or .gtf", process_gff)
    if fields is not None:
        fields = [fields] if isinstance(fields, str) else list(fields)
        invalid = [f for f in fields if f not in data.columns]
        if invalid:
            raise ValueError(f"Invalid field(s): {', '.join(invalid)}")
        columns = [*fields, "filename"] if Path(path).is_dir() else fields
        data = data[columns]
    return data


def process_gff(path: Path) -> pd.DataFrame:
    """Read one GFF file.

    Parameters
    ----------
    path : pathlib.Path
        The file.

    Returns
    -------
    pandas.DataFrame
        The table of `read_gff`, without ``filename``.
    """
    lines = read_lines(path)
    # Stop at the FASTA section (starts with ##FASTA or a line starting with >)
    for i, line in enumerate(lines):
        if re.match(r"^(##FASTA|>)", line):
            lines = lines[:i]
            break
    # Drop comment and directive lines; '#' inside a line is kept as data
    lines = [line for line in lines if not line.startswith("#")]
    data = read_tsv(lines, FIELDS)
    # The attributes as written, not as type_convert() read them
    rows = [line.split("\t") for line in lines if line != ""]
    attributes = [row[8] if len(row) > 8 else "" for row in rows]

    ids, keys, values = [], [], []
    for i, info in enumerate(attributes, start=1):
        for key, value in parse_attributes(info):
            ids.append(str(i))
            keys.append(key)
            values.append(value)
    data = data.drop(columns="attributes")
    if ids:
        record_ids, metadata = widen_attributes(ids, keys, values)
        position = {record: i for i, record in enumerate(record_ids)}
        index = [position.get(str(i)) for i in range(1, len(data) + 1)]
        rows_of = metadata.reindex([-1 if i is None else i for i in index]).reset_index(
            drop=True
        )
        data = pd.concat([data, rows_of], axis=1)
    return data


def parse_attributes(info: str | None) -> list[tuple[str, str | None]]:
    """Split the attributes of one row into key-value pairs.

    Parameters
    ----------
    info : str or None
        The ninth column of the row.

    Returns
    -------
    list of tuple
        ``(key, value)`` pairs; a key without a value has ``None``.
    """
    # No attributes: missing column, empty or "."
    if info is None or info == "NA" or info.strip() in ("", "."):
        return []
    pairs = [p.lstrip() for p in info.split(";")]
    pairs = [p for p in pairs if p != ""]

    # GTF and GFF2: key "value" pairs, in rows without any '=' (B48)
    if not any("=" in p for p in pairs) and any(re.search(r"\s", p) for p in pairs):
        result = []
        for pair in (p.strip() for p in pairs):
            if re.search(r"\s", pair):
                key = re.sub(r"\s.*$", "", pair, flags=re.S)
                value = re.sub(r"^\S+\s+", "", pair)
                value = re.sub(r'^"(.*)"$', r"\1", value, flags=re.S)
                result.append((key, value))
            else:
                result.append((pair, None))
        return result

    # GFF3: split each pair at the first '='; a key without '=' gets NA
    result = []
    for pair in pairs:
        if "=" in pair:
            key, value = pair.split("=", 1)
            result.append((key.strip(), url_decode(value)))
        else:
            result.append((pair.strip(), None))
    return result


def url_decode(value: str) -> str:
    """Decode ``%XX`` escapes as R's ``utils::URLdecode()``, as UTF-8.

    A ``%`` that is not followed by two hex digits is kept.

    Parameters
    ----------
    value : str
        The value.

    Returns
    -------
    str
        The decoded value.
    """
    if not re.search(r"%[0-9A-Fa-f]{2}", value):
        return value
    value = re.sub(r"%(?![0-9A-Fa-f]{2})", "%25", value)
    return unquote_to_bytes(value).decode("utf-8", errors="replace")
