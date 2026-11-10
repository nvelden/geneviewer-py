"""FASTA reader, ported from R's ``read_fasta()`` (``R/read_fasta.R``)."""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ._common import (
    STRING,
    bind_rows,
    file_path_sans_ext,
    list_files,
    location_span,
    read_lines,
    warn,
    widen_attributes,
)

# Key-value pairs: "[key=value]" (the value can contain spaces and brackets;
# it ends at the "]" before the next "[key=" or the end), or without brackets
# "key=value" running up to the next " key=" or the end
_BRACKETED = re.compile(r"\[\w+=.*?\](?=\s*\[\w+=|\s*$)", re.ASCII)
_UNBRACKETED = re.compile(r"(\w+)=(.*?)(?=\s+\w+=|$)", re.ASCII)


def read_fasta(
    fasta_path: str | os.PathLike[str],
    sequence: bool = True,
    keys: Sequence[str] | None = None,
    file_extension: str = "fasta",
) -> pd.DataFrame:
    """Read FASTA files into one table, a row per sequence.

    The keys of each header become columns: ``[key=value]`` pairs, as NCBI
    writes them, or ``key=value`` pairs, as in UniProt headers. A key
    repeated in one header keeps its values joined with ``","``. With a
    ``location`` key, the table gets ``start`` and ``end``, the full span of
    the location, and ``strand`` when a location is on the complement
    strand. Sequences without any key are left out, with a warning; a file
    that cannot be read is skipped, with a warning.

    Parameters
    ----------
    fasta_path : str or path-like
        A FASTA file (also gzip-compressed), or a folder of them.
    sequence : bool, default True
        Whether to keep the sequences, in the column ``sequence``.
    keys : sequence of str, optional
        The columns to keep.
    file_extension : str, default "fasta"
        The extension of the files read from a folder.

    Returns
    -------
    pandas.DataFrame
        One row per sequence, with ``cluster``, the file's name without its
        extension.

    Raises
    ------
    FileNotFoundError
        When the path does not exist.
    """
    path = Path(fasta_path)
    if path.is_dir():
        files = list_files(path, f"\\.{file_extension}$")
    elif path.exists():
        files = [path]
    else:
        raise FileNotFoundError(
            "Provided path does not exist or is not a valid FASTA file or directory."
        )
    tables = []
    for file in files:
        try:
            table = _read_fasta_file(file, sequence, keys)
        except Exception as error:  # noqa: BLE001 (R's tryCatch)
            warn(f"An error occurred while processing file '{file}': {error}")
            continue
        if table is not None:
            tables.append(table)
    return bind_rows(tables)


def _read_fasta_file(
    file: Path, sequence: bool, keys: Sequence[str] | None
) -> pd.DataFrame | None:
    cluster_name = file_path_sans_ext(file)
    headers, sequences = read_fasta_records(file)

    bracketed = [bool(_BRACKETED.search(h)) for h in headers]
    matches = [
        list((_BRACKETED if b else _UNBRACKETED).finditer(h))
        for h, b in zip(headers, bracketed, strict=True)
    ]
    no_matches = [i for i, m in enumerate(matches) if not m]
    if len(sequences) == len(no_matches):
        warn(
            f"No keys found in '{cluster_name}'. Make sure the fasta headers are "
            "in the right format. Skipping file."
        )
        return None
    if no_matches:
        warn(
            f"No keys found in '{cluster_name}' for seq(s) at index(es): "
            f"{', '.join(str(i + 1) for i in no_matches)}. They will be removed."
        )
        keep = [i for i in range(len(headers)) if i not in set(no_matches)]
        sequences = [sequences[i] for i in keep]
        matches = [matches[i] for i in keep]
        bracketed = [bracketed[i] for i in keep]

    ids, names, values = [], [], []
    for i, (found, in_brackets) in enumerate(zip(matches, bracketed, strict=True)):
        for match in found:
            pair = match.group(0)
            if in_brackets:
                pair = pair[1:-1]
            key, _, value = pair.partition("=")
            ids.append(str(i + 1))
            names.append(key)
            values.append(value.strip())
    _, df = widen_attributes(ids, names, values)
    # 'NA' as a value is NA
    df = df.mask((df == "NA").fillna(False))

    if sequence:
        df["sequence"] = pd.array(sequences, dtype=STRING)

    if "location" in df.columns:
        spans = [
            location_span(None if pd.isna(v) else v) for v in df["location"].tolist()
        ]
        df["start"] = [s[0] for s in spans]
        df["end"] = [s[1] for s in spans]
        complement = [
            not pd.isna(v) and "complement" in v for v in df["location"].tolist()
        ]
        if any(complement):
            df["strand"] = pd.array(
                ["complement" if c else "forward" for c in complement], dtype=STRING
            )

    df["cluster"] = pd.array([cluster_name] * len(df), dtype=STRING)

    if keys is not None:
        keys = [keys] if isinstance(keys, str) else list(keys)
        missing = [k for k in keys if k not in df.columns]
        if missing:
            warn(
                "The following keys are not column names in the dataframe: "
                f"{', '.join(missing)}."
            )
        df = df[[k for k in dict.fromkeys(keys) if k in df.columns]]
    return df


def read_fasta_records(path: str | os.PathLike[str]) -> tuple[list[str], list[str]]:
    """Read the headers and sequences of a FASTA file.

    Comment lines (``;``) and empty lines are left out; the lines of a
    sequence are joined, upper-cased, and characters that are not one-letter
    codes are removed (letters, ``*``, ``+``, ``.`` and ``-`` are kept), with
    a warning that counts them.

    Parameters
    ----------
    path : str or path-like
        The file.

    Returns
    -------
    tuple of list of str
        The headers, without ``>``, and the sequences.

    Raises
    ------
    ValueError
        When the first line is not a header.
    """
    numbered = [
        (n, line)
        for n, line in enumerate(read_lines(path), start=1)
        if not line.startswith(";") and line != ""
    ]
    if numbered and not numbered[0][1].startswith(">"):
        raise ValueError(
            f'reading FASTA file {path}: ">" expected at beginning of line '
            f"{numbered[0][0]}"
        )
    headers: list[str] = []
    parts: list[list[str]] = []
    for _, line in numbered:
        if line.startswith(">"):
            headers.append(line[1:])
            parts.append([])
        else:
            parts[-1].append(line)
    raw = ["".join(p) for p in parts]
    sequences = [re.sub(r"[^A-Za-z*+.\-]", "", s, flags=re.ASCII).upper() for s in raw]
    invalid = sum(
        len(r.encode("utf-8")) - len(s.encode("utf-8"))
        for r, s in zip(raw, sequences, strict=True)
    )
    if invalid > 0:
        warn(
            f"reading FASTA file {path}: ignored {invalid} invalid one-letter "
            "sequence codes"
        )
    return headers, sequences
