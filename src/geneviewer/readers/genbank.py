"""GenBank reader, ported from R (``R/read_genbank.R``, ``genbank_to_fasta``)."""

from __future__ import annotations

import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from ._common import (
    STRING,
    bind_rows,
    dir_files,
    file_path_sans_ext,
    location_span,
    make_names,
    make_unique,
    read_lines,
    warn,
)

PATTERN = r"\.gbk$|\.gb$"

Feature = dict[str, str | None]
"""One feature: its qualifiers by name, and ``region``, its location.

A qualifier repeated in one feature (``/db_xref``) keeps each value, the
later ones as ``db_xref.1``, ``db_xref.2``, as R names them in
`gbk_features_to_df`.
"""


def read_gbk(
    path: str | os.PathLike[str],
    sections: Sequence[str] | None = None,
    features: Sequence[str] | None = None,
    origin: bool = True,
) -> dict[str, dict[str, Any]]:
    """Read GenBank files.

    Parameters
    ----------
    path : str or path-like
        A GenBank file (``.gbk`` or ``.gb``), or a folder of them.
    sections : sequence of str, optional
        The sections to read, e.g. ``["LOCUS", "FEATURES"]``; by default
        all of them.
    features : sequence of str, optional
        The feature keys to read from ``FEATURES``, e.g. ``["CDS"]``; by
        default all of them.
    origin : bool, default True
        Whether to read the sequence (``ORIGIN``).

    Returns
    -------
    dict
        One entry per file, by its name without the extension. Each is a
        dict of sections: text, except ``FEATURES``, a dict of feature key
        to a list of features (`Feature`), and ``ORIGIN``, the sequence.
        Only the first record of a file with several is read correctly; a
        warning says so.

    Raises
    ------
    FileNotFoundError
        When the path does not exist, or a folder has no GenBank files.
    """
    path = Path(path)
    if path.is_dir():
        files = dir_files(path, PATTERN, "gbk")
        return {
            re.sub(PATTERN, "", file.name): process_gbk(
                file, sections, features, origin, file.name
            )
            for file in files
        }
    if path.exists():
        name = re.sub(PATTERN, "", path.name)
        return {name: process_gbk(path, sections, features, origin, name)}
    raise FileNotFoundError("The specified path does not exist.")


def _setdiff(x: Sequence[str], y: Sequence[str]) -> list[str]:
    # R's setdiff(): unique values of x not in y
    return [v for v in dict.fromkeys(x) if v not in set(y)]


def process_gbk(
    file: Path,
    sections: Sequence[str] | None = None,
    features: Sequence[str] | None = None,
    origin: bool = True,
    basename: str | None = None,
) -> dict[str, Any]:
    """Read one GenBank file.

    Parameters
    ----------
    file : pathlib.Path
        The file.
    sections, features, origin
        As for `read_gbk`.
    basename : str, optional
        The name used in warnings.

    Returns
    -------
    dict
        The sections of the file.
    """
    lines = read_lines(file)
    if sum(line == "//" for line in lines) > 1:
        warn(
            f"Multiple sequence records found in {basename}. Only the first "
            "record will be processed."
        )
    actual_sections = gbk_get_section_keys(lines)
    section_keys = list(sections) if sections is not None else actual_sections
    if not origin:
        section_keys = _setdiff(section_keys, ["ORIGIN"])
    missing = _setdiff(section_keys, actual_sections)
    if missing:
        warn("The following sections were not found in the file: " + ", ".join(missing))
    values: dict[str, Any] = (
        gbk_get_sections(lines, section_keys) if section_keys else {}
    )

    feature_lines = values.get("FEATURES") or []
    actual_features = gbk_get_feature_keys(feature_lines)
    feature_keys = list(features) if features is not None else actual_features
    missing = _setdiff(feature_keys, actual_features)
    if missing:
        warn(
            "The following features were not found in the FEATURES section: "
            + ", ".join(missing)
        )
    values["FEATURES"] = {
        key: gbk_get_features(feature_lines, key) for key in feature_keys
    }
    return values


def gbk_get_feature_keys(lines: Sequence[str]) -> list[str]:
    """Return the feature keys of the ``FEATURES`` section, in order.

    Parameters
    ----------
    lines : sequence of str
        The lines of the section.

    Returns
    -------
    list of str
        The keys: words starting in column 6 with a letter.
    """
    keys = [
        re.sub(r"^ {5}(\S+).*", r"\1", line, flags=re.S)
        for line in lines
        if re.match(r"^ {5}[a-zA-Z]\w*", line)
    ]
    return list(dict.fromkeys(keys))


def gbk_get_section_keys(lines: Sequence[str]) -> list[str]:
    """Return the section keys of a GenBank file, in order.

    Parameters
    ----------
    lines : sequence of str
        The lines of the file.

    Returns
    -------
    list of str
        The keys: words starting in the first three columns.
    """
    found = []
    for line in lines:
        match = re.match(r"^ {0,2}([A-Za-z]+)(?=\s|$)", line)
        if match:
            found.append(match.group(0))
    return [key.strip() for key in dict.fromkeys(found)]


def gbk_get_sections(lines: Sequence[str], keys: Sequence[str]) -> dict[str, Any]:
    """Return the text of each section.

    Parameters
    ----------
    lines : sequence of str
        The lines of the file.
    keys : sequence of str
        The sections.

    Returns
    -------
    dict
        Each section found: its text on one line, without its key, except
        ``FEATURES`` (its lines) and ``ORIGIN`` (the letters of the sequence).
    """
    # Not the sequence lines of ORIGIN: from base 1,000,001 their position
    # has no more than two spaces before it
    section_start = re.compile(r"^ {0,2}(?![0-9])([^ ]+)(?=\s|$)")
    starts = [i for i, line in enumerate(lines) if section_start.match(line)]
    sections: dict[str, Any] = {}
    for key in keys:
        found = [i for i, line in enumerate(lines) if re.match(r"^\s{0,2}" + key, line)]
        if not found:
            continue
        later = [i for i in starts if i > found[0]]
        end = later[0] if later else len(lines)
        section = lines[found[0] : end]
        if key == "FEATURES":
            sections[key] = list(section)
            continue
        text = "\n".join(section)
        text = re.sub(r"^\s{0,2}" + key + r"\s+", "", text, count=1)
        sections[key] = re.sub(r"\s+", " ", text).strip()
    if sections.get("ORIGIN") is not None:
        sections["ORIGIN"] = re.sub(r"[^a-zA-Z]", "", sections["ORIGIN"])
    return sections


_QUALIFIER = re.compile(r'(?:^|(?<= ))/([^\s/="]+)=("(?:[^"]|"")*"|[^\s"]*)')


def gbk_process_features(text: str, key: str) -> Feature:
    """Return the qualifiers and the region of one feature.

    Parameters
    ----------
    text : str
        The feature's lines, on one line.
    key : str
        The feature key.

    Returns
    -------
    Feature
        The qualifiers, then ``region``.
    """
    names, values = [], []
    for match in _QUALIFIER.finditer(text):
        qualifier = match.group(0)
        name = re.sub(r"^/([^=]+)=.*$", r"\1", qualifier, flags=re.S)
        value = re.sub(r"^/[^=]+=", "", qualifier)
        value = re.sub(r'^"(.*)"$', r"\1", value, flags=re.S)
        names.append(name)
        values.append(value.replace('""', '"').strip())
    region_match = re.match(r".*?" + key + r"\s+(.*?)\s*/.*", text, flags=re.S)
    region = region_match.group(1) if region_match else None
    if "region" in names:
        values[names.index("region")] = region
    else:
        names.append("region")
        values.append(region)
    # Remove spaces from translation
    if "translation" in names:
        i = names.index("translation")
        values[i] = re.sub(r"[^a-zA-Z]", "", values[i])
    return dict(zip(make_unique(names), values, strict=True))


def gbk_get_features(lines: Sequence[str], key: str) -> list[Feature] | None:
    """Return every feature of one key.

    Parameters
    ----------
    lines : sequence of str
        The lines of the ``FEATURES`` section.
    key : str
        The feature key, e.g. ``"CDS"``.

    Returns
    -------
    list of Feature or None
        The features; ``None`` when there are none.
    """
    starts = [
        i for i, line in enumerate(lines) if re.match(r"^ {5}" + key + r"(\s|$)", line)
    ]
    if not starts:
        return None
    all_starts = [i for i, line in enumerate(lines) if re.match(r"^ {5}[a-zA-Z]", line)]
    result = []
    for start in starts:
        later = [i for i in all_starts if i > start]
        end = later[0] if later else len(lines)
        text = re.sub(r"\s+", " ", "\n".join(lines[start:end])).strip()
        result.append(gbk_process_features(text, key))
    return result


def gbk_features_to_df(
    gbk_list: Mapping[str, Mapping[str, Any]],
    feature: str = "CDS",
    keys: Sequence[str] | None = None,
    process_region: bool = True,
) -> pd.DataFrame:
    """Make a table of one kind of feature of GenBank files read with `read_gbk`.

    Parameters
    ----------
    gbk_list : dict
        From `read_gbk`.
    feature : str, default "CDS"
        The feature key.
    keys : sequence of str, optional
        The qualifiers to keep, e.g. ``["protein_id", "region"]``; by
        default all of them.
    process_region : bool, default True
        Whether to add ``strand`` (``"complement"`` or ``"forward"``),
        ``start`` and ``end`` from ``region``. ``start`` and ``end`` are the
        full span of the location, also of spliced genes listed from right
        to left.

    Returns
    -------
    pandas.DataFrame
        One row per feature, with ``cluster``, the file's name. Qualifier
        names are made syntactic as in R (``function`` becomes
        ``function.``, ``db-xref`` becomes ``db.xref``). Rows and columns
        that are all missing are left out.
    """
    tables = []
    for cluster_name, item in gbk_list.items():
        found = item.get("FEATURES") or {}
        if feature not in found:
            warn(f"Feature {feature} not found in cluster {cluster_name}")
            continue
        feature_list = found[feature] or []
        if keys is not None and not any(
            any(k in f for k in keys) for f in feature_list
        ):
            warn(f"None of the specified keys found in {cluster_name}")
        rows = []
        for f in feature_list:
            if keys is not None:
                names: list[str | None] = [k if k in f else None for k in keys]
                values = [f.get(k) for k in keys]
            else:
                names, values = list(f), list(f.values())
            rows.append(
                pd.DataFrame(
                    {
                        n: pd.array([v], dtype=STRING)
                        for n, v in zip(make_names(names), values, strict=True)
                    }
                )
            )
        df = bind_rows(rows)
        all_na = (
            df.isna().all(axis=1) if len(df.columns) else pd.Series([True] * len(df))
        )
        df["cluster"] = pd.array([cluster_name] * len(df), dtype=STRING)
        df = df[~all_na.to_numpy()].reset_index(drop=True)
        df = df[[c for c in df.columns if df[c].isna().sum() < len(df)]]
        if process_region and "region" in df.columns:
            regions = [None if pd.isna(r) else r for r in df["region"].tolist()]
            df["strand"] = pd.array(
                [
                    None
                    if r is None
                    else ("complement" if "complement" in r.lower() else "forward")
                    for r in regions
                ],
                dtype=STRING,
            )
            spans = [location_span(r) for r in regions]
            df["start"] = [s[0] for s in spans]
            df["end"] = [s[1] for s in spans]
        tables.append(df)
    return bind_rows(tables)


def genbank_to_fasta(
    path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = tempfile.gettempdir(),  # noqa: B008
) -> Path:
    """Write the sequences of a GenBank file to a FASTA file.

    Each record becomes a FASTA record with the header ``>VERSION
    DEFINITION`` and the sequence of ``ORIGIN`` in capitals.

    Parameters
    ----------
    path : str or path-like
        The GenBank file.
    output_dir : str or path-like or None, optional
        The folder of the FASTA file, which gets the input's name with the
        extension ``.fasta``. By default the temporary folder; ``None``
        writes it next to the input.

    Returns
    -------
    pathlib.Path
        The FASTA file.

    Raises
    ------
    FileNotFoundError
        When the file does not exist.
    ValueError
        When the file has no sequence, or its records do not end in ``//``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError("File does not exist")
    lines = read_lines(path)
    ends = [i for i, line in enumerate(lines) if line.startswith("//")]
    if not any(line.startswith("ORIGIN") for line in lines):
        raise ValueError("No sequence detected in GenBank file.")
    if not ends:
        raise ValueError("The sequence block is not properly formatted or missing.")
    starts = [0, *(i + 1 for i in ends[:-1])]
    records = []
    for r, (start, end) in enumerate(zip(starts, ends, strict=True), start=1):
        record = lines[start : end + 1]
        definition = [
            re.sub(r"[^a-zA-Z0-9]+$", "", re.sub(r"^DEFINITION\s+", "", line))
            for line in record
            if line.startswith("DEFINITION")
        ]
        version = [
            re.sub(r"^VERSION\s+", "", line)
            for line in record
            if line.startswith("VERSION")
        ]
        origin = [i for i, line in enumerate(record) if line.startswith("ORIGIN")]
        if not origin:
            warn(f"Record {r} has no sequence and is skipped.")
            continue
        if not definition or not version:
            continue
        sequence = "".join(record[origin[0] + 1 : len(record) - 1])
        sequence = re.sub(r"\s+|[0-9]+", "", sequence).upper()
        records.append(f">{version[0]} {definition[0]}\n{sequence}")
    folder = path.parent if output_dir is None else Path(output_dir)
    fasta = folder / f"{file_path_sans_ext(path)}.fasta"
    fasta.write_text("\n".join(records) + "\n", encoding="utf-8")
    return fasta


def read_gbk_cds(path: str | os.PathLike[str], start: str, end: str) -> pd.DataFrame:
    """Read the CDS of GenBank files for a chart, as R's ``GC_chart()`` does.

    Parameters
    ----------
    path : str or path-like
        A GenBank file, or a folder of them.
    start, end : str
        The columns of the start and end; rows without them are left out.

    Returns
    -------
    pandas.DataFrame
        ``protein_id``, ``region``, ``translation``, ``cluster``, ``strand``,
        ``start`` and ``end`` of each CDS.
    """
    data = gbk_features_to_df(
        read_gbk(path), feature="CDS", keys=["protein_id", "region", "translation"]
    )
    keep = data[start].notna() & data[end].notna()
    return data[keep.to_numpy()].reset_index(drop=True)
