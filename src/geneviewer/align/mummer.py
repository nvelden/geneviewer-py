"""MUMmer alignment between clusters: ``mummer_alignment`` (``R/mummer_alignment.R``).

Runs MUMmer's ``nucmer`` or ``promer``, ``delta-filter`` and ``show-coords``
as R does, without a shell: paths are passed as they are, and the options
strings are split as a shell would split them.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ..readers._common import list_files, read_lines, read_tsv, warn
from ..readers.genbank import genbank_to_fasta

_INPUT = r"\.gbk$|\.gb$|\.fasta$"

NUCMER_COLUMNS = (
    "start1",
    "end1",
    "start2",
    "end2",
    "length1",
    "length2",
    "identity",
    "tag1",
    "tag2",
)
PROMER_COLUMNS = (
    "start1",
    "end1",
    "start2",
    "end2",
    "length1",
    "length2",
    "identity",
    "similarity",
    "stop",
    "frame1",
    "frame2",
    "tag1",
    "tag2",
)


def _run(command: list[str], stdout: Path | None, what: str) -> None:
    try:
        if stdout is None:
            done = subprocess.run(command, check=False)
        else:
            with open(stdout, "wb") as out:
                done = subprocess.run(command, stdout=out, check=False)
    except OSError:
        raise RuntimeError(f"Failed to run {what}.") from None
    if done.returncode != 0:
        raise RuntimeError(f"Failed to run {what}.")


def mummer(
    reference: str | os.PathLike[str],
    query: str | os.PathLike[str],
    maptype: str = "many-to-many",
    seqtype: str = "protein",
    mummer_options: str = "",
    filter_options: str = "",
    output_dir: str | os.PathLike[str] | None = None,
) -> Path:
    """Align two FASTA files with MUMmer and write the coordinates file.

    ``nucmer`` (nucleotide) or ``promer`` (protein), then ``delta-filter``
    (``-m`` many-to-many, ``-1`` one-to-one) and ``show-coords -H -T -r``
    (and ``-k`` for protein).

    Parameters
    ----------
    reference, query : path
        The FASTA files.
    maptype : str
        ``"many-to-many"`` or ``"one-to-one"``.
    seqtype : str
        ``"protein"`` or ``"nucleotide"``.
    mummer_options, filter_options : str
        More options for ``nucmer``/``promer`` and ``delta-filter``.
    output_dir : path, optional
        Where the files go. Default: the query's folder.

    Returns
    -------
    pathlib.Path
        The coordinates file,
        ``<reference>_<query>_<seqtype>_<maptype>.coords``.
    """
    if seqtype not in ("protein", "nucleotide"):
        raise ValueError(
            "Invalid sequence type. Please choose either 'protein' or 'nucleotide'."
        )
    command = "promer" if seqtype == "protein" else "nucmer"
    executable = shutil.which(command)
    if executable is None:
        raise RuntimeError(f"{command} is not installed or not in the PATH.")
    if maptype not in ("many-to-many", "one-to-one"):
        raise ValueError(
            "Invalid mapping type. Please choose either 'many-to-many' or 'one-to-one'."
        )
    reference, query = Path(reference), Path(query)
    if not query.exists():
        raise FileNotFoundError(f"Query file does not exist: {query}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference file does not exist: {reference}")

    folder = Path(output_dir) if output_dir is not None else query.parent
    prefix = folder / f"{reference.stem}_{query.stem}"
    maptype_name = maptype.replace("-", "_")
    _run(
        [
            executable,
            f"--prefix={prefix}_{seqtype}",
            *shlex.split(mummer_options),
            str(reference),
            str(query),
        ],
        None,
        command,
    )
    delta = Path(f"{prefix}_{seqtype}.delta")
    filtered = Path(f"{prefix}_{seqtype}_{maptype_name}.delta")
    coords = Path(f"{prefix}_{seqtype}_{maptype_name}.coords")
    flag = "-1" if maptype == "one-to-one" else "-m"
    _run(
        ["delta-filter", flag, *shlex.split(filter_options), str(delta)],
        filtered,
        "delta-filter",
    )
    options = ["-H", "-T", "-r"] + (["-k"] if seqtype == "protein" else [])
    _run(["show-coords", *options, str(filtered)], coords, "show-coords")
    return coords


def parse_coords(
    path: str | os.PathLike[str],
    reference: str,
    query: str,
    col_names: Sequence[str],
) -> pd.DataFrame:
    """Read a ``show-coords -T`` file (R's ``parse_coords``).

    Parameters
    ----------
    path : path
        The coordinates file.
    reference, query : str
        The cluster names, as ``cluster1`` and ``cluster2``.
    col_names : sequence of str
        The column names.

    Returns
    -------
    pandas.DataFrame
        The alignments; an empty table (without the cluster columns) for an
        empty file.
    """
    if not Path(path).exists():
        raise FileNotFoundError("The specified path does not exist.")
    lines = read_lines(path)
    if not lines:
        return pd.DataFrame({name: pd.Series(dtype="boolean") for name in col_names})
    data = read_tsv(lines, col_names)
    data["cluster1"] = pd.array([reference] * len(data), dtype="string")
    data["cluster2"] = pd.array([query] * len(data), dtype="string")
    return data


def parse_nucmer(
    path: str | os.PathLike[str], reference: str, query: str
) -> pd.DataFrame:
    """Read the coordinates of a ``nucmer`` alignment.

    Parameters
    ----------
    path : path
        The coordinates file.
    reference, query : str
        The cluster names.

    Returns
    -------
    pandas.DataFrame
        The alignments.
    """
    return parse_coords(path, reference, query, NUCMER_COLUMNS)


def parse_promer(
    path: str | os.PathLike[str], reference: str, query: str
) -> pd.DataFrame:
    """Read the coordinates of a ``promer`` alignment.

    Parameters
    ----------
    path : path
        The coordinates file.
    reference, query : str
        The cluster names.

    Returns
    -------
    pandas.DataFrame
        The alignments.
    """
    return parse_coords(path, reference, query, PROMER_COLUMNS)


def mummer_alignment(
    path: str | os.PathLike[str],
    cluster: Sequence[str] | None = None,
    maptype: str = "many-to-many",
    seqtype: str = "protein",
    mummer_options: str = "",
    filter_options: str = "",
    remove_files: bool = True,
    output_dir: str | os.PathLike[str] = tempfile.gettempdir(),  # noqa: B008
) -> pd.DataFrame | None:
    """Align neighbouring clusters with MUMmer, for links between them.

    Each cluster is a file in ``path`` named after it: ``<cluster>.fasta``,
    or a GenBank file (``.gbk``, ``.gb``) that is written to FASTA first.
    Clusters are aligned in order, each with the next. Needs MUMmer
    (``nucmer`` or ``promer``, ``delta-filter`` and ``show-coords``) in the
    PATH.

    The files are copied into a new folder in ``output_dir``, where MUMmer
    writes its files; with ``remove_files`` the folder is removed after.
    Files already in ``path`` or ``output_dir`` are never changed.

    Parameters
    ----------
    path : path
        The folder of sequence files. No spaces (MUMmer's scripts need it).
    cluster : list of str, optional
        The clusters, in order. Default: every file, by name.
    maptype : str
        ``"many-to-many"`` (default) or ``"one-to-one"``.
    seqtype : str
        ``"protein"`` (default, ``promer``) or ``"nucleotide"`` (``nucmer``).
    mummer_options : str
        More options for ``nucmer``/``promer``.
    filter_options : str
        More options for ``delta-filter``.
    remove_files : bool
        Whether to remove MUMmer's files afterwards. Default True.
    output_dir : path
        Where the working folder goes. Default: the temporary folder.

    Returns
    -------
    pandas.DataFrame or None
        One row per alignment: ``start1``, ``end1``, ``start2``, ``end2``,
        ``length1``, ``length2``, ``identity`` (and for protein
        ``similarity``, ``stop``, ``frame1``, ``frame2``), ``tag1``, ``tag2``,
        ``cluster1`` and ``cluster2``. ``None`` when no pair aligned. A pair
        that fails gives a warning and no rows.
    """
    path = Path(path)
    output_dir = Path(output_dir)
    if not path.is_dir():
        raise FileNotFoundError(
            "The specified directory does not exist. Please check the file path."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    if " " in str(path):
        raise ValueError("MUMmer requires a directory path without spaces.")
    if " " in str(output_dir):
        raise ValueError("MUMmer requires an output directory path without spaces.")

    all_files = list_files(path, _INPUT)
    if not all_files:
        raise FileNotFoundError(
            "No files found in the specified directory. Check the path or file types."
        )
    names = [re.sub(_INPUT, "", f.name) for f in all_files]
    if cluster is None:
        cluster = list(dict.fromkeys(names))
    cluster = [str(c) for c in ([cluster] if isinstance(cluster, str) else cluster)]
    if len(cluster) < 2:
        raise ValueError("At least two clusters must be provided.")
    spaced = [c for c in cluster if " " in c]
    if spaced:
        raise ValueError(
            "MUMmer requires file names without spaces: " + ", ".join(spaced)
        )

    # Every call has its own working folder, also when the input and output
    # folders are the same; files of the user are never overwritten or removed
    work_dir = Path(tempfile.mkdtemp(prefix="geneviewer-mummer-", dir=output_dir))
    try:
        files = []
        for file, name in zip(all_files, names, strict=True):
            if name in cluster:
                files.append(Path(shutil.copy(file, work_dir / file.name)))

        fasta = {
            f.name[: -len(".fasta")]: f for f in files if f.name.endswith(".fasta")
        }
        missing = [c for c in cluster if c not in fasta]
        if missing:
            for file in files:
                stem = re.sub(r"\.gbk$|\.gb$", "", file.name)
                if stem != file.name and stem in missing:
                    genbank_to_fasta(file, output_dir=work_dir)
            fasta = {
                c: work_dir / f"{c}.fasta"
                for c in cluster
                if (work_dir / f"{c}.fasta").exists()
            }
            missing = [c for c in cluster if c not in fasta]
            if missing:
                raise FileNotFoundError(
                    f"Sequence file(s) for: {', '.join(missing)} missing."
                )

        links = []
        for reference, query in zip(cluster[:-1], cluster[1:], strict=True):
            try:
                coords = mummer(
                    fasta[reference],
                    fasta[query],
                    maptype=maptype,
                    seqtype=seqtype,
                    mummer_options=mummer_options,
                    filter_options=filter_options,
                )
                parse = parse_nucmer if seqtype == "nucleotide" else parse_promer
                links.append(parse(coords, reference, query))
            except Exception as error:  # noqa: BLE001 - as R's tryCatch()
                warn(f"Error in aligning {reference} vs {query}: {error}")
    finally:
        if remove_files:
            shutil.rmtree(work_dir, ignore_errors=True)

    if not links:
        return None
    return pd.concat(links, ignore_index=True)
