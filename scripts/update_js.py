"""Vendor a release of the geneviewer JavaScript library.

Usage, from anywhere::

    python scripts/update_js.py 0.3.10
    python scripts/update_js.py --dist <geneviewer-js>/dist 0.3.10.9000

Downloads the files the Python package ships from the GitHub release
``v<version>`` of nvelden/geneviewer-js, checks them against the release's
``SHA256SUMS``, and writes them to ``src/geneviewer/static/``, with the
version in ``static/VERSION`` and the release's ``SHA256SUMS``. Nothing in
the package changes until every file has been checked.

With ``--dist``, the files come from a build of geneviewer-js (``npm run
build``, which writes ``dist/SHA256SUMS``) instead of a release, for
development before a release. Give a development version, e.g. 0.3.10.9000;
vendor a release before the package is released.

``tests/unit/test_vendored.py`` checks the vendored files against the
vendored ``SHA256SUMS``. Written after the R package's ``dev/update_js.R``.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

FILES = (
    "geneviewer.min.js",
    "geneviewer.css",
    "themes.json",
    "spec-v1.schema.json",
    # Not used yet; kept for the later Jupyter (anywidget) work
    "geneviewer.standalone.mjs",
)
RELEASE = "https://github.com/nvelden/geneviewer-js/releases/download/v{version}"
VERSION_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+|-[0-9A-Za-z.]+)?$")
STATIC = Path(__file__).resolve().parent.parent / "src" / "geneviewer" / "static"


def parse_sums(text: str) -> dict[str, str]:
    """Read a ``SHA256SUMS`` file.

    Parameters
    ----------
    text : str
        The file's text: one ``<hash>  <file>`` line per file.

    Returns
    -------
    dict of str to str
        The expected hash of each file, by file name.
    """
    sums = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        match = re.match(r"^([0-9a-f]{64}) +\*?(.+)$", line)
        if not match:
            raise ValueError(f"Not a SHA256SUMS line: {line!r}")
        sums[match.group(2)] = match.group(1)
    return sums


def sha256(path: Path) -> str:
    """Return the SHA-256 hash of a file, in hex.

    Parameters
    ----------
    path : Path
        The file.

    Returns
    -------
    str
        The hash.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(version: str, dist: Path | None, into: Path) -> None:
    """Download or copy the files and ``SHA256SUMS`` into a folder.

    Parameters
    ----------
    version : str
        The release to download, without the ``v``.
    dist : Path or None
        A geneviewer-js ``dist`` folder to copy from instead.
    into : Path
        The folder to write to.
    """
    release = RELEASE.format(version=version)
    for file in (*FILES, "SHA256SUMS"):
        if dist is not None:
            print(f"Copying {dist / file}")
            if not (dist / file).is_file():
                sys.exit(f"Not found: {dist / file}")
            shutil.copyfile(dist / file, into / file)
            continue
        url = f"{release}/{file}"
        print(f"Downloading {url}")
        request = urllib.request.Request(url, headers={"User-Agent": "update_js.py"})
        try:
            with urllib.request.urlopen(request) as response:
                (into / file).write_bytes(response.read())
        except OSError as error:
            sys.exit(f"Download failed: {url} ({error})")


def check(folder: Path) -> None:
    """Check the files in a folder against its ``SHA256SUMS``.

    Parameters
    ----------
    folder : Path
        The folder holding the files and ``SHA256SUMS``.
    """
    expected = parse_sums((folder / "SHA256SUMS").read_text(encoding="utf-8"))
    missing = [file for file in FILES if file not in expected]
    if missing:
        sys.exit(f"SHA256SUMS does not list: {', '.join(missing)}")
    bad = [file for file in FILES if sha256(folder / file) != expected[file]]
    if bad:
        sys.exit(f"Checksum mismatch: {', '.join(bad)}")
    print(f"Checksums OK: {', '.join(FILES)}")


def main(argv: list[str] | None = None) -> None:
    """Vendor the library; see the module docstring.

    Parameters
    ----------
    argv : list of str, optional
        The command-line arguments, without the program name.
    """
    parser = argparse.ArgumentParser(description="Vendor a geneviewer-js release.")
    parser.add_argument("version", help="the release, e.g. 0.3.10")
    parser.add_argument("--dist", type=Path, help="copy from this dist folder")
    args = parser.parse_args(argv)
    if not VERSION_RE.match(args.version):
        parser.error(f"not a version: {args.version}")
    if args.dist is not None and not args.dist.is_dir():
        parser.error(f"no such folder: {args.dist}")

    with tempfile.TemporaryDirectory(prefix="geneviewer-js-") as tmp:
        download = Path(tmp)
        fetch(args.version, args.dist, download)
        check(download)
        STATIC.mkdir(parents=True, exist_ok=True)
        for file in (*FILES, "SHA256SUMS"):
            shutil.copyfile(download / file, STATIC / file)
    (STATIC / "VERSION").write_text(f"{args.version}\n", encoding="utf-8")

    source = f" from {args.dist}" if args.dist is not None else ""
    print(f"Vendored geneviewer-js {args.version}{source} in {STATIC}")


if __name__ == "__main__":
    main()
