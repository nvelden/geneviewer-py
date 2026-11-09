import hashlib
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "src" / "geneviewer" / "static"
FILES = (
    "geneviewer.min.js",
    "geneviewer.css",
    "themes.json",
    "spec-v1.schema.json",
    "geneviewer.standalone.mjs",
)


def read_sums() -> dict[str, str]:
    sums = {}
    for line in (STATIC / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = re.match(r"^([0-9a-f]{64}) +\*?(.+)$", line).groups()
            sums[name] = digest
    return sums


def test_version_file():
    version = (STATIC / "VERSION").read_text(encoding="utf-8").strip()
    assert re.match(r"^[0-9]+\.[0-9]+\.[0-9]+", version)


def test_files_match_sha256sums():
    sums = read_sums()
    for file in FILES:
        digest = hashlib.sha256((STATIC / file).read_bytes()).hexdigest()
        assert digest == sums[file], f"{file} does not match SHA256SUMS"


def test_wheel_includes_static_files(tmp_path):
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            "--no-deps",
            "--no-build-isolation",
            "--quiet",
            "--wheel-dir",
            str(tmp_path),
            str(ROOT),
        ],
        check=True,
    )
    (wheel,) = tmp_path.glob("*.whl")
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
    for file in (*FILES, "SHA256SUMS", "VERSION", "geneviewer-shiny.js"):
        assert f"geneviewer/static/{file}" in names, f"{file} is not in the wheel"
