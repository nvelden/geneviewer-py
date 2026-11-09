"""Alignment helpers, with the values of R and pwalign."""

import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

from geneviewer import mummer_alignment, protein_blast
from geneviewer.align.blast import (
    align_pairs,
    compute_identity,
    substitution_scores,
    synteny_score,
)

pytest.importorskip("Bio")

MUMMER_INPUTS = Path(__file__).resolve().parents[1] / "readers" / "inputs" / "mummer"


@pytest.mark.parametrize(
    ("n", "match", "mismatch"),
    [
        # pwalign: qualitySubstitutionMatrices(c(0, 1), n)["22", "22", ]
        (2, 0.981794907143049, -5.31737343780391),
        (4, 1.98175613057263, -5.89928563552943),
        (20, 4.30366789817117, -6.2390401421848),
        (23, 4.50530134187931, -6.24887761222764),
    ],
)
def test_substitution_scores_as_pwalign(n, match, mismatch):
    assert substitution_scores(n) == pytest.approx((match, mismatch), rel=1e-12)


def test_alignments_leave_out_end_gaps():
    # As pwalign's pattern() and subject() of these alignments
    assert align_pairs(["WWWMKVLAAA", "MKVLAAA"], ["MKVLAAA", "WMKVLAAAW"]) == [
        ("MKVLAAA", "MKVLAAA"),
        ("MKVLAAA", "MKVLAAA"),
    ]
    assert align_pairs([], []) == []
    with pytest.raises(ValueError, match="same length"):
        align_pairs(["A"], [])


def test_compute_identity():
    result = compute_identity(["MKV-LP", "AAAA", "M--K"], ["MRV-IG", "AAAA", "MK-K"])
    # Both gaps do not count: 2 of 5 identical, 4 of 5 similar (K/R, L/I)
    assert result["identity"].tolist() == [40.0, 100.0, pytest.approx(66.7)]
    assert result["similarity"].tolist() == [80.0, 100.0, pytest.approx(66.7)]
    with pytest.raises(ValueError, match="same length"):
        compute_identity(["AA"], ["A"])


def test_synteny_score():
    # Pairs (1,2) and (3,4) are neighbours in both; 5 is not
    assert synteny_score([1, 2, 3, 4, 5], [1, 2, 3, 4, 9], [50] * 5) == 2 + 1.25
    assert synteny_score([1], [7], [100]) == 0.5


def test_protein_blast_one_cluster():
    # R crashes here (pwalign gets nothing to align); each gene is its own
    # hit, and a hit only one gene has is "No Hit"
    data = pd.DataFrame(
        {
            "protein_id": ["a", "b"],
            "cluster": "A",
            "start": [1, 10],
            "end": [5, 15],
            "translation": ["MKV", "MKL"],
        }
    )
    result = protein_blast(data, query="A")
    assert result["identity"].tolist() == [100.0, 100.0]
    assert result["BlastP"].tolist() == ["No Hit", "No Hit"]


def test_protein_blast_needs_biopython(monkeypatch):
    monkeypatch.setitem(sys.modules, "Bio.Align", None)
    with pytest.raises(ImportError, match=r"geneviewer-py\[align\]"):
        protein_blast(pd.DataFrame(), query="A")


def test_mummer_failing_pair_warns(monkeypatch, tmp_path):
    # Each pair that fails gives a warning; the input files stay as they were
    folder = tmp_path / "in"
    shutil.copytree(MUMMER_INPUTS, folder)
    before = sorted(p.name for p in folder.iterdir())
    with pytest.warns(UserWarning, match="Error in aligning B vs C: Invalid sequence"):
        result = mummer_alignment(
            folder, cluster=["B", "C"], seqtype="dna", output_dir=folder
        )
    assert result is None
    assert sorted(p.name for p in folder.iterdir()) == before


def test_mummer_paths_without_spaces(monkeypatch, tmp_path):
    folder = tmp_path / "with space"
    folder.mkdir()
    with pytest.raises(ValueError, match="directory path without spaces"):
        mummer_alignment(folder)
    # A relative path: this checkout's path can have spaces
    monkeypatch.chdir(MUMMER_INPUTS.parent)
    with pytest.raises(ValueError, match="file names without spaces"):
        mummer_alignment("mummer", cluster=["A", "B 2"], output_dir=tmp_path)
