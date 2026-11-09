"""Parity case: the chart of fixture cov-transcript-group-styles."""

from _inputs import before_transcript

from geneviewer import GC


def build(data):
    # BRCA1_splice_variants, rebuilt from the result until the datasets ship
    # (M7.8)
    variants = before_transcript(data, "ensembl_transcript_id")
    return GC(variants, height="500px").transcript(
        transcript="ensembl_transcript_id",
        strand="strand",
        type="type",
        style_by="ensembl_transcript_id",
        group_style_exons={"ENST00000357654": {"fill": "red"}},
        group_style_introns={"ENST00000357654": {"stroke": "red"}},
        group_style_utrs={"ENST00000357654": {"fill": "pink"}},
        opacity=0.8,
    )
