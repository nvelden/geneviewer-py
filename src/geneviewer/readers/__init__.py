"""File readers: GenBank, GFF, FASTA and BED, ported from the R package.

They return the tables R's readers return: the same columns, in the same
order, with the same types (text as ``string``, R's integers as ``Int64``,
R's numbers as ``float64``, logicals as ``boolean``).
"""

from .bed import read_bed
from .fasta import read_fasta
from .genbank import gbk_features_to_df, genbank_to_fasta, read_gbk
from .gff import read_gff

__all__ = [
    "gbk_features_to_df",
    "genbank_to_fasta",
    "read_bed",
    "read_fasta",
    "read_gbk",
    "read_gff",
]
