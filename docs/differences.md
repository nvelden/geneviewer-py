# Differences from the R package

Every known difference in results between this package and the R package,
with its cause.

## `protein_blast()`: identity and similarity of some hits

**What differs.** For a small share of protein pairs, `identity` and
`similarity` differ from R's by a few tenths of a percent, rarely more. When
two hits of a gene are close, the best hit can differ too, or a tie can
appear in one package and not in the other (a tie gives the gene a row per
hit).

Measured on the 148 proteins of `erythromycin_BlastP` (5 clusters), against
R 0.3.0 (`m6-js-library`) with pwalign 1.0.0:

- query `BGC0000055`, `identity = 30` (the default): `identity`, `BlastP`,
  `score` and the row order are R's; `similarity` differs in 7 of 148 rows,
  by 0.1 to 0.8 (pinned in `tests/readers`, fixture `blast-erythromycin`);
- query `BGC0000054`, `identity = 50`: the same, with `similarity` different
  in 4 rows;
- query `BGC0000055`, six query genes, `identity = 20`: among weak hits
  (20 to 30 % identity), 3 genes get another best hit and one gene gets an
  extra tied hit (154 rows, R 153);
- each of the 2875 alignments behind the first case: the same score as
  pwalign's (pwalign computes in single precision, so up to 1e-5 lower);
  identity or similarity differ in 492 (17 %), by at most 4.0, more than 1.0
  in 25.

**Cause.** R aligns with `pwalign::pairwiseAlignment()`, Python with
Biopython's `PairwiseAligner`, with the same scores (Phred quality 22 for
every letter over the alphabet of all sequences, gaps of `10 + 4k`; see
`geneviewer.align.blast.aligner`). Where several alignments share the best
score, each library picks one by its own rules. pwalign fills its matrices
from the ends of the sequences, in single precision, preferring a
substitution to a deletion to an insertion. Python aligns the reversed
sequences, which makes the two agree on 60 % of the pairs instead of 0.3 %;
matching pwalign everywhere would need its algorithm ported cell by cell,
which is too slow in Python (1.5 billion cells for the example above).

Identity is counted on the alignment as pwalign returns it: without the
gaps at either end (Python does the same, so this is not a difference).

## `protein_blast()` with one cluster

R stops with a crash (a segmentation fault in pwalign) when the data has
only the query cluster, because there is nothing to align. Python returns
the data with each gene as its own hit; as in R, a hit only one gene has is
`"No Hit"`.
