# Differences from the R package

Every known difference in results between this package and the R package,
with its cause.

## `protein_blast()`: identity and similarity of some hits

**What differs.** For some protein pairs, `identity` and `similarity` differ
from R's, mostly by less than 1; among hits of 30 % identity or more, by at
most 4.4 in the data measured.
When two hits of a gene are close, the best hit can differ too, or a tie can
appear in one package and not in the other (a tie gives the gene a row per
hit). The order of the clusters was the same as R's in every run below.

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

Measured on the 418 proteins of 11 abyssomicin clusters (MIBiG GenBank
files, which the alignment was not tuned on), 4 queries, with the same R and
pwalign:

- `identity = 30` (the default), per query: `identity` differs in 0 to 5
  rows (at most 3.6), `similarity` in 15 to 33 rows (at most 4.4), and the
  best hit in 0 to 2 rows;
- `identity = 50`: `identity` differs in 0 to 2 rows, `similarity` in 2 to
  19, and the best hit in one row of one query;
- `identity = 20`: every query has 1 to 9 genes with another best hit or
  another number of tied hits;
- each of the 10920 alignments of the first query: the same alignment as
  pwalign's in 46 % (60 % for the erythromycin example).

On all five queries of `erythromycin_BlastP` at `identity = 30`: `identity`
differs in 0 to 2 of 148 rows (at most 1.9), `similarity` in 7 to 12 (at
most 2.8), and the best hit in none.

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

## `mummer_alignment()`: none known

With the same MUMmer, the tables are R's (checked on four phage genomes,
nucleotide and protein, both map types). MUMmer 3's `promer` itself can give
slightly different alignments from one run to the next on the same files, in
R and in Python alike.
