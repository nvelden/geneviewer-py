# geneviewer-py 0.1.0

The first release: the charts of the geneviewer R package and JavaScript
library from Python, in Shiny for Python.

- `GC()` makes a chart from a pandas DataFrame of genes, or from a GenBank
  file or folder. Its methods are the `GC_` functions of the R package
  (`.labels()`, `.links()`, `.transcript()`, `.theme()`, ...), with the same
  arguments in snake_case. A chart in Python and one in R with the same
  arguments give the same chart.
- `geneviewer.shiny`: `output_geneviewer()` and `@render_geneviewer` show a
  chart in a Shiny for Python app, its events are inputs such as
  `input.chart_select()`, and `zoom_to()` moves the zoom slider from the
  server.
- File readers: `read_gbk()`, `gbk_features_to_df()`, `read_gff()`,
  `read_fasta()`, `read_bed()` and `genbank_to_fasta()`, with the results of
  the R package's readers.
- Alignment: `protein_blast()` (with the `align` extra) and
  `mummer_alignment()` (with MUMmer installed).
- `geneviewer.datasets`: the 11 example datasets of the R package.
- `geneviewer_theme()` and `geneviewer_palette()` give the settings of the
  built-in themes and the colours of the palettes.

The package draws with version 0.3.10 of the geneviewer JavaScript library,
which it includes.
