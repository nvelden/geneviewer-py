# geneviewer-py (development version)

- `normalize()` spaces the genes over the whole span of their cluster, as
  documented. It shared the total gap between as many gaps as there are
  genes, one too many, so every gap was too small and the cluster ended
  before its last gene did: genes at 10-40, 90-120 and 130-160 became 10-40,
  60-90 and 110-140, and are now 10-40, 70-100 and 130-160. Genes that are
  already evenly spaced stay where they are, and a given `gap` is the gap
  drawn; it was (n - 1) / n of it for n genes. Charts made with
  `normalize()` change. The R package has the same fix.
- `links(group=)` works on data with a column named `cluster` when the
  clusters of the chart are another column, e.g. `GC(data,
  cluster="genome")`; it stopped with "cannot reindex on an axis with
  duplicate labels". That column is left out of the links. Links by that
  column itself (`group="cluster"`) stop with a message to rename it.
- `read_fasta()` gives every sequence with a `location` a `strand`, also
  in files whose locations are all forward. Before, such a file had no
  `strand` column, and in a folder with other files its rows got `None`,
  which `GC(strand="strand")` rejects.
- `read_gff()` and `read_bed()` give the column types of the R package
  0.2.1, so the files of a folder can be combined when a column holds
  numbers in one file and text in another (e.g. chromosomes `1` and `X`).
  `seqid`, `source`, `type`, `strand` and `phase` of GFF files, and
  `chrom`, `name`, `strand` and `itemRgb` of BED files are always text;
  `score` is always numbers, `NaN` for a score of `"."`. A score that is not
  a number becomes `NaN`, with a warning.
- Links made with `links(group=)` before `normalize()` move with their
  genes; they stayed at the old positions. Links given as `data` (e.g.
  alignments) keep their positions, with a warning to add them after
  `normalize()`.
- `zoom(range=)` and `zoom(min=)` are checked in the `mode` of an earlier
  `zoom()` when the call gives none: `.zoom(mode="absolute").zoom(range=[100,
  2000])` stopped with the error for a relative range.
- `show=False` in `labels()`, `coordinates()`, `scale_bar()`,
  `cluster_title()`, `cluster_label()` and `cluster_footer()` hides what an
  earlier call shows, e.g. `.coordinates().coordinates(show=False,
  cluster=2)`; it changed nothing. Hidden labels need no label column. With
  nothing shown before, the chart is unchanged, as before.
- A chart of data without rows (or only rows without a start or end) takes
  the chart's methods: `labels("name")` and `legend()` stopped with "column
  not found in data", and `normalize()` and `transcript()` with an internal
  error. The chart has no clusters. (Drawing it needs the next geneviewer
  JavaScript library.)
- `genbank_to_fasta()` writes records without a `VERSION` line, named by
  their `ACCESSION`, else their `LOCUS` name, else the file name; without a
  `DEFINITION`, the header is the name only. They were left out without a
  message, so the FASTA file could be empty.
- `read_gbk()` reads features whose key does not start with a letter, such
  as `5'UTR`, `3'UTR` and `-10_signal`; their lines were added to the
  feature before them. Qualifiers without a value, such as `/pseudo`, are
  kept, with the value `""`, and a feature without qualifiers gets its
  `region`.
- `read_gbk()` reads the whole sequence (`ORIGIN`) of records of more than
  a million bases; it stopped at base 1,000,020, without a message.

# geneviewer-py 0.2.0

- Charts in Jupyter notebooks, with the new `jupyter` extra
  (`pip install "geneviewer-py[jupyter]"`): a cell that ends with a chart
  draws it. It works where anywidget does: JupyterLab, Jupyter Notebook, VS
  Code, Google Colab and marimo.
- `chart.widget()` gives the chart as a widget (`geneviewer.jupyter.GCWidget`)
  to keep in a variable. Its attributes hold the chart's events, as the
  inputs of a Shiny app do (`widget.select`, `widget.zoom`, ...),
  `widget.on("select", callback)` calls a function at each event,
  `widget.zoom_to()` moves the zoom slider, and setting `widget.chart` draws
  another chart in the same place.

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
