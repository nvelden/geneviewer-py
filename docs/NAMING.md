# Names

How the R package's functions and arguments are named in Python, and the
names that read badly in Python, for review (decision D-d of the M7 plan:
such names change in the library and R first, then here, before 1.0).

## The rule

- Each `GC_*` function is a method of the chart, named after it without
  `GC_`, in snake_case: `GC_clusterTitle()` is `.cluster_title()`.
  `GC_chart()` is the constructor `GC()`. Every method returns a new chart.
- Arguments keep their R names, in snake_case: `itemStyle` is `item_style`,
  `marker_size` stays `marker_size`. The camelCase name works too
  (`itemStyle=`). Defaults are R's.
- Options passed through `**options` (R's `...`) go into the spec in
  camelCase: `font_size="10px"` and `fontSize="10px"` both give
  `fontSize`.
- Names inside dicts are not changed: `title_font={"fontSize": "16px"}`.
- In the spec, every name is as R writes it, so a Python chart and an R
  chart with the same arguments give the same spec (tested on 76 fixtures,
  `tests/parity`).

## R functions and their Python names

| R | Python | Notes |
|---|---|---|
| `GC_chart()` | `GC()` | No `elementId`: it sets the id of R's htmlwidget |
| `GC_align()` | `.align()` | Checks only; the library aligns |
| `GC_annotation()` | `.annotation()` | |
| `GC_chartOutput()` | `geneviewer.shiny.output_geneviewer()` | `outputId` is `id`, as in Shiny for Python's outputs |
| `renderGC_chart()` | `@geneviewer.shiny.render_geneviewer` | A decorator; the function's name is the output id |
| `GC_cluster()` | `.cluster()` | |
| `GC_clusterFooter()` | `.cluster_footer()` | |
| `GC_clusterLabel()` | `.cluster_label()` | |
| `GC_clusterTitle()` | `.cluster_title()` | |
| `GC_color()` | `.color()` | |
| `GC_coordinates()` | `.coordinates()` | |
| `GC_element_style()` | `.element_style()` | `cluster`, `part`, `end`, `reset`, `group` are keyword-only, as in R (after `...`) |
| `GC_events()` | `.events()` | |
| `GC_genes()` | `.genes()` | |
| `GC_grid()` | `.grid()` | |
| `GC_labels()` | `.labels()` | |
| `GC_legend()` | `.legend()` | |
| `GC_links()` | `.links()` | |
| `GC_normalize()` | `.normalize()` | |
| `GC_pattern()` | `.pattern()` | |
| `GC_scale()` | `.scale()` | |
| `GC_scaleBar()` | `.scale_bar()` | |
| `GC_sequence()` | `.sequence()` | |
| `GC_theme()` | `.theme()` | |
| `GC_title()` | `.title()` | |
| `GC_toolbar()` | `.toolbar()` | |
| `GC_tooltip()` | `.tooltip()` | |
| `GC_trackMouse()` | `.track_mouse()` | |
| `GC_transcript()` | `.transcript()` | |
| `GC_zoom()` | `.zoom()` | |
| `GC_zoom_to()` | `geneviewer.shiny.zoom_to()` | `outputId` is `id` |
| `geneviewer_theme()` | `geneviewer_theme()` | |
| `geneviewer_palette()` | `geneviewer_palette()` | |

Arguments whose spec name is not the camelCase of their Python name (the
method renames them, as R does):

| Method | Python argument | Spec name |
|---|---|---|
| `.legend()` | `legend_options`, `legend_text_options` | `markerStyle`, `textStyle` |
| `.sequence()` | `sequence_style` | `lineStyle` |
| `.scale()` | `scale_breaks`, `scale_break_threshold`, `scale_break_padding` | `autoBreaks`, `autoBreakThreshold`, `breakPadding` |
| `.scale_bar()` | `scale_bar_line_style`, `scale_bar_tick_style`; option `scale_bar_unit` | `lineStyle`, `tickStyle`, `unit` |
| `.links()` | `show_links`; option `label_adjustment_options` | `show`, `labelAdjustment` |
| `.labels()` | option `label_adjustment_options` | `labelAdjustment` |
| `.transcript()` | `style_utrs`, `item_style_utrs`, `group_style_utrs` | `styleUTRs`, `itemStyleUTRs`, `groupStyleUTRs` |
| `.element_style()` | `row_id` | `rowID` |
| all | `cluster` | `clusters`, always an array (not in `.element_style()`, nor where it passes through `**options`) |

## Other differences from R's API

- `GC()` takes a pandas DataFrame, or a GenBank file or folder as a path,
  as R's `GC_chart()` does.
- Readers return R's column types as pandas types: text as `string`,
  integers as `Int64`, numbers as `float64`, logicals as `boolean`; `NA`
  is `<NA>` (`NaN` in `float64`). `read_bed()`'s `blockSizes` and
  `blockStarts` are lists of floats.
- `read_gbk()` gives each feature as a dict of its qualifiers and `region`.
  R gives a named character vector, which can repeat a name (two
  `/db_xref`); in Python the later ones are `db_xref.1`, `db_xref.2`, the
  names `gbk_features_to_df()` gives their columns in both languages.
- `genbank_to_fasta()` returns the FASTA file as a `pathlib.Path`; R
  returns the path invisibly.
- `get_introns()` gives `start` and `end` as floats, as R's numbers.
- `protein_blast()` has no `parallel` argument (deprecated and ignored in
  R). It needs Biopython (`pip install "geneviewer-py[align]"`); its
  results can differ a little from R's, see `docs/differences.md`.
- `mummer_alignment()` returns `None` when no pair aligned, as R returns
  `NULL`. It runs MUMmer without a shell; the options strings are split as
  a shell would. Paths with spaces are refused, as in R (MUMmer's scripts
  need it).
- `.links(group=)` finds the links as R does, in R's row order (R's
  `merge()` reorders the genes of a group; Python follows it). Without
  `group` or `data`, it says "Either 'group' or 'data' must be provided.";
  R stops with an internal error there.
- R's `match.arg()` accepts the start of a choice (`"gen"` for `"genes"` in
  `GC_element_style()`); Python needs the whole name.
- Dates in data: a pandas datetime column whose values all fall at
  midnight is written as dates (`"2020-01-02"`, as R's `Date`); other
  datetimes as `"2020-01-02T03:04:05Z"` in UTC, as R's `POSIXct`.
- Messages name the Python arguments and values (`True`, `None`,
  `adjust_labels`); otherwise they are R's.
