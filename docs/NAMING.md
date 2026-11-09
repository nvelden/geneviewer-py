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
| `GC_chartOutput()` | `geneviewer.shiny.output_geneviewer()` | M7.4 |
| `renderGC_chart()` | `@geneviewer.shiny.render_geneviewer` | M7.4 |
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
| `GC_links()` | `.links()` | `group=` without `data=` in M7.6 |
| `GC_normalize()` | `.normalize()` | M7.6 |
| `GC_pattern()` | `.pattern()` | |
| `GC_scale()` | `.scale()` | |
| `GC_scaleBar()` | `.scale_bar()` | |
| `GC_sequence()` | `.sequence()` | |
| `GC_theme()` | `.theme()` | |
| `GC_title()` | `.title()` | |
| `GC_toolbar()` | `.toolbar()` | |
| `GC_tooltip()` | `.tooltip()` | |
| `GC_trackMouse()` | `.track_mouse()` | |
| `GC_transcript()` | `.transcript()` | M7.5/M7.6 (introns) |
| `GC_zoom()` | `.zoom()` | |
| `GC_zoom_to()` | `geneviewer.shiny.zoom_to()` | M7.4 |
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

- `GC()` takes a pandas DataFrame; R's `GC_chart()` also reads a folder of
  GenBank files given as a path. Python gets that with its GenBank reader
  (M7.5).
- R's `match.arg()` accepts the start of a choice (`"gen"` for `"genes"` in
  `GC_element_style()`); Python needs the whole name.
- Dates in data: a pandas datetime column whose values all fall at
  midnight is written as dates (`"2020-01-02"`, as R's `Date`); other
  datetimes as `"2020-01-02T03:04:05Z"` in UTC, as R's `POSIXct`.
- Messages name the Python arguments and values (`True`, `None`,
  `adjust_labels`); otherwise they are R's.
