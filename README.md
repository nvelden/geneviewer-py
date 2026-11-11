# geneviewer for Python

[![PyPI](https://img.shields.io/pypi/v/geneviewer-py)](https://pypi.org/project/geneviewer-py/)
[![CI](https://github.com/nvelden/geneviewer-py/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/nvelden/geneviewer-py/actions/workflows/ci.yml)

Gene cluster charts for Python and Shiny for Python, drawn by the
[geneviewer](https://github.com/nvelden/geneviewer-js) JavaScript library.
It is the Python counterpart of the
[geneviewer R package](https://github.com/nvelden/geneviewer).

Documentation: <https://nvelden.github.io/geneviewer-docs/python/>

<p align="center"><a href="https://nvelden.github.io/geneviewer-shinylive/"><img src="https://raw.githubusercontent.com/nvelden/geneviewer/main/man/figures/geneviewer_promo.gif" width="100%" alt="geneviewer: load GenBank, FASTA or CSV files, add titles, labels and a legend, zoom in with the zoom bar, drag labels, switch the theme and export the SVG"/></a></p>

## Install

```sh
pip install "geneviewer-py[jupyter]"
```

The extras are `jupyter` (charts in notebooks), `shiny` (charts in Shiny for
Python) and `align` (`protein_blast()`, with Biopython). The import name is
`geneviewer`.

## A first chart

In a Jupyter notebook, a cell that ends with a chart draws it:

```python
from geneviewer import GC, datasets

genes = datasets.load("ophA_clusters")

(
    GC(genes, cluster="cluster", group="class", height="300px")
    .labels("name")
    .links(group="class")
    .legend(position="bottom")
)
```

`GC()` takes a table with a start and an end column, and each method adds
or changes a part of the chart, as the `GC_` functions of the R package do.
[docs/NAMING.md](https://github.com/nvelden/geneviewer-py/blob/main/docs/NAMING.md) lists the Python name of each R function
and argument.

To read what happens in the chart, keep it as a widget. Its attributes hold
the chart's events, and `zoom_to()` moves its zoom slider:

```python
widget = (
    GC(genes, cluster="cluster", group="class")
    .zoom()
    .events(select=True, zoom=True)
    .widget()
)
widget  # draws the chart

# In a later cell
widget.select  # the gene clicked last
widget.zoom_to([0, 0.5])
```

## Shiny for Python

With the `shiny` extra (`pip install "geneviewer-py[shiny]"`), charts are
shown in a Shiny for Python app. Save this as `app.py` and run
`shiny run app.py`:

```python
from shiny import App, ui

from geneviewer import GC, datasets
from geneviewer.shiny import output_geneviewer, render_geneviewer

genes = datasets.load("ophA_clusters")

app_ui = ui.page_fluid(output_geneviewer("chart", height="300px"))


def server(input, output, session):
    @render_geneviewer
    def chart():
        return (
            GC(genes, cluster="cluster", group="class", height="300px")
            .labels("name")
            .links(group="class")
            .legend(position="bottom")
        )


app = App(app_ui, server)
```

A chart's events are inputs named after the output, such as
`input.chart_select()` for a click on a gene, and `zoom_to()` moves the
zoom slider from the server. The example app of this repository shows both,
with three charts:

```sh
git clone https://github.com/nvelden/geneviewer-py
shiny run geneviewer-py/examples/shiny/app.py
```

## Example datasets

```python
from geneviewer import datasets

datasets.names()  # the 11 datasets of the R package
genes = datasets.load("ophA_clusters")  # a pandas DataFrame
```

## Reading files

`read_gbk()`, `read_gff()`, `read_fasta()` and `read_bed()` read a file, or
every file of a folder, into a table. `GC()` also takes the path of a
GenBank file or folder and draws its coding sequences:

```python
import geneviewer as gv

chart = gv.GC("BGC0000001.gbk").labels("protein_id")

records = gv.read_gbk("BGC0000001.gbk")
cds = gv.gbk_features_to_df(records, feature="CDS", keys=["region", "gene", "product"])
chart = gv.GC(cds, group="product")

genes = gv.read_gff("annotation.gff3")
proteins = gv.read_fasta("proteins.fasta")
```

## Alignment

`protein_blast()` compares the proteins of a query cluster with those of
the other clusters and gives each gene its best hit, to colour and link
genes by. It needs the `align` extra.

```python
import geneviewer as gv

hits = gv.protein_blast("folder/of/genbank/files", query="BGC0000055")
chart = gv.GC(hits, cluster="cluster", group="BlastP").links(
    group="BlastP", measure="identity"
)
```

`mummer_alignment()` aligns neighbouring clusters with MUMmer, which must
be installed, and gives the table `links(data=...)` draws:

```python
alignment = gv.mummer_alignment("folder/of/fasta/files", seqtype="nucleotide")
```

Results can differ a little from the R package's;
[docs/differences.md](https://github.com/nvelden/geneviewer-py/blob/main/docs/differences.md) lists where and why.

## Development

```sh
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,shiny]"
ruff check
ruff format --check
pytest
```

The parity tests compare chart specs with the R package's fixtures in a
geneviewer-js checkout. Set `GENEVIEWER_JS` to that checkout to run them;
without it they are skipped.

Browser tests run the example app, `examples/shiny/app.py`, in playwright's
Chromium:

```sh
python -m playwright install chromium
pytest tests/browser -m browser
```

## License

MIT
