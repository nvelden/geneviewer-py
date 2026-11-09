# geneviewer for Python

Gene cluster charts for Python and Shiny for Python, drawn by the
[geneviewer](https://github.com/nvelden/geneviewer-js) JavaScript library.
It is the Python counterpart of the
[geneviewer R package](https://github.com/nvelden/geneviewer).

Work in progress; not yet on PyPI.

```sh
pip install geneviewer-py
```

```python
import geneviewer
```

## Development

```sh
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
ruff check
ruff format --check
pytest
```

The parity tests compare chart specs with the R package's fixtures in a
geneviewer-js checkout. Set `GENEVIEWER_JS` to that checkout to run them;
without it they are skipped.

Browser tests need playwright's Chromium:

```sh
python -m playwright install chromium
pytest tests/browser -m browser
```

## License

MIT
