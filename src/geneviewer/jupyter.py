"""Charts in Jupyter notebooks.

With the ``jupyter`` extra (``pip install "geneviewer-py[jupyter]"``), a
chart that a cell ends with is drawn below the cell::

    from geneviewer import GC, datasets

    GC(datasets.load("ophA_clusters"), cluster="cluster", group="class")

`GCWidget` is the widget that draws it; ``chart.widget()`` makes one. Keep
it in a variable to read what happens in the chart and to change it::

    widget = GC(genes, group="class").zoom().events(select=True, zoom=True).widget()
    widget          # draws the chart
    widget.select   # in a later cell: the gene clicked last

**Events.** The events the chart's ``.events()`` switches on set
attributes of the widget, the inputs of a Shiny app by another name:

- ``click``: the element clicked (``.events(click=True)``, the default).
- ``dblclick``: the element double-clicked (``dblclick=True``).
- ``select``: the element selected last (``select=True``).
- ``selection``: every selected element, as a list (``select=True``).
- ``drag``: an element moved by dragging, and its offset (``drag=True``).
- ``moves``: the elements a reset, undo or redo moved (``drag=True``).
- ``zoom``: the window of the zoom slider (``zoom=True``).
- ``rendered``: the chart was drawn (``rendered=True``).

Each is ``None`` until its event happens. ``widget.on(name, callback)``
calls a function at each event, also when the same event happens twice.
The payload of each event is described at
https://nvelden.github.io/geneviewer-docs/api#events.

The widget works where anywidget does: JupyterLab, Jupyter Notebook, VS
Code, Google Colab and marimo.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from functools import cache
from importlib import resources
from typing import Any

try:
    import anywidget
    import traitlets
except ImportError as error:  # pragma: no cover
    raise ImportError(
        "geneviewer.jupyter needs anywidget: pip install 'geneviewer-py[jupyter]'"
    ) from error

from ._checks import zoom_window
from .chart import GC

__all__ = ["EVENTS", "GCWidget"]

#: The events the widget reports, as its attributes of the same name: those
#: the Shiny binding sends as ``input.<id>_<event>()``.
EVENTS = (
    "click",
    "dblclick",
    "select",
    "selection",
    "drag",
    "moves",
    "zoom",
    "rendered",
)


@cache
def _module() -> str:
    """Return the widget's JavaScript: the library, then the widget's code.

    anywidget loads one module, without imports of files beside it, so the
    library (``geneviewer.standalone.mjs``, with D3) and the widget's code
    (``geneviewer-jupyter.js``) are joined. The library's build shortens its
    names: the widget's code gets ``newChart`` under the name the library's
    export statement gives it.
    """
    static = resources.files("geneviewer") / "static"
    library = (static / "geneviewer.standalone.mjs").read_text(encoding="utf-8")
    widget = (static / "geneviewer-jupyter.js").read_text(encoding="utf-8")
    exports = re.findall(r"export\s*\{([^}]*)\}", library)
    names = dict(
        reversed(re.split(r"\s+as\s+", part.strip()))
        if " as " in part
        else (part.strip(), part.strip())
        for part in ",".join(exports).split(",")
        if part.strip()
    )
    if "newChart" not in names:
        raise RuntimeError("geneviewer.standalone.mjs does not export newChart.")
    library = re.sub(r"\n//# sourceMappingURL=.*\s*$", "\n", library)
    return (
        library + "\n" + widget.replace("__GENEVIEWER_NEW_CHART__", names["newChart"])
    )


class GCWidget(anywidget.AnyWidget):
    """A chart as a Jupyter widget.

    A cell that ends with the widget draws the chart. The widget's
    attributes hold what happened in the chart last (see the events of
    this module), `zoom_to` moves the zoom slider, and setting `chart`
    draws another chart in the same place.

    Parameters
    ----------
    chart : GC
        The chart to draw.

    Attributes
    ----------
    chart : GC
        The chart drawn. Set it to draw another one in its place.
    spec : dict
        The spec of the chart, as the browser has it.
    click, dblclick, select, selection, drag, moves, zoom, rendered : dict, list or None
        The last event of that name, ``None`` until it happens.

    Examples
    --------
    >>> from geneviewer import GC, datasets
    >>> from geneviewer.jupyter import GCWidget
    >>> genes = datasets.load("ophA_clusters")
    >>> widget = GCWidget(GC(genes, cluster="cluster", group="class").zoom())
    >>> widget.spec["fields"]
    {'cluster': 'cluster', 'group': 'class'}
    >>> widget.zoom_to([0, 0.5])
    """

    _esm = _module()

    spec = traitlets.Dict().tag(sync=True)
    click = traitlets.Any(None, allow_none=True)
    dblclick = traitlets.Any(None, allow_none=True)
    select = traitlets.Any(None, allow_none=True)
    selection = traitlets.Any(None, allow_none=True)
    drag = traitlets.Any(None, allow_none=True)
    moves = traitlets.Any(None, allow_none=True)
    zoom = traitlets.Any(None, allow_none=True)
    rendered = traitlets.Any(None, allow_none=True)

    def __init__(self, chart: GC, **kwargs: Any) -> None:
        _check_chart(chart)
        super().__init__(spec=chart.to_spec(), **kwargs)
        self._chart = chart
        self._callbacks: dict[str, list[Callable[[Any], None]]] = {}
        self.on_msg(self._message)

    @property
    def chart(self) -> GC:
        """The chart drawn; set it to draw another one in its place."""
        return self._chart

    @chart.setter
    def chart(self, chart: GC) -> None:
        _check_chart(chart)
        self._chart = chart
        self.spec = chart.to_spec()

    def _message(self, _widget: Any, content: Any, _buffers: Any) -> None:
        # An event of the chart, sent by geneviewer-jupyter.js
        if not isinstance(content, dict) or content.get("type") != "event":
            return
        name = content.get("name")
        if name not in EVENTS:
            return
        detail = content.get("detail")
        setattr(self, name, detail)
        for callback in list(self._callbacks.get(name, ())):
            callback(detail)

    def on(
        self, event: str, callback: Callable[[Any], None], remove: bool = False
    ) -> None:
        """Call a function at each event of a kind.

        Unlike ``observe()`` on the attribute, the function is also called
        when the same event happens twice, such as a second click on the
        same gene.

        Parameters
        ----------
        event : str
            The event: "click", "dblclick", "select", "selection", "drag",
            "moves", "zoom" or "rendered".
        callback : callable
            Called with the event's payload.
        remove : bool
            Whether to stop calling the function instead. Default False.

        Raises
        ------
        ValueError
            When there is no such event.
        """
        if event not in EVENTS:
            raise ValueError(
                f"There is no event {event!r}. The events are: "
                + ", ".join(EVENTS)
                + "."
            )
        callbacks = self._callbacks.setdefault(event, [])
        if remove:
            if callback in callbacks:
                callbacks.remove(callback)
        elif callback not in callbacks:
            callbacks.append(callback)

    def zoom_to(self, range: Any = None) -> None:  # noqa: A002 (R's argument name)
        """Set the zoom window of the chart, R's ``GC_zoom_to()``.

        Moves the zoom slider of ``.zoom()`` without drawing the chart
        again. With ``.events(zoom=True)``, the chart reports the new
        window in ``zoom``, with ``action`` ``"set"``.

        Parameters
        ----------
        range : sequence of two numbers or None
            The window, as ``[start, end]`` in the units of the ``mode`` of
            ``.zoom()``: shares of the clusters (``"relative"``) or base
            pairs (``"absolute"``). ``None`` shows everything. The
            ``range`` of ``zoom`` can be passed back as it is.

        Raises
        ------
        ValueError
            When ``range`` is not valid.
        """
        self.send({"type": "zoom", "range": zoom_window(range)})


def _check_chart(chart: Any) -> None:
    if not isinstance(chart, GC):
        raise TypeError(
            f"A GCWidget draws a chart made with GC(), not {type(chart).__name__}."
        )
