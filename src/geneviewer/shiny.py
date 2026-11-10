"""Charts in Shiny for Python.

Put a chart in the UI with `output_geneviewer`, draw it with
`render_geneviewer`, and move its zoom slider from the server with
`zoom_to`::

    from shiny import App, ui
    from geneviewer import GC
    from geneviewer.shiny import output_geneviewer, render_geneviewer

    app_ui = ui.page_fluid(output_geneviewer("chart", height="300px"))

    def server(input, output, session):
        @render_geneviewer
        def chart():
            return GC(genes, group="class").labels("name").events(select=True)

    app = App(app_ui, server)

**Events.** The events the chart's ``.events()`` switches on set inputs
named after the output, ``input.<id>_<event>()``:

- ``<id>_click``: the element clicked (``.events(click=True)``, the default).
- ``<id>_dblclick``: the element double-clicked (``dblclick=True``).
- ``<id>_select``: the element selected last (``select=True``).
- ``<id>_selection``: every selected element, as a list (``select=True``).
- ``<id>_drag``: an element moved by dragging, and its offset (``drag=True``).
- ``<id>_moves``: the elements a reset, undo or redo moved (``drag=True``).
- ``<id>_zoom``: the window of the zoom slider (``zoom=True``).
- ``<id>_rendered``: the chart was drawn (``rendered=True``).

Each input is set again when the same event happens twice. The payload of
each event is described at https://nvelden.github.io/geneviewer-docs/api#events.
"""

from __future__ import annotations

from functools import cache
from importlib import resources
from typing import Any

try:
    from htmltools import HTMLDependency, Tag, css, div
    from shiny.module import resolve_id
    from shiny.render.renderer import Jsonifiable, Renderer
    from shiny.session import Session, get_current_session
except ImportError as error:  # pragma: no cover
    raise ImportError(
        "geneviewer.shiny needs Shiny for Python: pip install 'geneviewer-py[shiny]'"
    ) from error

from . import _r
from ._checks import zoom_window
from .chart import GC

__all__ = ["EVENTS", "output_geneviewer", "render_geneviewer", "zoom_to"]

#: The events the binding sends to Shiny as ``input.<id>_<event>()``, as R's
#: binding does (not ``dragmove``).
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
def _dependency() -> HTMLDependency:
    static = resources.files("geneviewer") / "static"
    version = (static / "VERSION").read_text(encoding="utf-8").strip()
    return HTMLDependency(
        "geneviewer",
        version,
        source={"subdir": str(static)},
        script=[{"src": "geneviewer.min.js"}, {"src": "geneviewer-shiny.js"}],
        stylesheet={"href": "geneviewer.css"},
    )


def _css_size(value: str | float | int) -> str:
    # As R's htmlwidgets: a number is pixels
    if _r.is_number(value):
        return f"{value}px"
    return str(value)


def output_geneviewer(
    id: str, width: str | float = "100%", height: str | float = "400px"
) -> Tag:
    """Create a chart output in the UI, R's ``GC_chartOutput()``.

    Draw the chart in the server with `render_geneviewer`, under the same
    id. A fixed pixel width set by ``GC(width=)`` or ``.grid(width=)`` takes
    precedence over the output's width; put the output in a container with
    ``overflow-x: auto`` to scroll a chart wider than the page.

    Parameters
    ----------
    id : str
        The output id.
    width, height : str or float
        The size of the output, as CSS (``"100%"``, ``"400px"``, ``"auto"``),
        or a number of pixels.

    Returns
    -------
    htmltools.Tag
        The output, with the library's JavaScript and styles.
    """
    return div(
        _dependency(),
        id=resolve_id(id),
        class_="geneviewer-output",
        style=css(width=_css_size(width), height=_css_size(height)),
    )


class render_geneviewer(Renderer[GC]):  # noqa: N801 (Shiny's renderer naming)
    """Draw a chart in an output, R's ``renderGC_chart()``.

    Decorate a function that returns a chart (`geneviewer.GC`); its name is
    the output id. Each time the function runs again, the output draws the
    new chart in place of the old one. Returning ``None`` removes the
    chart.

    In Shiny Express, the decorated function also places the output, as
    ``output_geneviewer(id)``.
    """

    def auto_output_ui(self) -> Tag:
        """Return the output for Shiny Express.

        Returns
        -------
        htmltools.Tag
            ``output_geneviewer()`` with this renderer's id.
        """
        return output_geneviewer(self.output_id)

    async def transform(self, value: GC) -> Jsonifiable:
        """Turn the chart into the value the output binding draws.

        Parameters
        ----------
        value : GC
            The chart the function returned.

        Returns
        -------
        dict
            ``{"spec": value.to_spec()}``.
        """
        if not isinstance(value, GC):
            raise TypeError(
                "A render_geneviewer function must return a chart made with GC(), "
                f"not {type(value).__name__}."
            )
        return {"spec": value.to_spec()}


def zoom_to(
    id: str,
    range: Any = None,  # noqa: A002 (R's argument name)
    session: Session | None = None,
) -> None:
    """Set the zoom window of a chart in a running app, R's ``GC_zoom_to()``.

    Moves the zoom slider of ``.zoom()`` without drawing the chart again
    from Python. With ``.events(zoom=True)``, the chart reports the new
    window as ``input.<id>_zoom()``, with ``action`` ``"set"``.

    Parameters
    ----------
    id : str
        The output id of the chart, as in `output_geneviewer`.
    range : sequence of two numbers or None
        The window, as ``[start, end]`` in the units of the ``mode`` of
        ``.zoom()``: shares of the clusters (``"relative"``) or base pairs
        (``"absolute"``). ``None`` shows everything. The ``range`` of
        ``input.<id>_zoom()`` can be passed back as it is.
    session : shiny.Session, optional
        The session; by default the current one.

    Raises
    ------
    ValueError
        When ``id`` or ``range`` is not valid.
    RuntimeError
        Outside a Shiny session.
    """
    if not isinstance(id, str):
        raise ValueError("'id' must be the output id of the chart.")
    window = zoom_window(range)
    session = session or get_current_session()
    if session is None:
        raise RuntimeError("zoom_to() works in a Shiny session.")
    # Sent as Shiny's own helpers send their custom messages (e.g.
    # ui.update_dark_mode()): at once, without making the caller await it
    session._send_message_sync(
        {"custom": {"geneviewer-zoom": {"id": session.ns(id), "range": window}}}
    )
