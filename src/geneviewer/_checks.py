"""Argument checks shared by the chart methods and themes.

Ported from R: ``check_align``, ``check_spacing``, ``break_padding_bp`` and
``get_relative_height`` (``R/utils.R``), ``check_group_style`` and
``check_pattern`` (``R/chart_genes.R``), ``check_toolbar`` and
``check_zoom`` (``R/chart_interaction.R``). The messages name the Python
arguments.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Iterable, Mapping
from typing import Any

from . import _r
from .spec import TOOLBAR_BUTTONS, camel_to_snake


def quoted(choices: Iterable[str]) -> str:
    """Return choices for an error message: ``"a", "b", "c"``.

    Parameters
    ----------
    choices : Iterable of str
        The choices.

    Returns
    -------
    str
        The quoted choices.
    """
    return ", ".join(f'"{choice}"' for choice in choices)


def is_str(value: Any) -> bool:
    """Return whether a value is a single string (R's character of length 1).

    Parameters
    ----------
    value : Any
        The value.

    Returns
    -------
    bool
        True for a ``str``.
    """
    return isinstance(value, str)


def is_strings(value: Any) -> bool:
    """Return whether a value is a string or a non-empty list of strings.

    Parameters
    ----------
    value : Any
        The value.

    Returns
    -------
    bool
        True for R's character vectors.
    """
    values = _r.as_vector(value)
    return len(values) > 0 and all(isinstance(v, str) for v in values)


def check_align(align: Any, per_cluster: bool = False) -> None:
    """Check an ``align`` argument: "left", "center" or "right".

    Parameters
    ----------
    align : str, list of str or None
        The argument.
    per_cluster : bool
        Whether ``align`` may have one value per cluster.
    """
    if align is None:
        return
    values = _r.as_vector(align)
    valid = is_strings(align) and all(v in ("left", "center", "right") for v in values)
    if not valid or (not per_cluster and len(values) != 1):
        raise ValueError('align must be "left", "center" or "right"')


def check_spacing(spacing: Any, per_cluster: bool = False) -> None:
    """Check a ``spacing`` argument: pixels between a title and its subtitle.

    Parameters
    ----------
    spacing : float, list of float or None
        The argument.
    per_cluster : bool
        Whether ``spacing`` may have one value per cluster.
    """
    if spacing is None:
        return
    values = _r.as_vector(spacing)
    valid = len(values) > 0 and all(
        _r.is_number(v) and math.isfinite(v) for v in values
    )
    if not valid or (not per_cluster and len(values) != 1):
        raise ValueError(
            "spacing must be numbers"
            if per_cluster
            else "spacing must be a single number"
        )


def check_group_style(group_style: Any, name: str = "group_style") -> None:
    """Check that a group style is a dict of style dicts, one per group.

    Parameters
    ----------
    group_style : dict or None
        The argument.
    name : str
        The argument's name, for the message.
    """
    if group_style is None:
        return
    if (
        not isinstance(group_style, Mapping)
        or any(str(key) == "" for key in group_style)
        or not all(isinstance(value, Mapping) for value in group_style.values())
    ):
        raise ValueError(f"{name} must be a dict of style dicts, one per group")


def check_style_by(state: Mapping[str, Any], style_by: Any) -> None:
    """Check that ``style_by`` is a single column of the data, or "geneTrack".

    Parameters
    ----------
    state : Mapping
        From ``data.chart_fields()``.
    style_by : str or None
        The argument.
    """
    if style_by is None:
        return
    columns = state["columns"] or []
    if not isinstance(style_by, str) or not (
        style_by in columns or style_by == "geneTrack"
    ):
        raise ValueError(
            'style_by must be a single column name of the data, or "geneTrack"'
        )


PATTERNS = (
    "none",
    "solid",
    "horizontal",
    "vertical",
    "diagonal",
    "antidiagonal",
    "grid",
    "crosshatch",
    "dots",
    "zigzag",
    "checkerboard",
    "rings",
    "bricks",
    "plus",
    "wideStripes",
    "inverseDots",
)


def check_pattern(pattern: Mapping[str, Any]) -> None:
    """Check the patterns of ``pattern()`` and ``theme(pattern=...)``.

    Parameters
    ----------
    pattern : Mapping
        The arguments of ``pattern()``.
    """
    patterns = pattern.get("patterns")
    if isinstance(patterns, Mapping):
        patterns = list(patterns.values())
    unknown = [p for p in _r.as_vector(patterns) if p not in PATTERNS]
    if unknown:
        raise ValueError(
            f"Unknown pattern(s): {', '.join(map(str, unknown))}. "
            f"Available patterns: {', '.join(PATTERNS)}"
        )


def _python_name(name: str) -> str:
    return camel_to_snake(name)


def check_toolbar(
    options: Mapping[str, Any], label: Callable[[str], str] = _python_name
) -> None:
    """Check the options of the toolbar.

    Parameters
    ----------
    options : Mapping
        The options given, by their camelCase names.
    label : callable
        How an option is named in a message.
    """

    def one_of(name: str, choices: tuple[str, ...]) -> None:
        value = options.get(name)
        if value is not None and not (is_str(value) and value in choices):
            raise ValueError(f"'{label(name)}' must be one of {quoted(choices)}.")

    one_of("show", ("hover", "always", "never"))
    one_of("position", ("topright", "topleft", "bottomright", "bottomleft"))
    buttons = options.get("buttons")
    if buttons is not None and not (
        is_strings(buttons) and all(b in TOOLBAR_BUTTONS for b in _r.as_vector(buttons))
    ):
        raise ValueError(f"'{label('buttons')}' must be of {quoted(TOOLBAR_BUTTONS)}.")
    size = options.get("size")
    if size is not None and not (
        _r.is_number(size) and not _r.is_na(size) and size > 0
    ):
        raise ValueError(
            f"'{label('size')}' must be a single positive number of pixels."
        )
    colours = ("color", "hoverColor", "tooltipBackground", "tooltipColor")
    for name in (*colours, "tooltipBorder", "tooltipShadow", "fontFamily"):
        value = options.get(name)
        if value is not None and not is_str(value):
            kind = "colour" if name in colours else "string"
            raise ValueError(f"'{label(name)}' must be a single {kind}.")
    radius = options.get("tooltipRadius")
    if radius is not None and not (
        is_str(radius)
        or (_r.is_number(radius) and not _r.is_na(radius) and radius >= 0)
    ):
        raise ValueError(
            f"'{label('tooltipRadius')}' must be a number of pixels or a CSS length."
        )
    icons = options.get("icons")
    if icons is not None and not (
        isinstance(icons, Mapping)
        and all(name in TOOLBAR_BUTTONS for name in icons)
        and all(is_str(icon) and icon != "" for icon in icons.values())
    ):
        raise ValueError(
            f"'{label('icons')}' must name Font Awesome icons for "
            f'{quoted(TOOLBAR_BUTTONS)}, e.g. {{"reset": "arrows-to-circle"}}.'
        )


def check_zoom(options: Mapping[str, Any], mode: str | None = None) -> None:
    """Check the options of the zoom slider.

    Parameters
    ----------
    options : Mapping
        The options given (not ``None``), by their Python names.
    mode : str, optional
        The mode that ``range`` and ``min`` are in: of these options, else
        of an earlier ``zoom()``.
    """
    if mode is None:
        mode = options.get("mode")
    for name in ("show", "overview"):
        value = options.get(name)
        if value is not None and not _r.is_bool(value):
            raise ValueError(f"'{name}' must be True or False.")

    def one_of(name: str, choices: tuple[str, ...]) -> None:
        value = options.get(name)
        if value is not None and not (is_str(value) and value in choices):
            raise ValueError(f"'{name}' must be one of {quoted(choices)}.")

    one_of("mode", ("relative", "absolute"))
    one_of("update", ("live", "release"))
    window = options.get("range")
    if window is not None:
        values = _r.as_vector(window)
        if (
            len(values) != 2
            or not all(_r.is_number(v) and math.isfinite(v) for v in values)
            or values[0] == values[1]
        ):
            raise ValueError("'range' must be two different numbers: [start, end].")
        if mode != "absolute" and (min(values) < 0 or max(values) > 1):
            raise ValueError(
                "A relative 'range' is shares of the clusters, between 0 and 1, "
                'e.g. [0.25, 0.5]. Use mode="absolute" for base pairs.'
            )
    for name in ("min", "height"):
        value = options.get(name)
        if value is not None and not (
            _r.is_number(value) and not _r.is_na(value) and value > 0
        ):
            raise ValueError(f"'{name}' must be a single positive number.")
    if options.get("min") is not None and mode != "absolute" and options["min"] > 1:
        raise ValueError("A relative 'min' is a share of the slider, between 0 and 1.")
    color = options.get("color")
    if color is not None and not is_str(color):
        raise ValueError("'color' must be a single colour.")
    for name in ("style", "window_style", "handle_style"):
        value = options.get(name)
        if value is not None and not isinstance(value, Mapping):
            raise ValueError(f"'{name}' must be a dict of CSS styles.")


def break_padding_bp(padding: str) -> float:
    """Return the base pairs of a padding given as a string, e.g. "500bp".

    Parameters
    ----------
    padding : str
        The padding.

    Returns
    -------
    float
        The base pairs.
    """
    if not re.fullmatch(r"\s*[0-9]*\.?[0-9]+\s*bp\s*", padding):
        raise ValueError(
            "scale_break_padding must be a percentage, e.g. 1, "
            'or base pairs, e.g. "500bp"'
        )
    return float(re.sub(r"\s*bp\s*$", "", padding))


def get_relative_height(base_height: Any, relative_height: Any) -> Any:
    """Return a height relative to a base height (R's ``get_relative_height``).

    Parameters
    ----------
    base_height : float or str
        The base height, e.g. ``"400px"``.
    relative_height : float or str
        A height in pixels or a percentage of the base height.

    Returns
    -------
    float or str
        The height in pixels, or the percentage when both are percentages.
    """
    if (
        is_str(base_height)
        and base_height.endswith("%")
        and is_str(relative_height)
        and relative_height.endswith("%")
    ):
        return relative_height
    if is_str(base_height) and base_height.endswith("px"):
        base_height = float(base_height.replace("px", ""))
    if is_str(relative_height) and relative_height.endswith("%"):
        return base_height * float(relative_height.replace("%", "")) / 100
    if is_str(relative_height) and relative_height.endswith("px"):
        relative_height = relative_height.replace("px", "")
    return float(relative_height)


def zoom_window(range: Any) -> list[float] | None:  # noqa: A002 (R's argument name)
    """Check a zoom window, R's checks of ``GC_zoom_to()``.

    Parameters
    ----------
    range : sequence of two numbers or None
        The window.

    Returns
    -------
    list of float or None
        The window, sorted, or None for everything.

    Raises
    ------
    ValueError
        When it is not None or two numbers.
    """
    if range is None:
        return None
    values = _r.as_vector(range)
    if len(values) != 2 or not all(
        _r.is_number(v) and math.isfinite(v) for v in values
    ):
        raise ValueError("'range' must be None or two numbers: [start, end].")
    return sorted(float(v) for v in values)
