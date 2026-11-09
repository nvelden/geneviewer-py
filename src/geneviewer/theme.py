"""Built-in themes and palettes, and the checks of theme settings.

Ported from R's ``R/theme.R``. The themes and palettes are those of the
vendored library's ``themes.json``; their option names are the spec's.
"""

from __future__ import annotations

import copy
import functools
import json
from collections.abc import Mapping
from importlib import resources
from typing import Any

from . import _r
from ._checks import check_pattern, check_toolbar, is_str, is_strings
from .spec import THEME_COMPONENTS, camel_case_names


@functools.cache
def _themes_json() -> dict[str, Any]:
    path = resources.files("geneviewer") / "static" / "themes.json"
    return json.loads(path.read_text(encoding="utf-8"))


def geneviewer_theme(name: str = "default") -> dict[str, Any]:
    """Return the settings of a built-in theme.

    To see them, or to start a theme of your own. The options of components
    have the names of the chart spec, e.g. ``cornerRadius`` and
    ``textStyle``; ``GC.theme()`` takes these and the method's argument
    names (``corner_radius``, ``legendTextOptions``).

    Parameters
    ----------
    name : str
        The theme: "default", "classic", "dark", "vintage", "infographic",
        "publication" or "bw".

    Returns
    -------
    dict
        A new copy of the theme's settings, to change and pass to
        ``GC.theme()``, as ``theme=`` or as ``**settings``.

    Examples
    --------
    >>> settings = geneviewer_theme("classic")
    >>> settings["textColor"] = "navy"
    """
    themes = _themes_json()["themes"]
    if not is_str(name) or name not in themes:
        raise ValueError(f"name must be one of: {_quoted(themes)}")
    return copy.deepcopy(themes[name])


def geneviewer_palette(name: str | None = None) -> list[str]:
    """Return the colours of a geneviewer palette, or the names of all palettes.

    Besides the d3 colour schemes (e.g. "schemeSet2"), the gene colours can
    be set by the name of a geneviewer palette, in ``GC.theme(palette=...)``
    and in the ``color_scheme`` of ``GC.genes()``, ``GC.transcript()`` and
    ``GC.color()``:

    - "okabeIto": 8 colours that stay distinct with colour blindness (Okabe
      & Ito, also known as the Wong palette), the palette of the
      "publication" theme.
    - "classic": 9 clear colours, the palette of the "classic" theme.
    - "vivid": 9 bright colours for dark backgrounds.
    - "dark": 11 soft colours for dark backgrounds, the palette of the
      "dark" theme.
    - "vintage": 10 muted colours, the palette of the "vintage" theme.
    - "infographic": 15 strong colours, the palette of the "infographic"
      theme.
    - "greys": 6 greys, from dark to light.
    - "blues": 7 blues, from dark to light.

    Parameters
    ----------
    name : str, optional
        The palette. Default: the names of all palettes.

    Returns
    -------
    list of str
        The colours of the palette, or the names of all palettes.
    """
    palettes = _themes_json()["palettes"]
    if name is None:
        return list(palettes)
    if not is_str(name) or name not in palettes:
        raise ValueError(f"name must be one of: {_quoted(palettes)}")
    return list(palettes[name])


def _quoted(names: Any) -> str:
    return ", ".join(f'"{name}"' for name in names)


# Global settings of a theme, and the annotation types it styles
THEME_GLOBAL_SETTINGS = (
    "textScale",
    "fontFamily",
    "textColor",
    "lineColor",
    "background",
    "palette",
    "pattern",
)
THEME_ANNOTATION_TYPES = (
    "text",
    "line",
    "arrow",
    "textMarker",
    "symbol",
    "rectangle",
    "promoter",
    "terminator",
    "arc",
)
PATTERN_ARGS = ("patterns", "color", "background", "size", "strokeWidth")


def resolve_theme(theme: Any) -> Any:
    """Check the ``theme`` of ``GC.theme()``: a built-in theme or settings.

    Parameters
    ----------
    theme : str or dict
        The name of a built-in theme, or theme settings.

    Returns
    -------
    dict
        The theme's settings.
    """
    if is_str(theme):
        return geneviewer_theme(theme)
    if not isinstance(theme, Mapping):
        raise ValueError(
            "theme must be the name of a built-in theme or a dict of theme settings"
        )
    check_theme(theme)
    return theme


def check_theme(settings: Mapping[str, Any]) -> None:
    """Check theme settings, by their camelCase names.

    Parameters
    ----------
    settings : Mapping
        The settings.
    """
    if not settings:
        return
    valid = (*THEME_GLOBAL_SETTINGS, *THEME_COMPONENTS, "annotations")
    unknown = [name for name in settings if name not in valid]
    if unknown:
        raise ValueError(
            f"Unknown theme setting: {', '.join(map(str, unknown))}. "
            f"Use one of: {', '.join(valid)}"
        )
    for name in (*THEME_COMPONENTS, "annotations"):
        value = settings.get(name)
        # None removes the component from the theme
        if value is not None and not isinstance(value, Mapping):
            raise ValueError(f"Theme setting '{name}' must be a dict of options")
    scale = settings.get("textScale")
    if scale is not None and not (
        _r.is_number(scale) and not _r.is_na(scale) and scale > 0
    ):
        raise ValueError("text_scale must be a single positive number")
    for name in ("fontFamily", "textColor", "lineColor", "background"):
        if name in settings and not is_str(settings[name]):
            raise ValueError(f"Theme setting '{name}' must be a single string")
    palette = settings.get("palette")
    if palette is not None and not is_strings(palette):
        raise ValueError(
            "palette must be the name of a d3 colour scheme or a list of colours"
        )
    pattern = settings.get("pattern")
    if pattern is not None:
        if not isinstance(pattern, Mapping) or not all(
            name in PATTERN_ARGS for name in pattern
        ):
            raise ValueError(
                "pattern must be a dict with the arguments of pattern(): "
                + ", ".join(PATTERN_ARGS)
            )
        check_pattern(pattern)
    if settings.get("toolbar") is not None:
        check_toolbar(camel_case_names(settings["toolbar"]), label=str)
    annotations = settings.get("annotations") or {}
    unknown = [name for name in annotations if name not in THEME_ANNOTATION_TYPES]
    if unknown:
        raise ValueError(
            f"Unknown annotation type in theme: {', '.join(map(str, unknown))}. "
            f"Use one of: {', '.join(THEME_ANNOTATION_TYPES)}"
        )
