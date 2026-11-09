"""The chart spec, ported from R's ``R/spec.R``.

``GC()`` starts the spec, and every method adds the arguments it was called
with, under its key and with the names of the spec (spec v1,
``schema/NAMING.md`` in geneviewer-js). Arguments are recorded as given: no
defaults are filled in and nothing is recycled or expanded. The functions
here have the names of their R counterparts without the ``spec_`` prefix
(``spec_args`` is ``args``, ``add_spec`` is ``add``), so the two can be read
side by side.

The functions do not change the spec they are given: each returns a new
one.
"""

from __future__ import annotations

import copy
import re
from collections.abc import Iterable, Mapping
from typing import Any

import pandas as pd

from . import _r


class _Missing:
    """The default of every argument: not given."""

    def __repr__(self) -> str:
        return "<default>"

    def __bool__(self) -> bool:
        return False


_MISSING: Any = _Missing()


class Args(dict):
    """The arguments given in a call, for the spec.

    A dict of argument names (Python names) and values. ``dots`` names the
    arguments passed through ``**options``, which keep their names apart
    from camelCase (R's ``attr(args, "dots")``).
    """

    def __init__(self, *args: Any, dots: Iterable[str] = (), **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.dots: set[str] = set(dots)

    def copy(self) -> Args:
        """Return a shallow copy, with its ``dots``."""
        return Args(self, dots=self.dots)


class Call:
    """A call of a chart method: what was given, and the value of each argument.

    Parameters
    ----------
    spec : Args
        The arguments given, for the spec (see `args`).
    values : dict
        Every argument's value: as given, or its default.
    options : dict
        The arguments passed through ``**options``, by their spec name.
    """

    def __init__(self, spec: Args, values: dict[str, Any], options: dict[str, Any]):
        self.spec = spec
        self.values = values
        self.options = options

    def __getitem__(self, name: str) -> Any:
        return self.values[name]


def snake_to_camel(name: str) -> str:
    """Return an option name of the methods as in the spec (R's ``snake_to_camel``).

    ``normal_color`` becomes ``normalColor``; camelCase names stay as they are.

    Parameters
    ----------
    name : str
        The name.

    Returns
    -------
    str
        The name in camelCase.
    """
    return re.sub(r"_([a-z])", lambda match: match.group(1).upper(), name)


def camel_to_snake(name: str) -> str:
    """Return the snake_case name of a camelCase name, e.g. ``itemStyle``.

    Parameters
    ----------
    name : str
        The name.

    Returns
    -------
    str
        The name in snake_case: ``item_style``, and ``row_id`` for ``rowID``.
    """
    return re.sub(r"(?<=[a-z0-9])([A-Z]+)", lambda m: "_" + m.group(1).lower(), name)


def camel_case_names(options: Any) -> Any:
    """Return a dict with its names in camelCase (R's ``camel_case_names``).

    Parameters
    ----------
    options : Any
        A dict of options; other values are returned as they are.

    Returns
    -------
    Any
        The dict with camelCase names.
    """
    if isinstance(options, Mapping):
        return {snake_to_camel(name): value for name, value in options.items()}
    return options


def _is_empty_list(value: Any) -> bool:
    # R's list(): an empty dict, list or tuple
    return isinstance(value, (dict, list, tuple)) and len(value) == 0


def args(values: Mapping[str, Any], defaults: Mapping[str, Any]) -> Call:
    """Return the arguments given in a call, for the spec (R's ``spec_args``).

    Call it first in a method, with ``locals()``, before any argument is
    changed.

    Parameters
    ----------
    values : Mapping
        The arguments of the call, from ``locals()``: ``_MISSING`` for those
        not given, and ``options`` for those passed through ``**options``.
        ``self`` is left out.
    defaults : Mapping
        The default of each argument (``None`` for R's ``NULL``).

    Returns
    -------
    Call
        ``spec``: the arguments given, without empty lists. A ``None`` is
        kept only where the argument's default is not ``None`` (e.g.
        ``pattern(background=None)``); a ``None`` passed through
        ``**options`` is left out. ``values``: the value of every argument.
        A camelCase name of an argument in ``**options`` (``itemStyle``) is
        taken as the argument (``item_style``).
    """
    values = {name: value for name, value in values.items() if name != "self"}
    options = dict(values.pop("options", None) or {})
    for name in list(options):
        snake = camel_to_snake(name)
        if snake != name and snake in defaults:
            if values.get(snake, _MISSING) is not _MISSING:
                raise TypeError(f"got both '{snake}' and '{name}'")
            values[snake] = options.pop(name)

    spec = Args()
    for name, value in values.items():
        if value is _MISSING:
            continue
        if value is None:
            if defaults.get(name) is not None:
                spec[name] = None
            continue
        if not _is_empty_list(value):
            spec[name] = value
    for name, value in options.items():
        if value is None or _is_empty_list(value):
            continue
        spec[name] = value
        spec.dots.add(name)

    resolved = {
        name: defaults.get(name) if value is _MISSING else value
        for name, value in values.items()
    }
    return Call(spec, resolved, {snake_to_camel(k): v for k, v in options.items()})


def entry(value: Mapping[str, Any]) -> dict[str, Any]:
    """Return an entry of the spec: a dict, also when empty (R's ``spec_entry``).

    Parameters
    ----------
    value : Mapping
        The entry.

    Returns
    -------
    dict
        The entry as a dict.
    """
    return dict(value)


def is_entry(value: Any) -> bool:
    """Return whether a spec value is one entry, not an array of entries.

    Parameters
    ----------
    value : Any
        The value of a spec key.

    Returns
    -------
    bool
        True for a dict.
    """
    return isinstance(value, Mapping)


def names(
    arg_names: Iterable[str],
    rename: Mapping[str, str] | None = None,
    clusters: bool = True,
) -> list[str]:
    """Return the spec names of arguments (R's ``spec_names``).

    Parameters
    ----------
    arg_names : Iterable of str
        The argument names, in Python's snake_case.
    rename : Mapping, optional
        Spec names of single arguments, by argument name.
    clusters : bool
        Whether ``cluster`` becomes ``clusters``.

    Returns
    -------
    list of str
        The spec names.
    """
    rename = rename or {}
    out = []
    for name in arg_names:
        if name in rename:
            out.append(rename[name])
        elif clusters and name == "cluster":
            out.append("clusters")
        else:
            out.append(snake_to_camel(name))
    return out


# The renderer reads scaleBarUnit and labelAdjustmentOptions; spec v1 names
# them scaleBar.unit and labelAdjustment (decided 2026-11-05)
DOT_RENAMES: dict[str, dict[str, str]] = {
    "scaleBar": {"scaleBarUnit": "unit"},
    "labels": {"labelAdjustmentOptions": "labelAdjustment"},
    "links": {"labelAdjustmentOptions": "labelAdjustment"},
}


def dot_names(key: str, arg_names: Iterable[str]) -> list[str]:
    """Return the spec names of options passed through ``**options``.

    They are written in camelCase, and the spec renames some of them
    (``DOT_RENAMES``; R's ``spec_dot_names``).

    Parameters
    ----------
    key : str
        The key of the spec, e.g. ``"scaleBar"``.
    arg_names : Iterable of str
        The names of the options.

    Returns
    -------
    list of str
        The spec names.
    """
    renames = DOT_RENAMES.get(key, {})
    out = []
    for name in arg_names:
        name = snake_to_camel(name)
        out.append(renames.get(name, name))
    return out


# The options of each component in a theme have the names of the component
# in the spec (theme.legend.textStyle, decided 2026-11-05); the settings a
# user passes take the argument names too, camelCased first
THEME_RENAMES: dict[str, dict[str, str]] = {
    "legend": {"legendOptions": "markerStyle", "legendTextOptions": "textStyle"},
    "sequence": {"sequenceStyle": "lineStyle"},
    "scaleBar": {
        "scaleBarLineStyle": "lineStyle",
        "scaleBarTickStyle": "tickStyle",
        "scaleBarUnit": "unit",
    },
    "links": {"showLinks": "show", "labelAdjustmentOptions": "labelAdjustment"},
    "labels": {"labelAdjustmentOptions": "labelAdjustment"},
    "scale": {
        "scaleBreaks": "autoBreaks",
        "scaleBreakThreshold": "autoBreakThreshold",
        "scaleBreakPadding": "breakPadding",
    },
}

# The components of a theme, and its global settings
THEME_COMPONENTS = (
    "chart",
    "grid",
    "cluster",
    "title",
    "clusterTitle",
    "clusterFooter",
    "clusterLabel",
    "sequence",
    "genes",
    "transcript",
    "labels",
    "coordinates",
    "scale",
    "scaleBar",
    "legend",
    "links",
    "tooltip",
    "trackMouse",
    "toolbar",
)


def _copy(value: Any) -> Any:
    # A deep copy; data frames are copied with pandas
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if isinstance(value, Mapping):
        out = {name: _copy(item) for name, item in value.items()}
        if isinstance(value, Args):
            return Args(out, dots=value.dots)
        return out
    if isinstance(value, list):
        return [_copy(item) for item in value]
    return copy.deepcopy(value)


def _as_list(value: Any) -> list[Any]:
    # R's as.list() of a vector: a scalar becomes a list of one
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, pd.Series):
        return value.tolist()
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        result = value.tolist()
        return result if isinstance(result, list) else [result]
    return [value]


def merge_settings(x: Mapping[str, Any], y: Mapping[str, Any], keep_null: bool) -> dict:
    """Merge settings: dicts at every level, other values replace.

    R's ``merge_settings``: dicts are merged at every level; other values,
    lists included, replace the old value; ``None`` removes it, or is kept
    as ``None`` when ``keep_null`` is true.

    Parameters
    ----------
    x : Mapping
        The settings.
    y : Mapping
        The settings merged into them.
    keep_null : bool
        Whether a ``None`` in ``y`` is kept or removes the setting.

    Returns
    -------
    dict
        The merged settings, a new dict.
    """
    out = dict(x)
    for name, value in y.items():
        if value is None:
            if keep_null:
                out[name] = None
            else:
                out.pop(name, None)
        elif isinstance(out.get(name), Mapping) and isinstance(value, Mapping):
            out[name] = merge_settings(out[name], value, keep_null)
        else:
            out[name] = value
    return out


def add(
    spec: Mapping[str, Any],
    key: str,
    call_args: Args,
    rename: Mapping[str, str] | None = None,
    mode: str = "append",
    clusters: bool = True,
) -> dict[str, Any]:
    """Add the arguments of a call to the spec (R's ``add_spec``).

    Parameters
    ----------
    spec : Mapping
        The spec.
    key : str
        The key of the spec, e.g. ``"labels"``.
    call_args : Args
        The arguments, from `args`.
    rename : Mapping, optional
        Spec names of single arguments, by argument name.
    mode : {"append", "items", "merge", "replace"}
        How a call adds to earlier calls: ``"append"`` (an entry; more calls
        make an array of entries), ``"items"`` (always an array),
        ``"merge"`` (merged into the one entry: dicts are merged, other
        values replace the old ones and ``None`` removes them) or
        ``"replace"``.
    clusters : bool
        Whether the ``cluster`` argument becomes ``clusters``.

    Returns
    -------
    dict
        The new spec.
    """
    if mode not in ("append", "items", "merge", "replace"):
        raise ValueError(f"unknown mode: {mode}")
    # Options passed through **options keep their names, in camelCase
    new_args: dict[str, Any] = {}
    for name, value in call_args.items():
        if name in call_args.dots:
            (spec_name,) = dot_names(key, [name])
        else:
            (spec_name,) = names([name], rename, clusters)
        new_args[spec_name] = _copy(value)
    for name in ("clusters", "selection"):
        if name in new_args:
            new_args[name] = _as_list(new_args[name])
    new_args = entry(new_args)

    old = spec.get(key)
    if mode == "append":
        if old is None:
            new: Any = new_args
        else:
            new = ([old] if is_entry(old) else list(old)) + [new_args]
    elif mode == "items":
        new = list(old or []) + [new_args]
    elif mode == "merge":
        new = new_args if old is None else merge_settings(old, new_args, False)
    else:
        new = new_args
    out = dict(spec)
    out[key] = new
    return out


# The buttons of the toolbar, in their order (as in the library)
TOOLBAR_BUTTONS = ("undo", "redo", "reset", "zoom", "save")


def start(data: Any, call_args: Mapping[str, Any]) -> dict[str, Any]:
    """Start the spec of a chart (R's ``start_spec``).

    Parameters
    ----------
    data : Any
        The data of the chart, as given.
    call_args : Mapping
        The arguments of ``GC()``, from `args`.

    Returns
    -------
    dict
        The spec.
    """
    spec: dict[str, Any] = {"version": 1, "data": data}
    fields = {
        name: call_args[name]
        for name in ("start", "end", "cluster", "group", "strand")
        if name in call_args
    }
    if fields:
        spec["fields"] = fields
    for name in ("width", "height"):
        if call_args.get(name) is not None:
            spec[name] = call_args[name]
    if call_args.get("style"):
        spec["style"] = _copy(call_args["style"])
    # save_button=False is the toolbar without the save button
    if call_args.get("save_button") is False:
        spec["toolbar"] = {"buttons": [b for b in TOOLBAR_BUTTONS if b != "save"]}
    return spec


def flag(call_args: Args, name: str) -> Args:
    """Record a logical given for a column argument as ``show`` (R's ``spec_flag``).

    E.g. ``labels(False)``.

    Parameters
    ----------
    call_args : Args
        The arguments, from `args`.
    name : str
        The column argument.

    Returns
    -------
    Args
        The arguments, with ``show`` instead of the logical.
    """
    value = call_args.get(name)
    if _r.is_bool(value):
        call_args = call_args.copy()
        del call_args[name]
        call_args["show"] = bool(value)
    return call_args


def add_zoom(spec: Mapping[str, Any], call_args: Args) -> dict[str, Any]:
    """Add a ``zoom()`` call to the spec (R's ``add_spec_zoom``).

    Calls are merged; a new mode does not keep the range or minimum of the
    old one.

    Parameters
    ----------
    spec : Mapping
        The spec.
    call_args : Args
        The arguments, from `args`.

    Returns
    -------
    dict
        The new spec.
    """
    if call_args.get("range") is not None:
        call_args = call_args.copy()
        call_args["range"] = _as_list(call_args["range"])
    old = spec.get("zoom")
    mode = call_args.get("mode")
    if old is not None and mode is not None and mode != old.get("mode"):
        old = {k: v for k, v in old.items() if k not in ("range", "min")}
        spec = {**spec, "zoom": entry(old)}
    return add(spec, "zoom", call_args, mode="merge")


def add_element_style(
    spec: Mapping[str, Any], call_args: Args, component: str, reset: bool
) -> dict[str, Any]:
    """Add an ``element_style()`` call to the spec (R's ``add_spec_element_style``).

    The spec has one entry per call, in order. ``reset=True`` removes the
    earlier entries for the element instead.

    Parameters
    ----------
    spec : Mapping
        The spec.
    call_args : Args
        The arguments, from `args`.
    component : str
        The component of the element.
    reset : bool
        Whether to remove the element's entries.

    Returns
    -------
    dict
        The new spec.
    """
    call_args = call_args.copy()
    call_args.pop("reset", None)
    if reset:
        ids = {"rowID": "row_id", "end": "end", "cluster": "cluster"}
        ids.update(part="part", group="group")
        target = {name: call_args.get(arg) for name, arg in ids.items()}

        def as_character(value: Any) -> list[str | None]:
            return [_r.as_character(v) for v in _r.as_vector(value)]

        def same(item: Mapping[str, Any]) -> bool:
            if item.get("component", "genes") != component:
                return False
            return all(
                as_character(item.get(name)) == as_character(target[name])
                for name in ids
            )

        items = [item for item in spec.get("elementStyle") or [] if not same(item)]
        out = dict(spec)
        if items:
            out["elementStyle"] = items
        else:
            out.pop("elementStyle", None)
        return out
    return add(
        spec,
        "elementStyle",
        call_args,
        rename={"row_id": "rowID"},
        mode="items",
        clusters=False,
    )


def add_theme(spec: Mapping[str, Any], call_args: Mapping[str, Any]) -> dict[str, Any]:
    """Add a ``theme()`` call to the spec (R's ``add_spec_theme``).

    A call with a theme starts again from it; a call without one changes the
    theme set before, merged as in ``theme()``. The options of components
    are named as in the spec's components (``legendTextOptions`` becomes
    ``textStyle``, ``corner_radius`` becomes ``cornerRadius``).

    Parameters
    ----------
    spec : Mapping
        The spec.
    call_args : Mapping
        The arguments, from `args`, with camelCase names.

    Returns
    -------
    dict
        The new spec.
    """
    settings = {name: _copy(value) for name, value in call_args.items()}
    base = settings.pop("theme", None)
    settings = theme_names(settings)
    if isinstance(base, Mapping):
        # A dict of theme settings is the theme itself
        new = merge_settings(theme_names(base), settings, keep_null=False)
    elif base is not None:
        new = {"base": base, **settings}
    else:
        old = spec.get("theme")
        if isinstance(old, str):
            old = {"base": old}
        new = merge_settings(old or {}, settings, keep_null=True)
    return {**spec, "theme": entry(new)}


def theme_names(settings: Mapping[str, Any]) -> dict[str, Any]:
    """Return theme settings with their component options named as in the spec.

    R's ``spec_theme_names``: options of components in camelCase, renamed by
    ``THEME_RENAMES``; an option given under both names, e.g. a built-in
    theme's ``textStyle`` and a ``legendTextOptions`` added to it, becomes
    one option, the later one merged in.

    Parameters
    ----------
    settings : Mapping
        Theme settings.

    Returns
    -------
    dict
        The settings, a new dict.
    """
    out = dict(settings)
    for name in out:
        value = out[name]
        if value is None or name not in (*THEME_COMPONENTS, "annotations"):
            continue
        if name == "annotations":
            if isinstance(value, Mapping):
                out[name] = {k: camel_case_names(v) for k, v in value.items()}
            continue
        options = camel_case_names(value)
        if not isinstance(options, Mapping):
            continue
        renames = THEME_RENAMES.get(name, {})
        merged: dict[str, Any] = {}
        for option, option_value in options.items():
            merged = merge_settings(
                merged, {renames.get(option, option): option_value}, keep_null=True
            )
        out[name] = merged
    return out


# Columns sent to the browser only when an option names them, as in a
# tooltip formatter "{translation}": they hold long sequences the chart does
# not draw
SEQUENCE_COLUMNS = ("translation", "sequence")


def drop_sequence_columns(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Return the spec without the sequence columns no option names.

    R's ``drop_sequence_columns``.

    Parameters
    ----------
    spec : Mapping
        The spec.

    Returns
    -------
    dict
        The spec, with the data without those columns.
    """
    data = spec.get("data")
    if not isinstance(data, pd.DataFrame):
        return dict(spec)
    strings = option_strings({k: v for k, v in spec.items() if k != "data"})
    unused = [
        column
        for column in SEQUENCE_COLUMNS
        if column in data.columns and not any(column in s for s in strings)
    ]
    out = dict(spec)
    if unused:
        out["data"] = data.drop(columns=unused)
    return out


def option_strings(value: Any) -> list[str]:
    """Return all strings of the options, outside data frames.

    Parameters
    ----------
    value : Any
        Options.

    Returns
    -------
    list of str
        The strings.
    """
    if isinstance(value, pd.DataFrame):
        return []
    if isinstance(value, Mapping):
        return [s for item in value.values() for s in option_strings(item)]
    if isinstance(value, (list, tuple)):
        return [s for item in value for s in option_strings(item)]
    return [value] if isinstance(value, str) else []
