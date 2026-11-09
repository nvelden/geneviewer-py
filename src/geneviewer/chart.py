"""The ``GC`` chart: one method per R ``GC_*`` function.

Each method checks its arguments as the R function does and records the
arguments it was given in the chart spec (see `geneviewer.spec`). Charts
are immutable: every method returns a new chart.

Argument names are the R names in snake_case (``itemStyle`` is
``item_style``); camelCase names are accepted too. Options passed through
``**options`` go to the spec in camelCase, as in R's ``...``.
"""

from __future__ import annotations

import functools
import json
import math
import re
import warnings
from collections.abc import Mapping
from importlib import resources
from typing import Any

import pandas as pd

from . import _r, processing
from . import data as _data
from . import spec as _spec
from ._checks import (
    break_padding_bp,
    check_align,
    check_group_style,
    check_pattern,
    check_spacing,
    check_style_by,
    check_toolbar,
    check_zoom,
    get_relative_height,
    is_str,
    quoted,
)
from .spec import _MISSING
from .theme import check_theme, resolve_theme


def _warn(message: str) -> None:
    warnings.warn(message, UserWarning, stacklevel=3)


def _any(value: Any) -> bool:
    # R's any(): True when a logical (or one of a list of them) is True
    return any(bool(v) for v in _r.as_vector(value))


def _is_flag(value: Any) -> bool:
    return _r.is_bool(value)


@functools.cache
def _validator() -> Any:
    import jsonschema

    path = resources.files("geneviewer") / "static" / "spec-v1.schema.json"
    schema = json.loads(path.read_text(encoding="utf-8"))
    return jsonschema.Draft202012Validator(schema)


def validate_spec(spec: Mapping[str, Any]) -> None:
    """Check a spec against the JSON Schema of spec v1.

    Parameters
    ----------
    spec : Mapping
        The spec, as JSON values (from ``GC.to_spec()``).

    Raises
    ------
    ValueError
        When the spec is not valid, with the first error and where it is.
    """
    errors = sorted(_validator().iter_errors(spec), key=lambda e: list(e.path))
    if errors:
        error = errors[0]
        where = "/".join(str(part) for part in error.path) or "(root)"
        raise ValueError(f"The chart spec is not valid at {where}: {error.message}")


# The defaults of each method's arguments, as in R (None for NULL)
_CHART = {
    "start": "start",
    "end": "end",
    "cluster": None,
    "group": None,
    "strand": None,
    "width": "100%",
    "height": "400px",
    "style": {},
    "save_button": True,
}


class GC:
    """A gene cluster chart.

    Build a chart from a data frame of genes, one row per gene, then add
    components with the methods. Each method returns a new chart.
    ``to_spec()`` gives the chart spec that the geneviewer JavaScript
    library draws.

    Parameters
    ----------
    data : pandas.DataFrame
        The genes: one row per gene, with its start and end positions.
    start : str
        Column of start positions. Default "start".
    end : str
        Column of end positions. Default "end".
    cluster : str, optional
        Column of cluster names: one chart row per cluster.
    group : str, optional
        Column of groups, which set the gene colours and the legend.
    strand : str, optional
        Column of strands: 1, "forward", "sense" or "+" for the forward
        strand, and -1, 0, "reverse", "antisense", "complement" or "-" for
        the reverse strand. Default: from the order of start and end.
    width : str or float
        Width of the chart, e.g. "100%" or 500. Default "100%".
    height : str or float
        Height of the chart, e.g. "400px" or 300. Default "400px".
    style : dict
        CSS styles of the chart container, e.g.
        ``{"backgroundColor": "white"}``.
    save_button : bool
        Whether the toolbar has the button to save the chart as SVG or PNG.
        Default True.

    Examples
    --------
    >>> import pandas as pd
    >>> genes = pd.DataFrame({
    ...     "start": [10, 90, 130], "end": [40, 120, 160],
    ...     "name": ["a", "b", "c"], "group": ["A", "B", "A"],
    ... })
    >>> chart = GC(genes, group="group", height="200px").labels("name")
    >>> chart.to_spec()["labels"]
    {'label': 'name'}
    """

    def __init__(
        self,
        data: pd.DataFrame,
        start: str = _MISSING,
        end: str = _MISSING,
        cluster: str | None = _MISSING,
        group: str | None = _MISSING,
        strand: str | None = _MISSING,
        width: str | float = _MISSING,
        height: str | float = _MISSING,
        style: dict[str, Any] = _MISSING,
        save_button: bool = _MISSING,
    ) -> None:
        values = dict(locals())
        del values["data"]
        a = _spec.args(values, _CHART)
        if isinstance(data, str) or hasattr(data, "__fspath__"):
            raise NotImplementedError(
                "GC() from GenBank files is not available yet in this version; "
                "pass a pandas DataFrame"
            )
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data must be a pandas DataFrame")
        data = data.copy()
        # The spec holds the data as given; the library prepares it as
        # prepare_data() does, which checks it here
        _data.prepare_data(
            data, a["start"], a["end"], a["cluster"], a["group"], a["strand"]
        )
        self._spec: dict[str, Any] = _spec.start(data, a.spec)
        # The size of the chart's output, as R's widget has it
        self._width = a["width"]
        self._height = a["height"]

    def _with(self, spec: dict[str, Any], **size: Any) -> GC:
        chart = object.__new__(GC)
        chart._spec = spec
        chart._width = size.get("width", self._width)
        chart._height = size.get("height", self._height)
        return chart

    def __repr__(self) -> str:
        parts = [k for k in self._spec if k not in ("version", "data", "fields")]
        return (
            f"<GC chart: {len(self._spec['data'])} rows"
            + (f"; {', '.join(parts)}" if parts else "")
            + ">"
        )

    def to_spec(self, validate: bool = True) -> dict[str, Any]:
        """Return the chart spec, as JSON values.

        The sequence columns ``translation`` and ``sequence`` are left out
        of the data unless an option names them (e.g. a tooltip
        ``"{translation}"``).

        Parameters
        ----------
        validate : bool
            Whether to check the spec against the JSON Schema of spec v1.
            Default True.

        Returns
        -------
        dict
            The spec: dicts, lists, strings, numbers, booleans and None,
            ready for ``json.dumps()``. Data frames are dicts of columns.
        """
        out = _data.jsonable(_spec.drop_sequence_columns(self._spec))
        if validate:
            validate_spec(out)
        return out

    def _fields(self) -> dict[str, Any]:
        return _data.chart_fields(self._spec)

    def _state(self) -> dict[str, Any]:
        return _data.chart_state(self._spec)

    # chart_genes.R ---------------------------------------------------------

    def labels(
        self,
        label: str | bool | None = _MISSING,
        show: bool = _MISSING,
        cluster: Any = _MISSING,
        item_style: list[dict[str, Any]] = _MISSING,
        group_style: dict[str, dict[str, Any]] | None = _MISSING,
        style_by: str | None = _MISSING,
        **options: Any,
    ) -> GC:
        """Add labels to the genes.

        Parameters
        ----------
        label : str, bool or None
            Column of the labels. Default: the ``group`` of the chart.
            ``False`` hides the labels.
        show : bool
            Whether to show the labels. Default True.
        cluster : int, str or list, optional
            The clusters to label, by number (from 1) or name. Default: all.
        item_style : list of dict
            Styles of single labels, e.g. ``[{"index": 0, "fill": "red"}]``.
        group_style : dict of dict, optional
            Styles by group, applied to the labels of all genes of that
            group (the ``group`` of the chart, or the column ``style_by``),
            e.g. ``{"tnsA": {"fill": "white"}}``. ``item_style`` takes
            precedence.
        style_by : str, optional
            The column ``group_style`` is matched against, or "geneTrack"
            for the track of a gene with ``cluster(prevent_gene_overlap=)``.
            Default: the ``group`` of the chart.
        **options
            Further label options, e.g. ``font_size="10px"``.
            ``adjust_labels``: True (default) turns labels that overlap
            another one, "hide" hides them, False leaves them as they are,
            and "spread" puts the labels on one line below the lowest gene,
            with a line to each gene. ``leader_line`` sets the distance and
            the style of those lines, e.g. ``{"length": 30, "stroke":
            "grey"}``, or False for no lines. ``spread_position="top"`` puts
            the spread labels above the highest gene.

        Returns
        -------
        GC
            The chart with the labels.
        """
        a = _spec.args(locals(), _LABELS)
        state = self._fields()
        check_group_style(a["group_style"])
        check_style_by(state, a["style_by"])
        adjust = a.options.get("adjustLabels")
        if adjust is not None and not (
            (_is_flag(adjust)) or adjust in ("hide", "spread")
        ):
            raise ValueError(
                '\'adjust_labels\' must be True, False, "hide" or "spread".'
            )
        position = a.options.get("spreadPosition")
        if position is not None and position not in ("bottom", "top"):
            raise ValueError('\'spread_position\' must be "bottom" or "top".')
        leader = a.options.get("leaderLine")
        if (
            leader is not None
            and leader is not False
            and not isinstance(leader, Mapping)
        ):
            raise ValueError(
                '\'leader_line\' must be False or a dict, e.g. {"stroke": "grey"}.'
            )
        if isinstance(leader, Mapping) and leader.get("length") is not None:
            length = leader["length"]
            if not (_r.is_number(length) and math.isfinite(length)):
                raise ValueError(
                    "'leader_line[\"length\"]' must be a single number of pixels, "
                    "e.g. 30."
                )

        if not _any(a["show"]):
            return self
        label = a["label"]
        if _is_flag(label):
            label = None
        if label is None and state["group"] is None:
            raise ValueError("Please define labels")
        if label is None:
            label = state["group"]
        columns = state["columns"] or []
        if not all(v in columns for v in _r.as_vector(label)):
            raise ValueError("label column not found in data")
        _data.get_updated_clusters(state["clusters"], a["cluster"])
        return self._with(_spec.add(self._spec, "labels", _spec.flag(a.spec, "label")))

    def genes(
        self,
        group: str | bool | None = _MISSING,
        marker: str | None = _MISSING,
        marker_size: str | None = _MISSING,
        show: bool = _MISSING,
        color_scheme: str | None = _MISSING,
        custom_colors: Any = _MISSING,
        cluster: Any = _MISSING,
        item_style: list[dict[str, Any]] = _MISSING,
        group_style: dict[str, dict[str, Any]] | None = _MISSING,
        style_by: str | None = _MISSING,
        corner_radius: float | None = _MISSING,
        **options: Any,
    ) -> GC:
        """Set how the genes are drawn.

        Parameters
        ----------
        group : str, optional
            Column of groups, for the gene colours. Default: the ``group``
            of the chart.
        marker : str, optional
            The shape of the genes: "arrow", "boxarrow", "box", "cbox" or
            "rbox".
        marker_size : str, optional
            "small", "medium" or "large".
        show : bool
            Whether to show the genes. Default True.
        color_scheme : str, optional
            A d3 colour scheme (e.g. "schemeSet2") or a geneviewer palette
            (e.g. "okabeIto", see `geneviewer_palette`). Without
            ``cluster``, the legend uses the same colours.
        custom_colors : list or dict, optional
            Colours recycled over the groups, or a dict of colours by group;
            groups left out of a dict keep the colours of the scheme.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        item_style : list of dict
            Styles of single genes, e.g. ``[{"index": 2, "fill": "red"}]``.
        group_style : dict of dict, optional
            Styles by group, applied to every gene of that group, e.g.
            ``{"promoter": {"marker": "box", "fill": "grey"}}``.
            ``item_style`` takes precedence.
        style_by : str, optional
            The column ``group_style`` is matched against, or "geneTrack".
            Default: ``group``.
        corner_radius : float, optional
            Radius in pixels of the rounded corners of "arrow", "boxarrow",
            "box" and "rbox" markers. Default: sharp corners, and 5 for
            "rbox".
        **options
            Further gene options, e.g. ``stroke="black"``,
            ``stroke_width=1``, ``marker_height=12``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _GENES)
        state = self._fields()
        check_group_style(a["group_style"])
        check_style_by(state, a["style_by"])
        group = a["group"]
        if _is_flag(group):
            group = None
        if group is None and state["group"] is not None:
            group = state["group"]
        if group is not None and not all(
            v in (state["columns"] or []) for v in _r.as_vector(group)
        ):
            raise ValueError("group column not found in data")
        _data.get_updated_clusters(state["clusters"], a["cluster"])
        return self._with(_spec.add(self._spec, "genes", _spec.flag(a.spec, "group")))

    def transcript(
        self,
        transcript: str | None = _MISSING,
        type: str | None = _MISSING,
        strand: str | None = _MISSING,
        group: str | bool | None = _MISSING,
        selection: Any = _MISSING,
        show: bool = _MISSING,
        color_scheme: str | None = _MISSING,
        custom_colors: Any = _MISSING,
        style_exons: dict[str, Any] = _MISSING,
        style_introns: dict[str, Any] = _MISSING,
        style_utrs: dict[str, Any] = _MISSING,
        item_style_exons: list[dict[str, Any]] = _MISSING,
        item_style_introns: list[dict[str, Any]] = _MISSING,
        item_style_utrs: list[dict[str, Any]] = _MISSING,
        group_style_exons: dict[str, dict[str, Any]] | None = _MISSING,
        group_style_introns: dict[str, dict[str, Any]] | None = _MISSING,
        group_style_utrs: dict[str, dict[str, Any]] | None = _MISSING,
        style_by: str | None = _MISSING,
        label_options: dict[str, Any] = _MISSING,
        **options: Any,
    ) -> GC:
        """Draw the data as transcripts: exons, introns and UTRs.

        The first call makes the chart again with the transcripts as
        clusters, and adds an intron between each pair of exons to the
        data. Settings of the whole chart made before it (title, legend,
        theme, grid margins) are kept; settings per cluster must come after
        it.

        Parameters
        ----------
        transcript : str, optional
            Column of the transcript of each row. Default "transcript".
        type : str, optional
            Column of the feature type: "exon", or a UTR ("5_utr", "3_utr").
        strand : str, optional
            Column of strands, as for ``GC()``.
        group : str, optional
            Column of groups, for the colours. Default: the ``group`` of
            the chart.
        selection : int, str or list, optional
            The transcripts, by number (from 1) or name. Default: all.
        show : bool
            Whether to show the transcripts. Default True.
        color_scheme, custom_colors
            As for ``genes()``.
        style_exons, style_introns, style_utrs : dict
            Styles of all exons, introns or UTRs.
        item_style_exons, item_style_introns, item_style_utrs : list of dict
            Styles of single exons, introns or UTRs.
        group_style_exons, group_style_introns, group_style_utrs : dict of dict
            Styles by group, e.g. ``{"tx1": {"fill": "red"}}``.
        style_by : str, optional
            The column the group styles are matched against. Default:
            ``group``.
        label_options : dict
            Options of the transcript labels, e.g. ``{"fontSize": "12px"}``.
        **options
            Further options.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _TRANSCRIPT)
        state = self._fields()
        check_group_style(a["group_style_exons"], "group_style_exons")
        check_group_style(a["group_style_introns"], "group_style_introns")
        check_group_style(a["group_style_utrs"], "group_style_utrs")
        check_style_by(state, a["style_by"])
        new_spec = self._spec
        # The first call adds the introns to the data and draws the
        # transcripts as the clusters
        if self._spec.get("transcript") is None:
            data = self._state()["data"]
            data_tmp = data
            columns = list(data.columns)
            type_, transcript = a["type"], a["transcript"]
            if type_ is not None and not is_str(type_):
                raise ValueError("type must be a single column name")
            if type_ is not None and type_ not in columns:
                raise ValueError("type column not found in data")
            data = data.copy()
            if type_ is not None:
                data["type"] = data_tmp[type_]
            if transcript is not None and not is_str(transcript):
                raise ValueError("transcript must be a single column name")
            if transcript is not None and transcript not in columns:
                raise ValueError("transcript column not found in data")
            if transcript is not None:
                data["transcript"] = [_r.as_character(v) for v in data_tmp[transcript]]
            data = processing.get_introns(data)
            # The data of the chart with the transcripts as clusters, checked
            # as GC() checks it
            cluster = "transcript" if transcript is None else transcript
            prepared = _data.prepare_data(
                data, state["start"], state["end"], cluster, state["group"], a["strand"]
            )
            state["columns"] = list(prepared.columns)
            state["clusters"] = _data.state_clusters(prepared, cluster)
            new_spec = {**self._spec, "data": data}

        group = a["group"]
        if _is_flag(group):
            group = None
        if group is None and state["group"] is not None:
            group = state["group"]
        if group is not None and not all(
            v in (state["columns"] or []) for v in _r.as_vector(group)
        ):
            raise ValueError("group column not found in data")
        _data.get_updated_clusters(state["clusters"], a["selection"])
        return self._with(
            _spec.add(
                new_spec,
                "transcript",
                _spec.flag(a.spec, "group"),
                rename={
                    "style_utrs": "styleUTRs",
                    "item_style_utrs": "itemStyleUTRs",
                    "group_style_utrs": "groupStyleUTRs",
                },
            )
        )

    def color(
        self, color_scheme: str | None = _MISSING, custom_colors: Any = _MISSING
    ) -> GC:
        """Set the colours of the genes and the legend.

        Parameters
        ----------
        color_scheme : str, optional
            A d3 colour scheme (e.g. "schemeCategory10") or a geneviewer
            palette (e.g. "okabeIto", see `geneviewer_palette`).
        custom_colors : list or dict, optional
            Colours recycled over the groups, or a dict of colours by group.
            Groups left out of a dict keep the colours of ``color_scheme``,
            the palette of the theme or the default colours.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _COLOR)
        # color_scheme only colours the groups left out of a dict of colours
        if (
            a["color_scheme"] is not None
            and a["custom_colors"] is not None
            and not isinstance(a["custom_colors"], Mapping)
        ):
            _warn(
                "Both color_scheme and custom_colors are selected. Use a dict of "
                "custom_colors to set the colours of some groups, or choose only one."
            )
            return self
        return self._with(_spec.add(self._spec, "color", a.spec, mode="replace"))

    def pattern(
        self,
        patterns: Any = _MISSING,
        color: Any = _MISSING,
        background: Any = _MISSING,
        size: float = _MISSING,
        stroke_width: float = _MISSING,
    ) -> GC:
        """Fill the genes, exons and legend markers with a pattern per group.

        For black and white figures, or to show each group with both a
        colour and a pattern.

        Parameters
        ----------
        patterns : list or dict, optional
            Pattern names recycled over the groups, or a dict of pattern
            names by group (groups left out get "none"). Default: a pattern
            per group in a fixed order. Patterns: "none", "solid",
            "horizontal", "vertical", "diagonal", "antidiagonal", "grid",
            "crosshatch", "dots", "zigzag", "checkerboard", "rings",
            "bricks", "plus", "wideStripes" and "inverseDots".
        color : str, list or dict
            Colour of the pattern lines and dots: one, recycled over the
            groups, or by group. Default "black".
        background : str, list, dict or None
            Colour behind the pattern, as ``color``. Default "white"; None
            keeps the gene colours.
        size : float
            Width and height of a pattern tile in pixels. Default 6.
        stroke_width : float
            Width of the pattern lines in pixels. Default 1.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _PATTERN)
        check_pattern({"patterns": a["patterns"]})
        return self._with(_spec.add(self._spec, "pattern", a.spec, mode="replace"))

    def legend(
        self,
        group: str | bool | None = _MISSING,
        show: bool = _MISSING,
        background_color: str | None = _MISSING,
        order: list[Any] = _MISSING,
        position: str | None = _MISSING,
        align: str | None = _MISSING,
        style: dict[str, Any] = _MISSING,
        legend_options: dict[str, Any] = _MISSING,
        legend_text_options: dict[str, Any] = _MISSING,
        corner_radius: float | None = _MISSING,
        group_style: dict[str, dict[str, Any]] | None = _MISSING,
        **options: Any,
    ) -> GC:
        """Set the legend.

        Parameters
        ----------
        group : str or bool, optional
            Column of the groups in the legend. Default: the ``group`` of
            the chart. False hides the legend.
        show : bool
            Whether to show the legend. Default True.
        background_color : str, optional
            Background colour of the legend; ``style={"backgroundColor":
            ...}`` takes precedence. Default: the theme's, else transparent.
        order : list
            The order of the groups in the legend.
        position : str, optional
            "top" or "bottom". Default: the theme's, else "bottom".
        align : str, optional
            "left", "center" or "right". Default: the theme's, else "left".
        style : dict
            CSS styles of the legend.
        legend_options : dict
            Options of the markers (``markerStyle`` in the spec).
        legend_text_options : dict
            Options of the text (``textStyle`` in the spec).
        corner_radius : float, optional
            Radius in pixels of the rounded corners of the markers.
        group_style : dict of dict, optional
            Styles of the entry of a group, e.g. ``{"tnsA": {"fontStyle":
            "italic"}}``; styles in ``marker`` go to its marker.
        **options
            Further options, e.g. ``orientation="vertical"``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _LEGEND)
        state = self._fields()
        check_group_style(a["group_style"])
        group, show = a["group"], a["show"]
        if _is_flag(group):
            show, group = group, None
        if show is True and group is None and state["group"] is None:
            raise ValueError("Please define a group")
        if group is None and state["group"] is not None:
            group = state["group"]
        if group is not None and group not in (state["columns"] or []):
            raise ValueError("group column not found in data")
        check_align(a["align"])
        return self._with(
            _spec.add(
                self._spec,
                "legend",
                _spec.flag(a.spec, "group"),
                mode="merge",
                rename={
                    "legend_options": "markerStyle",
                    "legend_text_options": "textStyle",
                },
            )
        )

    # chart_layout.R --------------------------------------------------------

    def title(
        self,
        title: str | None = _MISSING,
        subtitle: str | None = _MISSING,
        style: dict[str, Any] = _MISSING,
        subtitle_font: dict[str, Any] = _MISSING,
        title_font: dict[str, Any] = _MISSING,
        show: bool = _MISSING,
        height: Any = _MISSING,
        align: str | None = _MISSING,
        spacing: float | None = _MISSING,
        **options: Any,
    ) -> GC:
        """Add a title and subtitle to the chart.

        Parameters
        ----------
        title, subtitle : str, optional
            The texts; HTML is allowed.
        style : dict
            CSS styles of the title container.
        subtitle_font, title_font : dict
            Font options, e.g. ``{"fontSize": "18px"}``.
        show : bool
            Whether to show the title. Default True.
        height : str or float, optional
            Height of the title, e.g. "40px". Default: from the text.
        align : str, optional
            "left", "center" or "right". Default "center".
        spacing : float, optional
            Pixels between the title and the subtitle.
        **options
            Further options, e.g. ``position``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _TITLE)
        if not _is_flag(a["show"]):
            raise ValueError("show must be a single True or False")
        check_align(a["align"])
        check_spacing(a["spacing"])
        if not a["show"] or (a["title"] is None and a["subtitle"] is None):
            return self
        return self._with(_spec.add(self._spec, "title", a.spec, mode="merge"))

    def cluster_title(
        self,
        title: Any = _MISSING,
        subtitle: Any = _MISSING,
        subtitle_font: dict[str, Any] = _MISSING,
        title_font: dict[str, Any] = _MISSING,
        show: bool = _MISSING,
        height: Any = _MISSING,
        cluster: Any = _MISSING,
        align: Any = _MISSING,
        spacing: Any = _MISSING,
        **options: Any,
    ) -> GC:
        """Add a title and subtitle above each cluster.

        Parameters
        ----------
        title, subtitle : str or list of str, optional
            The texts, one per cluster or recycled.
        subtitle_font, title_font : dict
            Font options.
        show : bool
            Whether to show the titles. Default True.
        height : str or float, optional
            Height of the titles.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        align : str or list of str, optional
            "left", "center" or "right", per cluster or for all.
        spacing : float or list of float, optional
            Pixels between title and subtitle, per cluster or for all.
        **options
            Further options, e.g. ``position``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _CLUSTER_TITLE)
        if not _any(a["show"]):
            return self
        check_align(a["align"], per_cluster=True)
        check_spacing(a["spacing"], per_cluster=True)
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(_spec.add(self._spec, "clusterTitle", a.spec))

    def sequence(
        self,
        show: bool = _MISSING,
        cluster: Any = _MISSING,
        y: float | None = _MISSING,
        sequence_style: dict[str, Any] = _MISSING,
        marker_style: dict[str, Any] = _MISSING,
        **options: Any,
    ) -> GC:
        """Set the line of each cluster's sequence.

        Parameters
        ----------
        show : bool
            Whether to show the line. Default True.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        y : float, optional
            Vertical position of the line.
        sequence_style : dict
            Style of the line (``lineStyle`` in the spec), e.g.
            ``{"stroke": "grey"}``.
        marker_style : dict
            Style of the markers of scale breaks.
        **options
            Further options.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _SEQUENCE)
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(
            _spec.add(
                self._spec, "sequence", a.spec, rename={"sequence_style": "lineStyle"}
            )
        )

    def grid(
        self,
        margin: dict[str, Any] | None = _MISSING,
        width: Any = _MISSING,
        height: Any = _MISSING,
        direction: str = _MISSING,
        cluster: Any = _MISSING,
    ) -> GC:
        """Set the layout of the clusters: margins, width and height.

        Parameters
        ----------
        margin : dict, optional
            Margins, e.g. ``{"left": 50, "top": "10%"}``. Top and bottom
            margins of the whole chart are padding of the chart.
        width : str or float, optional
            Width of the chart (without ``cluster``) or of the clusters.
        height : str or float, optional
            Height of the clusters, e.g. "30%" or 100.
        direction : str
            "column" or "row". Default "column".
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: the whole
            chart.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _GRID)
        _data.check_clusters(self._spec, a["cluster"])
        size: dict[str, Any] = {}
        margin, width = a["margin"], a["width"]
        # The width of the chart's output
        if width is not None and a["cluster"] is None:
            size["width"] = get_relative_height(self._width, width)
        # Top and bottom margins are padding of the whole chart
        if margin is not None and a["cluster"] is None:
            top = margin.get("top") if isinstance(margin, Mapping) else None
            bottom = margin.get("bottom") if isinstance(margin, Mapping) else None
            padding_top = 0 if top is None else get_relative_height(self._height, top)
            padding_bottom = (
                0 if bottom is None else get_relative_height(self._height, bottom)
            )
            if not (_r.is_number(padding_top) and _r.is_number(padding_bottom)):
                raise ValueError(
                    "Use pixel values for the top and bottom margins when the chart "
                    "height is a percentage."
                )
            # The library takes the padding from the cluster heights; the
            # output fills its container
            if padding_top > 0 or padding_bottom > 0:
                size["height"] = "100%"
        return self._with(_spec.add(self._spec, "grid", a.spec), **size)

    def scale(
        self,
        cluster: Any = _MISSING,
        start: Any = _MISSING,
        end: Any = _MISSING,
        padding: Any = _MISSING,
        hidden: bool = _MISSING,
        breaks: Any = _MISSING,
        tick_values: Any = _MISSING,
        reverse: bool = _MISSING,
        axis_position: str | None = _MISSING,
        axis_type: str = _MISSING,
        y: float | None = _MISSING,
        scale_breaks: bool = _MISSING,
        scale_break_threshold: float = _MISSING,
        scale_break_padding: float | str = _MISSING,
        ticks_count: int | None = _MISSING,
        ticks_format: str | None = _MISSING,
        tick_style: dict[str, Any] = _MISSING,
        text_style: dict[str, Any] = _MISSING,
        line_style: dict[str, Any] = _MISSING,
        **options: Any,
    ) -> GC:
        """Set the scale of the clusters, and add an axis.

        Parameters
        ----------
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        start, end : float or list of float, optional
            Start and end of the scale, per cluster or for all.
        padding : float or str, optional
            Padding around the genes, in percent or e.g. "500bp".
        hidden : bool
            Whether to hide the axis. Default False.
        breaks : list of dict or pandas.DataFrame
            Breaks in the scale, e.g. ``[{"start": 100, "end": 200}]``, or a
            data frame with ``start`` and ``end`` columns and, for breaks
            per cluster, a ``cluster`` column.
        tick_values : list of float, optional
            Positions of the ticks.
        reverse : bool
            Whether to reverse the scale. Default False.
        axis_position : str, optional
            "bottom" or "top".
        axis_type : str
            "position" or "range". Default "position".
        y : float, optional
            Vertical position of the axis.
        scale_breaks : bool
            Whether to break the scale where genes are far apart
            (``autoBreaks`` in the spec). Default False.
        scale_break_threshold : float
            Gaps longer than this percentage of the cluster get a break
            (``autoBreakThreshold``). Default 20.
        scale_break_padding : float or str
            Room around automatic breaks, in percent or e.g. "500bp"
            (``breakPadding``). Default 1.
        ticks_count : int, optional
            Number of ticks.
        ticks_format : str, optional
            d3 format of the tick labels, e.g. ",.0f".
        tick_style, text_style, line_style : dict
            Styles of the ticks, their text and the axis line.
        **options
            Further options.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _SCALE)
        # Fail early on a padding that is no percentage or base pairs
        if is_str(a["scale_break_padding"]):
            break_padding_bp(a["scale_break_padding"])
        clusters = _data.get_updated_clusters(self._fields()["clusters"], a["cluster"])
        breaks = a["breaks"]
        if (
            not a["scale_breaks"]
            and clusters
            and isinstance(breaks, pd.DataFrame)
            and not {"start", "end"} <= set(breaks.columns)
        ):
            raise ValueError("A data frame of breaks needs a start and an end column")
        return self._with(
            _spec.add(
                self._spec,
                "scale",
                a.spec,
                rename={
                    "scale_breaks": "autoBreaks",
                    "scale_break_threshold": "autoBreakThreshold",
                    "scale_break_padding": "breakPadding",
                },
            )
        )

    def scale_bar(
        self,
        show: bool = _MISSING,
        cluster: Any = _MISSING,
        scale_bar_line_style: dict[str, Any] = _MISSING,
        scale_bar_tick_style: dict[str, Any] = _MISSING,
        label_style: dict[str, Any] = _MISSING,
        **options: Any,
    ) -> GC:
        """Add a scale bar to the clusters.

        Parameters
        ----------
        show : bool
            Whether to show the scale bar. Default True.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        scale_bar_line_style, scale_bar_tick_style : dict
            Styles of the line and the ticks (``lineStyle`` and
            ``tickStyle`` in the spec).
        label_style : dict
            Style of the label.
        **options
            Further options, e.g. ``title="1 kb"``, ``scale_bar_unit=1000``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _SCALE_BAR)
        if not _any(a["show"]):
            return self
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(
            _spec.add(
                self._spec,
                "scaleBar",
                a.spec,
                rename={
                    "scale_bar_line_style": "lineStyle",
                    "scale_bar_tick_style": "tickStyle",
                },
            )
        )

    def cluster_label(
        self,
        title: Any = _MISSING,
        show: bool = _MISSING,
        width: Any = _MISSING,
        cluster: Any = _MISSING,
        position: str | None = _MISSING,
        align: Any = _MISSING,
        wrap_label: bool | None = _MISSING,
        wrap_options: dict[str, Any] = _MISSING,
        subtitle: Any = _MISSING,
        subtitle_font: dict[str, Any] = _MISSING,
        spacing: Any = _MISSING,
        **options: Any,
    ) -> GC:
        """Add a label beside each cluster.

        Parameters
        ----------
        title : str or list of str, optional
            The labels, one per cluster or recycled. Default: the cluster
            names.
        show : bool
            Whether to show the labels. Default True.
        width : str or float, optional
            Width of the label column.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        position : str, optional
            "left" or "right".
        align : str or list of str, optional
            "left", "center" or "right".
        wrap_label : bool, optional
            Whether to wrap long labels.
        wrap_options : dict
            Options of the wrapping.
        subtitle : str or list of str, optional
            A second line under each label.
        subtitle_font : dict
            Font options of the subtitle.
        spacing : float or list of float, optional
            Pixels between label and subtitle.
        **options
            Further options, e.g. ``font_weight="bold"``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _CLUSTER_LABEL)
        if not _any(a["show"]):
            return self
        check_align(a["align"], per_cluster=True)
        check_spacing(a["spacing"], per_cluster=True)
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(_spec.add(self._spec, "clusterLabel", a.spec))

    def cluster_footer(
        self,
        title: Any = _MISSING,
        subtitle: Any = _MISSING,
        show: bool = _MISSING,
        cluster: Any = _MISSING,
        subtitle_font: dict[str, Any] = _MISSING,
        title_font: dict[str, Any] = _MISSING,
        align: Any = _MISSING,
        spacing: Any = _MISSING,
        **options: Any,
    ) -> GC:
        """Add a footer below each cluster.

        Parameters
        ----------
        title, subtitle : str or list of str, optional
            The texts, one per cluster or recycled.
        show : bool
            Whether to show the footers. Default True.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        subtitle_font, title_font : dict
            Font options.
        align : str or list of str, optional
            "left", "center" or "right".
        spacing : float or list of float, optional
            Pixels between title and subtitle.
        **options
            Further options.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _CLUSTER_FOOTER)
        if not _any(a["show"]):
            return self
        check_align(a["align"], per_cluster=True)
        check_spacing(a["spacing"], per_cluster=True)
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(_spec.add(self._spec, "clusterFooter", a.spec))

    def coordinates(
        self,
        show: bool = _MISSING,
        tick_values_top: Any = _MISSING,
        tick_values_bottom: Any = _MISSING,
        ticks_format: str | None = _MISSING,
        tick_style: dict[str, Any] = _MISSING,
        text_style: dict[str, Any] = _MISSING,
        cluster: Any = _MISSING,
        group_style: dict[str, dict[str, Any]] | None = _MISSING,
        style_by: str | None = _MISSING,
        **options: Any,
    ) -> GC:
        """Show the start and end coordinates of the genes.

        Parameters
        ----------
        show : bool
            Whether to show the coordinates. Default True.
        tick_values_top, tick_values_bottom : list of float, optional
            Coordinates to show above or below the genes. Default: those of
            the genes.
        ticks_format : str, optional
            d3 format of the coordinates.
        tick_style, text_style : dict
            Styles of the ticks and their text.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        group_style : dict of dict, optional
            Styles of the coordinates of the genes of a group: ``show`` and
            ``tick`` and text styles, e.g. ``{"other": {"show": False}}``.
        style_by : str, optional
            The column ``group_style`` is matched against. Default: the
            ``group`` of the chart.
        **options
            Further options, e.g. ``rotate=-45``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _COORDINATES)
        state = self._fields()
        check_group_style(a["group_style"])
        check_style_by(state, a["style_by"])
        style_by = a["style_by"] if a["style_by"] is not None else state["group"]
        if a["group_style"] is not None and style_by is None:
            raise ValueError("group_style needs a group in GC() or a style_by column")
        if not _any(a["show"]):
            return self
        _data.get_updated_clusters(state["clusters"], a["cluster"])
        return self._with(_spec.add(self._spec, "coordinates", a.spec))

    def cluster(
        self,
        separate_strands: bool | None = _MISSING,
        strand_spacing: float | None = _MISSING,
        prevent_gene_overlap: bool | str | None = _MISSING,
        overlap_min: float | str | None = _MISSING,
        overlap_spacing: float | None = _MISSING,
        cluster: Any = _MISSING,
        style: dict[str, Any] = _MISSING,
        **options: Any,
    ) -> GC:
        """Set how genes are placed in their cluster.

        Parameters
        ----------
        separate_strands : bool, optional
            Whether to draw the strands on separate lines.
        strand_spacing : float, optional
            Pixels between the strands.
        prevent_gene_overlap : bool or str, optional
            True moves overlapping genes to other tracks; "nested" draws
            genes inside longer ones smaller, on the same track.
        overlap_min : float or str, optional
            Overlap genes may have before they move, in base pairs or as a
            share, e.g. "10%".
        overlap_spacing : float, optional
            Pixels between tracks.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        style : dict
            CSS styles of the cluster.
        **options
            Further options.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _CLUSTER)
        overlap = a["prevent_gene_overlap"]
        if overlap is not None and not (_is_flag(overlap) or overlap == "nested"):
            raise ValueError(
                "'prevent_gene_overlap' must be True, False or \"nested\"."
            )
        overlap_min = a["overlap_min"]
        if overlap_min is not None:
            valid = (_r.is_number(overlap_min) and overlap_min >= 0) or (
                is_str(overlap_min)
                and re.fullmatch(r"[0-9]+(\.[0-9]+)?%", overlap_min) is not None
                and float(overlap_min[:-1]) <= 100
            )
            if not valid:
                raise ValueError(
                    "'overlap_min' must be a number of base pairs or a share such "
                    'as "10%".'
                )
        if (
            a["separate_strands"] is True
            and overlap is not None
            and overlap is not False
        ):
            _warn(
                "Setting both 'separate_strands' and 'prevent_gene_overlap' to True "
                "is not supported. Resetting to None."
            )
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(_spec.add(self._spec, "cluster", a.spec))

    # chart_interaction.R ---------------------------------------------------

    def annotation(
        self, type: str = _MISSING, cluster: Any = _MISSING, **options: Any
    ) -> GC:
        """Add an annotation: text, a line, an arrow, a symbol and more.

        Parameters
        ----------
        type : str
            "text", "textMarker", "line", "arrow", "symbol", "rectangle",
            "promoter", "terminator" or "arc". Default "text".
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        **options
            The options of the annotation type, e.g. ``text="Here", x=100,
            y=70``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _ANNOTATION)
        if not is_str(a["type"]) or a["type"] not in ANNOTATION_TYPES:
            raise ValueError(f"'type' must be one of {quoted(ANNOTATION_TYPES)}.")
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(_spec.add(self._spec, "annotation", a.spec, mode="items"))

    def track_mouse(
        self, show: bool = _MISSING, cluster: Any = _MISSING, **options: Any
    ) -> GC:
        """Show a line and the position under the mouse.

        Parameters
        ----------
        show : bool
            Whether to track the mouse. Default True.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        **options
            Further options, e.g. ``cursor="cell"``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _TRACK_MOUSE)
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(_spec.add(self._spec, "trackMouse", a.spec))

    def events(
        self,
        click: bool = _MISSING,
        highlight: bool = _MISSING,
        rendered: bool = _MISSING,
        render_id: str | float | None = _MISSING,
        dblclick: bool = _MISSING,
        select: bool = _MISSING,
        legend_click: bool = _MISSING,
        drag: bool = _MISSING,
        zoom: bool = _MISSING,
    ) -> GC:
        """Make the chart report clicks, selections, drags and zooms to Shiny.

        The payload of each event is described at
        https://nvelden.github.io/geneviewer-docs/api#events.

        Parameters
        ----------
        click : bool
            Report clicks as the input ``<id>_click``. Default True.
        highlight : bool
            Outline the part under the mouse. Default False.
        rendered : bool
            Report each drawn update as ``<id>_rendered``. Default False.
        render_id : str or float, optional
            Identifier of the update; a new one marks a new load.
        dblclick : bool
            Report double-clicks as ``<id>_dblclick``. Default False.
        select : bool
            Select genes and labels: ``<id>_select`` and
            ``<id>_selection``. Default False.
        legend_click : bool
            Whether clicking a legend entry hides its genes. Default True.
        drag : bool
            Move labels, titles, scales and legend parts: ``<id>_drag`` and
            ``<id>_moves``. Default False.
        zoom : bool
            Report the window of the zoom slider as ``<id>_zoom``. Default
            False.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _EVENTS)
        for name in (
            "click",
            "highlight",
            "rendered",
            "dblclick",
            "select",
            "legend_click",
            "drag",
            "zoom",
        ):
            if not _is_flag(a[name]):
                raise ValueError("Event options must be True or False.")
        render_id = a["render_id"]
        if render_id is not None and not (
            is_str(render_id) or (_r.is_number(render_id) and math.isfinite(render_id))
        ):
            raise ValueError(
                "'render_id' must be None, a single string or a finite number."
            )
        return self._with(_spec.add(self._spec, "events", a.spec, mode="replace"))

    def toolbar(
        self,
        show: str | None = _MISSING,
        position: str | None = _MISSING,
        buttons: Any = _MISSING,
        size: float | None = _MISSING,
        color: str | None = _MISSING,
        hover_color: str | None = _MISSING,
        tooltip_background: str | None = _MISSING,
        tooltip_color: str | None = _MISSING,
        tooltip_border: str | None = _MISSING,
        tooltip_radius: float | str | None = _MISSING,
        tooltip_shadow: str | None = _MISSING,
        font_family: str | None = _MISSING,
        icons: dict[str, str] | None = _MISSING,
    ) -> GC:
        """Set the toolbar: when it shows, its buttons, icons and colours.

        Custom icons are Font Awesome icons by name, which need Font Awesome
        on the page.

        Parameters
        ----------
        show : str, optional
            "hover", "always" or "never".
        position : str, optional
            "topright", "topleft", "bottomright" or "bottomleft".
        buttons : list of str, optional
            The buttons, of "undo", "redo", "reset", "zoom" and "save".
        size : float, optional
            Size of the icons in pixels.
        color, hover_color : str, optional
            Colours of the icons.
        tooltip_background, tooltip_color, tooltip_border : str, optional
            Style of the tooltips and the save menu.
        tooltip_radius : float or str, optional
            Corner radius of the tooltips.
        tooltip_shadow : str, optional
            CSS box shadow of the tooltips.
        font_family : str, optional
            Font of the tooltips.
        icons : dict, optional
            Font Awesome icons by button, e.g. ``{"reset":
            "arrows-to-circle"}``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _TOOLBAR)
        options = {
            _spec.snake_to_camel(name): a[name]
            for name in _TOOLBAR
            if a[name] is not None
        }
        check_toolbar(options)
        call_args = a.spec.copy()
        if call_args.get("buttons") is not None:
            call_args["buttons"] = _r.as_vector(call_args["buttons"])
        return self._with(_spec.add(self._spec, "toolbar", call_args, mode="merge"))

    def zoom(
        self,
        show: bool = _MISSING,
        range: Any = _MISSING,
        mode: str | None = _MISSING,
        min: float | None = _MISSING,
        update: str | None = _MISSING,
        height: float | None = _MISSING,
        overview: bool | None = _MISSING,
        color: str | None = _MISSING,
        style: dict[str, Any] | None = _MISSING,
        window_style: dict[str, Any] | None = _MISSING,
        handle_style: dict[str, Any] | None = _MISSING,
    ) -> GC:
        """Add a zoom slider below the clusters.

        Parameters
        ----------
        show : bool
            Whether to show the slider. Default True.
        range : list of float, optional
            The window shown first: shares of the clusters, e.g. ``[0.25,
            0.5]``, or base pairs with ``mode="absolute"``.
        mode : str, optional
            "relative" (default): the same share of every cluster;
            "absolute": the same base pairs.
        min : float, optional
            The narrowest window.
        update : str, optional
            "live" or "release": draw while dragging, or when let go.
        height : float, optional
            Height of the slider in pixels.
        overview : bool, optional
            Whether the slider shows the genes of each cluster.
        color : str, optional
            Colour of the slider.
        style, window_style, handle_style : dict, optional
            CSS styles of the slider, its window and its handles.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _ZOOM)
        check_zoom({name: a[name] for name in _ZOOM if a[name] is not None})
        return self._with(_spec.add_zoom(self._spec, a.spec))

    def tooltip(
        self,
        formatter: str | bool = _MISSING,
        show: bool = _MISSING,
        cluster: Any = _MISSING,
        **options: Any,
    ) -> GC:
        """Set the tooltip of the genes.

        Parameters
        ----------
        formatter : str or bool
            The text, with column values in braces. Default "<b>Start:</b>
            {start}<br><b>End:</b> {end}". False hides the tooltip.
        show : bool
            Whether to show the tooltip. Default True.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        **options
            Further options, e.g. ``background_color="black"``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _TOOLTIP)
        _data.check_clusters(self._spec, a["cluster"])
        return self._with(
            _spec.add(self._spec, "tooltip", _spec.flag(a.spec, "formatter"))
        )

    # chart_links.R ---------------------------------------------------------

    def links(
        self,
        group: str | None = _MISSING,
        data: pd.DataFrame | None = _MISSING,
        value1: Any = _MISSING,
        value2: Any = _MISSING,
        cluster: Any = _MISSING,
        curve: bool | None = _MISSING,
        measure: str = _MISSING,
        show_links: bool = _MISSING,
        label: bool = _MISSING,
        label_position: str = _MISSING,
        normal_color: str | None = _MISSING,
        inverted_color: str | None = _MISSING,
        low_color: str | None = _MISSING,
        vmin: float | None = _MISSING,
        vmax: float | None = _MISSING,
        use_group_colors: bool = _MISSING,
        color_bar: bool | None = _MISSING,
        color_bar_options: dict[str, Any] = _MISSING,
        link_width: float | None = _MISSING,
        link_style: dict[str, Any] = _MISSING,
        label_style: dict[str, Any] = _MISSING,
        **options: Any,
    ) -> GC:
        """Add links between genes of neighbouring clusters.

        Parameters
        ----------
        group : str, optional
            Column matched between clusters: genes with the same value are
            linked.
        data : pandas.DataFrame, optional
            Links to draw instead, e.g. from ``protein_blast()``.
        value1, value2 : list, optional
            Pairs of group values to link.
        cluster : list, optional
            The clusters to link, at least two. Default: all.
        curve : bool, optional
            Whether links are curved.
        measure : str
            "identity", "similarity" or "none": what colours and labels the
            links. Default "identity".
        show_links : bool
            Whether to show the links (``show`` in the spec). Default True.
        label : bool
            Whether to label the links with the measure. Default True.
        label_position : str
            "gene" or "link". Default "gene".
        normal_color, inverted_color, low_color : str, optional
            Colours of links, of inverted links, and of low values.
        vmin, vmax : float, optional
            The range of the measure for the colours.
        use_group_colors : bool
            Whether links take the colour of their group. Default False.
        color_bar : bool, optional
            Whether to show a colour bar.
        color_bar_options : dict
            Options of the colour bar: ``position``, ``orientation``,
            ``labels``, ``legendSide``.
        link_width : float, optional
            Width of the links.
        link_style, label_style : dict
            Styles of the links and their labels.
        **options
            Further options.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _LINKS)
        measure = a["measure"]
        if not is_str(measure) or measure not in ("identity", "similarity", "none"):
            raise ValueError(
                "Invalid measure specified. Choose 'identity', 'similarity', or 'none'."
            )
        position = a["label_position"]
        if not is_str(position) or position not in ("gene", "link"):
            raise ValueError(
                f"'label_position' must be one of {quoted(('gene', 'link'))}."
            )
        bar_options = a["color_bar_options"] or {}
        for name, choices in _COLOR_BAR_CHOICES.items():
            value = bar_options.get(name) if isinstance(bar_options, Mapping) else None
            if value is not None and not (is_str(value) and value in choices):
                raise ValueError(
                    f"'color_bar_options[\"{name}\"]' must be one of {quoted(choices)}."
                )
        vmin, vmax = a["vmin"], a["vmax"]
        for v in (vmin, vmax):
            if v is not None and not (_r.is_number(v) and math.isfinite(v)):
                raise ValueError("'vmin' and 'vmax' must be single numbers.")
        if (vmin is not None or vmax is not None) and (
            (0 if vmin is None else vmin) >= (100 if vmax is None else vmax)
        ):
            raise ValueError("'vmin' must be smaller than 'vmax'.")

        if a["data"] is not None:
            links_data = a["data"]
        else:
            links_data = processing.get_links(
                self._state(), a["group"], a["value1"], a["value2"], a["cluster"]
            )
        # A valid alignment can have no hits; keep existing links
        if len(links_data) == 0:
            return self
        # The spec holds the links: given, or as get_links() made them. The
        # group (and value1, value2) also pick the genes labelled with the
        # measure.
        call_args = a.spec.copy()
        if a["data"] is not None:
            call_args.pop("group", None)
        call_args["data"] = links_data
        return self._with(
            _spec.add(self._spec, "links", call_args, rename={"show_links": "show"})
        )

    # chart_align.R ---------------------------------------------------------

    def align(self, id_column: str, id: Any, align: str = _MISSING) -> GC:
        """Align the clusters on a gene.

        Each cluster's scale is set so that the gene is at the same place
        in every cluster.

        Parameters
        ----------
        id_column : str
            Column of gene ids.
        id : str or float
            The id of the gene, once per cluster.
        align : str
            "left", "center" or "right": the side of the gene that is
            aligned. Default "left".

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _ALIGN)
        state = self._state()
        if state["cluster"] is None:
            _warn("Could not align selected gene. Please define cluster in GC().")
            return self
        if len(_r.as_vector(id)) != 1 or len(_r.as_vector(id_column)) != 1:
            raise ValueError("id and id_column must be single values")
        data = state["data"]
        if id_column not in data.columns:
            _warn(f"Column {id_column} not found in the data.")
            return self
        values = [_r.as_character(v) for v in data[id_column]]
        target = _r.as_character(id)
        if target not in values:
            _warn(f"ID {target} not found in the {id_column} column.")
            return self
        # The library sets the scale of each cluster from the alignment; it
        # needs the gene once per cluster
        genes = data[[v == target for v in values]]
        if genes[state["cluster"]].duplicated().any():
            _warn(
                f"The selected gene {target} appears multiple times in some clusters. "
                "Please make sure to use a unique gene id."
            )
        elif not (is_str(a["align"]) and a["align"] in ("left", "right", "center")):
            _warn(
                "Alignment should be 'left', 'right', or 'center'. "
                "Defaulting to 'left'."
            )
        return self._with(_spec.add(self._spec, "align", a.spec, mode="replace"))

    def normalize(
        self,
        group: str | None = _MISSING,
        cluster: Any = _MISSING,
        normalize: bool = _MISSING,
        preserve_gene_length: bool = _MISSING,
        gap: float | None = _MISSING,
        custom_gaps: dict[str, float] = _MISSING,
    ) -> GC:
        """Normalise the positions of the genes: even gaps between them.

        Parameters
        ----------
        group : str, optional
            Column of groups, for ``custom_gaps``. Default: the ``group`` of
            the chart.
        cluster : int, str or list, optional
            The clusters, by number (from 1) or name. Default: all.
        normalize : bool
            Whether to normalise the gaps. Default True.
        preserve_gene_length : bool
            Whether genes keep their length. Default True.
        gap : float, optional
            The gap between genes, as a share of the cluster.
        custom_gaps : dict
            Factors of the gap before genes, by group, e.g. ``{"tnsA": 2}``.

        Returns
        -------
        GC
            The chart, with the new positions in its data.
        """
        a = _spec.args(locals(), _NORMALIZE)
        return self._with(processing.normalize(self._spec, **a.values))

    # element_style.R -------------------------------------------------------

    def element_style(
        self,
        component: str = _MISSING,
        row_id: float | None = _MISSING,
        *,
        cluster: Any = _MISSING,
        part: str | None = _MISSING,
        end: str | None = _MISSING,
        reset: bool = _MISSING,
        group: Any = _MISSING,
        **options: Any,
    ) -> GC:
        """Edit or reset one element: a gene, label, title or legend entry.

        Parameters
        ----------
        component : str
            "genes", "labels", "coordinates", "clusterTitle",
            "clusterFooter", "clusterLabel", "legend", "title", "scaleBar"
            or "scale". Default "genes".
        row_id : float, optional
            The row (``rowID``, from 1) of a gene, label or coordinate.
        cluster : int or str, optional
            The cluster of a title, footer, label, scale or scale bar.
        part : str, optional
            "title" or "subtitle"; "label" or "marker" for the legend.
        end : str, optional
            "start" or "end", for a coordinate.
        reset : bool
            Remove the earlier edits of the element instead. Default False.
        group : str, optional
            The group of a legend entry.
        **options
            The element's new options, e.g. ``fill="blue"``, ``move_x=10``,
            ``text="New"``.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _ELEMENT_STYLE)
        component = a["component"]
        if component not in _ELEMENT_COMPONENTS:
            raise ValueError(f"component must be one of {quoted(_ELEMENT_COMPONENTS)}")
        state = self._state()
        row_id, cluster, part, end, group = (
            a["row_id"],
            a["cluster"],
            a["part"],
            a["end"],
            a["group"],
        )
        by_row = component in ("genes", "labels", "coordinates")
        if component == "title":
            if any(v is not None for v in (row_id, cluster, end, group)):
                raise ValueError("Use 'part' to identify a chart title or subtitle.")
            if part not in ("title", "subtitle"):
                raise ValueError('\'part\' must be "title" or "subtitle".')
        elif component == "legend":
            if any(v is not None for v in (row_id, cluster, end)):
                raise ValueError(
                    "Use 'group' and 'part' to identify a legend label or marker."
                )
            column = self._legend_group(state)
            groups = (
                [_r.as_character(v) for v in state["data"][column]]
                if column is not None and state["data"] is not None
                else []
            )
            if (
                len(_r.as_vector(group)) != 1
                or _r.is_na(group)
                or _r.as_character(group) not in groups
            ):
                raise ValueError("'group' must identify one group in the legend.")
            if part not in ("label", "marker"):
                raise ValueError('\'part\' must be "label" or "marker" for the legend.')
        elif group is not None:
            raise ValueError("'group' is only used for the legend.")
        if by_row:
            row_ids = [] if state["data"] is None else state["data"]["rowID"].tolist()
            if not (
                _r.is_number(row_id) and math.isfinite(row_id) and row_id in row_ids
            ):
                raise ValueError("'row_id' must identify one row in the chart's data.")
            if cluster is not None or part is not None:
                raise ValueError(
                    "Use 'row_id' to identify a gene, label or coordinate, "
                    "not 'cluster' or 'part'."
                )
            if end is not None and (
                component != "coordinates" or end not in ("start", "end")
            ):
                raise ValueError(
                    '\'end\' must be "start" or "end", and only for coordinates.'
                )
        elif component not in ("legend", "title"):
            if row_id is not None or end is not None:
                raise ValueError(
                    "Use 'cluster' and 'part' to identify a title or subtitle, "
                    "not 'row_id' or 'end'."
                )
            if len(_r.as_vector(cluster)) != 1 or _r.is_na(cluster):
                raise ValueError("'cluster' must be one cluster.")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                selected = _data.get_updated_clusters(
                    _data.chart_fields(self._spec)["clusters"], cluster
                )
            if selected is None:
                raise ValueError(
                    "'cluster' must be the name or number of a cluster of the chart."
                )
            if component in ("scaleBar", "scale"):
                if part is not None:
                    raise ValueError(
                        "Scale bars and axes are edited as a whole; "
                        "do not supply 'part'."
                    )
            elif part not in ("title", "subtitle"):
                raise ValueError('\'part\' must be "title" or "subtitle".')
        if not _is_flag(a["reset"]):
            raise ValueError("'reset' must be True or False.")
        style = a.options
        if any(name in ("index", "rowID", "part", "end") for name in style):
            raise ValueError(
                "Identify the element with 'row_id', 'cluster', 'part' or 'end', "
                "not with element options."
            )
        if a["reset"] and style:
            raise ValueError("Do not supply element options with reset=True.")
        if component in ("genes", "coordinates") and any(
            name in ("moveX", "moveY") for name in style
        ):
            raise ValueError(
                "Genes and coordinates cannot be moved; 'move_x' and 'move_y' are for "
                "labels, titles, scales and the legend."
            )
        return self._with(
            _spec.add_element_style(self._spec, a.spec, component, bool(a["reset"]))
        )

    def _legend_group(self, state: Mapping[str, Any]) -> str | None:
        legend = self._spec.get("legend") or {}
        group = legend.get("group") if isinstance(legend, Mapping) else None
        return group if is_str(group) else state["group"]

    # theme.R ---------------------------------------------------------------

    def theme(
        self,
        theme: str | dict[str, Any] | None = _MISSING,
        text_scale: float | None = _MISSING,
        font_family: str | None = _MISSING,
        text_color: str | None = _MISSING,
        line_color: str | None = _MISSING,
        background: str | None = _MISSING,
        palette: str | list[str] | None = _MISSING,
        pattern: dict[str, Any] | bool | None = _MISSING,
        **options: Any,
    ) -> GC:
        """Set the theme of the chart.

        A theme sets the default look of every part of the chart: fonts,
        text and line colours, the background, the gene colours and the
        options of each component. Options set with the other methods (e.g.
        ``labels(font_size=...)``) take precedence over the theme, whatever
        the order of the calls. Call it with a ``theme`` to start from that
        theme, or without one to change the current theme.

        Parameters
        ----------
        theme : str or dict, optional
            A built-in theme: "default", "classic", "dark", "vintage",
            "infographic", "publication" or "bw" (black and white, with a
            pattern per group); or theme settings (see `geneviewer_theme`).
            Default: change the current theme.
        text_scale : float, optional
            Text size, as a factor of the default font sizes.
        font_family : str, optional
            Font of all text.
        text_color : str, optional
            Colour of all text.
        line_color : str, optional
            Colour of the sequence line, the axis, the coordinates, the scale
            bar and the lines of spread labels.
        background : str, optional
            Background colour of the chart.
        palette : str or list of str, optional
            Gene colours: a d3 colour scheme, a geneviewer palette (see
            `geneviewer_palette`), or a list of colours.
        pattern : dict or bool, optional
            Patterns of the genes, with the arguments of ``pattern()``
            (``patterns``, ``color``, ``background``, ``size``,
            ``strokeWidth``). False removes the patterns of the theme.
        **options
            Options per component, as dicts named after the component:
            ``chart``, ``grid``, ``cluster``, ``title``, ``cluster_title``,
            ``cluster_footer``, ``cluster_label``, ``sequence``, ``genes``,
            ``transcript``, ``labels``, ``coordinates``, ``scale``,
            ``scale_bar``, ``legend``, ``links``, ``tooltip``,
            ``track_mouse``, ``toolbar`` and ``annotations``, e.g.
            ``labels={"fontSize": "10px"}``. None removes a component from
            the theme.

        Returns
        -------
        GC
            The chart.
        """
        a = _spec.args(locals(), _THEME)
        theme_args = {_spec.snake_to_camel(k): v for k, v in a.spec.items()}
        # None removes a component from the theme
        for name, value in a.options.items():
            if value is None:
                theme_args[name] = None
        if a["theme"] is not None:
            resolve_theme(a["theme"])
        settings_given = {
            _spec.snake_to_camel(name): a[name]
            for name in _THEME
            if name != "theme"
            and a[name] is not None
            and not (name == "pattern" and a[name] is False)
        }
        settings_given.update(a.options)
        check_theme(settings_given)
        return self._with(_spec.add_theme(self._spec, theme_args))


# The defaults of the methods' arguments, as in R (None for NULL)
_LABELS = {
    "label": None,
    "show": True,
    "cluster": None,
    "item_style": [],
    "group_style": None,
    "style_by": None,
}
_GENES = {
    "group": None,
    "marker": None,
    "marker_size": None,
    "show": True,
    "color_scheme": None,
    "custom_colors": None,
    "cluster": None,
    "item_style": [],
    "group_style": None,
    "style_by": None,
    "corner_radius": None,
}
_TRANSCRIPT = {
    "transcript": None,
    "type": None,
    "strand": None,
    "group": None,
    "selection": None,
    "show": True,
    "color_scheme": None,
    "custom_colors": None,
    "style_exons": {},
    "style_introns": {},
    "style_utrs": {},
    "item_style_exons": [],
    "item_style_introns": [],
    "item_style_utrs": [],
    "group_style_exons": None,
    "group_style_introns": None,
    "group_style_utrs": None,
    "style_by": None,
    "label_options": {},
}
_COLOR = {"color_scheme": None, "custom_colors": None}
_PATTERN = {
    "patterns": None,
    "color": "black",
    "background": "white",
    "size": 6,
    "stroke_width": 1,
}
_LEGEND = {
    "group": None,
    "show": True,
    "background_color": None,
    "order": [],
    "position": None,
    "align": None,
    "style": {},
    "legend_options": {},
    "legend_text_options": {},
    "corner_radius": None,
    "group_style": None,
}
_TITLE = {
    "title": None,
    "subtitle": None,
    "style": {},
    "subtitle_font": {},
    "title_font": {},
    "show": True,
    "height": None,
    "align": None,
    "spacing": None,
}
_CLUSTER_TITLE = {
    "title": None,
    "subtitle": None,
    "subtitle_font": {},
    "title_font": {},
    "show": True,
    "height": None,
    "cluster": None,
    "align": None,
    "spacing": None,
}
_SEQUENCE = {
    "show": True,
    "cluster": None,
    "y": None,
    "sequence_style": {},
    "marker_style": {},
}
_GRID = {
    "margin": None,
    "width": None,
    "height": None,
    "direction": "column",
    "cluster": None,
}
_SCALE = {
    "cluster": None,
    "start": None,
    "end": None,
    "padding": None,
    "hidden": False,
    "breaks": [],
    "tick_values": None,
    "reverse": False,
    "axis_position": None,
    "axis_type": "position",
    "y": None,
    "scale_breaks": False,
    "scale_break_threshold": 20,
    "scale_break_padding": 1,
    "ticks_count": None,
    "ticks_format": None,
    "tick_style": {},
    "text_style": {},
    "line_style": {},
}
_SCALE_BAR = {
    "show": True,
    "cluster": None,
    "scale_bar_line_style": {},
    "scale_bar_tick_style": {},
    "label_style": {},
}
_CLUSTER_LABEL = {
    "title": None,
    "show": True,
    "width": None,
    "cluster": None,
    "position": None,
    "align": None,
    "wrap_label": None,
    "wrap_options": {},
    "subtitle": None,
    "subtitle_font": {},
    "spacing": None,
}
_CLUSTER_FOOTER = {
    "title": None,
    "subtitle": None,
    "show": True,
    "cluster": None,
    "subtitle_font": {},
    "title_font": {},
    "align": None,
    "spacing": None,
}
_COORDINATES = {
    "show": True,
    "tick_values_top": None,
    "tick_values_bottom": None,
    "ticks_format": None,
    "tick_style": {},
    "text_style": {},
    "cluster": None,
    "group_style": None,
    "style_by": None,
}
_CLUSTER = {
    "separate_strands": None,
    "strand_spacing": None,
    "prevent_gene_overlap": None,
    "overlap_min": None,
    "overlap_spacing": None,
    "cluster": None,
    "style": {},
}
_ANNOTATION = {"type": "text", "cluster": None}
# The annotation types the library draws
ANNOTATION_TYPES = (
    "text",
    "textMarker",
    "line",
    "arrow",
    "symbol",
    "rectangle",
    "promoter",
    "terminator",
    "arc",
)
_TRACK_MOUSE = {"show": True, "cluster": None}
_EVENTS = {
    "click": True,
    "highlight": False,
    "rendered": False,
    "render_id": None,
    "dblclick": False,
    "select": False,
    "legend_click": True,
    "drag": False,
    "zoom": False,
}
_TOOLBAR = {
    "show": None,
    "position": None,
    "buttons": None,
    "size": None,
    "color": None,
    "hover_color": None,
    "tooltip_background": None,
    "tooltip_color": None,
    "tooltip_border": None,
    "tooltip_radius": None,
    "tooltip_shadow": None,
    "font_family": None,
    "icons": None,
}
_ZOOM = {
    "show": True,
    "range": None,
    "mode": None,
    "min": None,
    "update": None,
    "height": None,
    "overview": None,
    "color": None,
    "style": None,
    "window_style": None,
    "handle_style": None,
}
_TOOLTIP = {
    "formatter": "<b>Start:</b> {start}<br><b>End:</b> {end}",
    "show": True,
    "cluster": None,
}
_LINKS = {
    "group": None,
    "data": None,
    "value1": None,
    "value2": None,
    "cluster": None,
    "curve": None,
    "measure": "identity",
    "show_links": True,
    "label": True,
    "label_position": "gene",
    "normal_color": None,
    "inverted_color": None,
    "low_color": None,
    "vmin": None,
    "vmax": None,
    "use_group_colors": False,
    "color_bar": None,
    "color_bar_options": {},
    "link_width": None,
    "link_style": {},
    "label_style": {},
}
_COLOR_BAR_CHOICES = {
    "position": ("bottomright", "bottomleft", "topright", "topleft", "legend"),
    "orientation": ("vertical", "horizontal"),
    "labels": ("inline", "below"),
    "legendSide": ("left", "right"),
}
_ALIGN = {"id_column": None, "id": None, "align": "left"}
_NORMALIZE = {
    "group": None,
    "cluster": None,
    "normalize": True,
    "preserve_gene_length": True,
    "gap": None,
    "custom_gaps": {},
}
_ELEMENT_STYLE = {
    "component": "genes",
    "row_id": None,
    "cluster": None,
    "part": None,
    "end": None,
    "reset": False,
    "group": None,
}
_ELEMENT_COMPONENTS = (
    "genes",
    "labels",
    "coordinates",
    "clusterTitle",
    "clusterFooter",
    "clusterLabel",
    "legend",
    "title",
    "scaleBar",
    "scale",
)
_THEME = {
    "theme": None,
    "text_scale": None,
    "font_family": None,
    "text_color": None,
    "line_color": None,
    "background": None,
    "palette": None,
    "pattern": None,
}
