import math

import numpy as np
import pandas as pd
import pytest

from geneviewer import _r, spec
from geneviewer.spec import _MISSING


def call(defaults, **given):
    values = {name: _MISSING for name in defaults}
    options = {k: given.pop(k) for k in list(given) if k not in defaults}
    values.update(given)
    return spec.args({"self": None, **values, "options": options}, defaults)


# args ---------------------------------------------------------------------


def test_args_only_the_arguments_given():
    a = call({"label": None, "show": True, "item_style": []}, label="name")
    assert dict(a.spec) == {"label": "name"}
    assert a["show"] is True
    assert a["label"] == "name"


def test_args_none_kept_only_where_the_default_is_not_none():
    a = call({"background": "white", "patterns": None}, background=None, patterns=None)
    assert dict(a.spec) == {"background": None}


def test_args_empty_lists_left_out():
    a = call({"item_style": [], "style": {}}, item_style=[], style={})
    assert dict(a.spec) == {}


def test_args_options_keep_their_names_and_drop_none():
    a = call({"show": True}, fontSize="10px", adjust_labels=False, x=None)
    assert dict(a.spec) == {"fontSize": "10px", "adjust_labels": False}
    assert a.spec.dots == {"fontSize", "adjust_labels"}
    assert a.options == {"fontSize": "10px", "adjustLabels": False, "x": None}


def test_args_camel_case_names_of_arguments():
    a = call({"item_style": [], "row_id": None}, itemStyle=[{"index": 0}], rowID=2)
    assert dict(a.spec) == {"item_style": [{"index": 0}], "row_id": 2}
    assert a.spec.dots == set()
    with pytest.raises(TypeError, match="both"):
        call({"item_style": []}, item_style=[1], itemStyle=[2])


# names --------------------------------------------------------------------


def test_snake_and_camel():
    assert spec.snake_to_camel("normal_color") == "normalColor"
    assert spec.snake_to_camel("itemStyle") == "itemStyle"
    assert spec.camel_to_snake("legendTextOptions") == "legend_text_options"
    assert spec.camel_to_snake("rowID") == "row_id"


def test_names():
    assert spec.names(
        ["item_style", "cluster", "show_links"], {"show_links": "show"}
    ) == [
        "itemStyle",
        "clusters",
        "show",
    ]
    assert spec.names(["cluster"], clusters=False) == ["cluster"]


def test_dot_names():
    assert spec.dot_names("scaleBar", ["scale_bar_unit", "fontSize"]) == [
        "unit",
        "fontSize",
    ]
    assert spec.dot_names("labels", ["labelAdjustmentOptions"]) == ["labelAdjustment"]


# add ----------------------------------------------------------------------


def args_of(**values):
    return spec.Args(values)


def test_add_append():
    s = spec.add({}, "labels", args_of(label="name"))
    assert s["labels"] == {"label": "name"}
    s = spec.add(s, "labels", args_of(cluster=2))
    assert s["labels"] == [{"label": "name"}, {"clusters": [2]}]
    s = spec.add(s, "labels", args_of())
    assert s["labels"] == [{"label": "name"}, {"clusters": [2]}, {}]


def test_add_items():
    s = spec.add({}, "annotation", args_of(type="text"), mode="items")
    assert s["annotation"] == [{"type": "text"}]


def test_add_merge():
    s = spec.add(
        {}, "legend", args_of(position="top", style={"fill": "red"}), mode="merge"
    )
    s = spec.add(s, "legend", args_of(style={"stroke": "blue"}), mode="merge")
    assert s["legend"] == {
        "position": "top",
        "style": {"fill": "red", "stroke": "blue"},
    }


def test_add_merge_list_replaces_b64():
    s = spec.add({}, "toolbar", args_of(buttons=["undo"]), mode="merge")
    s = spec.add(s, "toolbar", args_of(buttons=["undo", "redo"]), mode="merge")
    assert s["toolbar"] == {"buttons": ["undo", "redo"]}


def test_add_merge_none_removes_later_but_is_kept_first():
    s = spec.add({}, "title", args_of(title="a", style=None), mode="merge")
    assert s["title"] == {"title": "a", "style": None}
    s = spec.add(s, "title", args_of(title="b", style=None), mode="merge")
    assert s["title"] == {"title": "b"}


def test_add_replace():
    s = spec.add({}, "color", args_of(color_scheme="a"), mode="replace")
    s = spec.add(s, "color", args_of(custom_colors=["red"]), mode="replace")
    assert s["color"] == {"customColors": ["red"]}


def test_add_clusters_and_selection_always_arrays():
    s = spec.add({}, "transcript", args_of(selection="t1", cluster=[1, 2]))
    assert s["transcript"] == {"selection": ["t1"], "clusters": [1, 2]}


def test_add_renames_and_dots():
    a = spec.Args({"scale_bar_line_style": {"stroke": "red"}, "scaleBarUnit": "kb"})
    a.dots.add("scaleBarUnit")
    s = spec.add({}, "scaleBar", a, rename={"scale_bar_line_style": "lineStyle"})
    assert s["scaleBar"] == {"lineStyle": {"stroke": "red"}, "unit": "kb"}


def test_add_does_not_change_the_spec_or_args():
    old = {"legend": {"style": {"fill": "red"}}}
    style = {"stroke": "blue"}
    new = spec.add(old, "legend", args_of(style=style), mode="merge")
    style["stroke"] = "green"
    assert old == {"legend": {"style": {"fill": "red"}}}
    assert new["legend"]["style"] == {"fill": "red", "stroke": "blue"}


# start, flag --------------------------------------------------------------


def test_start():
    data = pd.DataFrame({"start": [1], "end": [2]})
    s = spec.start(data, {"cluster": "c", "height": "200px", "style": {}})
    assert s["version"] == 1
    assert s["data"] is data
    assert s["fields"] == {"cluster": "c"}
    assert s["height"] == "200px"
    assert "style" not in s and "width" not in s and "toolbar" not in s


def test_start_without_save_button():
    s = spec.start(None, {"save_button": False})
    assert s["toolbar"] == {"buttons": ["undo", "redo", "reset", "zoom"]}


def test_flag():
    a = spec.flag(spec.Args(label=False, cluster=1), "label")
    assert dict(a) == {"cluster": 1, "show": False}
    a = spec.flag(spec.Args(label="name"), "label")
    assert dict(a) == {"label": "name"}


# zoom ---------------------------------------------------------------------


def test_add_zoom_new_mode_drops_range_and_min():
    s = spec.add_zoom({}, args_of(range=(0.2, 0.4), min=0.1))
    assert s["zoom"] == {"range": [0.2, 0.4], "min": 0.1}
    s = spec.add_zoom(s, args_of(range=[0.5, 0.6]))
    assert s["zoom"] == {"range": [0.5, 0.6], "min": 0.1}
    s = spec.add_zoom(s, args_of(mode="absolute"))
    assert s["zoom"] == {"mode": "absolute"}


# element style ------------------------------------------------------------


def test_add_element_style_and_reset():
    s = spec.add_element_style({}, args_of(row_id=1, fill="blue"), "genes", False)
    s = spec.add_element_style(
        s, args_of(component="labels", row_id=2, fill="green"), "labels", False
    )
    s = spec.add_element_style(
        s, args_of(component="labels", row_id=1, cluster=1), "labels", False
    )
    assert s["elementStyle"] == [
        {"rowID": 1, "fill": "blue"},
        {"component": "labels", "rowID": 2, "fill": "green"},
        {"component": "labels", "rowID": 1, "cluster": 1},
    ]
    r = spec.add_element_style(
        s, args_of(component="labels", row_id=2.0, reset=True), "labels", True
    )
    assert r["elementStyle"] == [s["elementStyle"][0], s["elementStyle"][2]]
    r = spec.add_element_style(r, args_of(row_id=1, reset=True), "genes", True)
    r = spec.add_element_style(
        r,
        args_of(component="labels", row_id=1, cluster="1", reset=True),
        "labels",
        True,
    )
    assert "elementStyle" not in r


# theme --------------------------------------------------------------------


def test_add_theme_changes_the_theme_set_before():
    # as R's test "theme calls without a theme change the one set before"
    s = spec.add_theme(
        {},
        {"theme": "dark", "links": {"normal_color": "red", "show_links": False}},
    )
    s = spec.add_theme(
        s,
        {
            "links": {"low_color": "white"},
            "pattern": False,
            "legend": {"legendTextOptions": {"fill": "red"}},
        },
    )
    assert s["theme"] == {
        "base": "dark",
        "links": {"normalColor": "red", "show": False, "lowColor": "white"},
        "pattern": False,
        "legend": {"textStyle": {"fill": "red"}},
    }


def test_add_theme_none_kept_and_new_theme_starts_again():
    s = spec.add_theme({}, {"theme": "classic"})
    assert s["theme"] == {"base": "classic"}
    s = spec.add_theme(s, {"labels": None, "textScale": 1.2})
    assert s["theme"] == {"base": "classic", "labels": None, "textScale": 1.2}
    s = spec.add_theme(s, {"theme": "dark"})
    assert s["theme"] == {"base": "dark"}


def test_add_theme_from_settings():
    base = {"textColor": "navy", "legend": {"textStyle": {"fill": "red"}}}
    s = spec.add_theme(
        {}, {"theme": base, "legend": {"legendTextOptions": {"fontSize": "9px"}}}
    )
    assert s["theme"] == {
        "textColor": "navy",
        "legend": {"textStyle": {"fill": "red", "fontSize": "9px"}},
    }


def test_theme_names_merges_both_names():
    settings = {
        "legend": {"textStyle": {"fill": "red"}, "legendTextOptions": {"fill": "blue"}},
        "annotations": {"text": {"font_size": "9px"}},
        "scale_bar": {"x": 1},
    }
    assert spec.theme_names(settings) == {
        "legend": {"textStyle": {"fill": "blue"}},
        "annotations": {"text": {"fontSize": "9px"}},
        "scale_bar": {"x": 1},
    }


def test_merge_settings():
    old = {"a": {"b": 1, "c": [1, 2]}, "d": 1}
    new = {"a": {"c": [3]}, "d": None}
    assert spec.merge_settings(old, new, keep_null=False) == {"a": {"b": 1, "c": [3]}}
    assert spec.merge_settings(old, new, keep_null=True) == {
        "a": {"b": 1, "c": [3]},
        "d": None,
    }
    assert old == {"a": {"b": 1, "c": [1, 2]}, "d": 1}


# sequence columns ---------------------------------------------------------


def test_drop_sequence_columns():
    data = pd.DataFrame({"start": [1], "translation": ["MK"], "sequence": ["AT"]})
    s = spec.drop_sequence_columns(
        {"data": data, "tooltip": {"formatter": "{sequence}"}}
    )
    assert list(s["data"].columns) == ["start", "sequence"]
    links = pd.DataFrame({"x": ["translation"]})
    s = spec.drop_sequence_columns({"data": data, "links": {"data": links}})
    assert list(s["data"].columns) == ["start"]


# R values -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("value", "text"),
    [
        (1, "1"),
        (1.0, "1"),
        (-2.5, "-2.5"),
        (100000.0, "1e+05"),
        (123456.0, "123456"),
        (100000, "100000"),
        (0.1, "0.1"),
        (1 / 3, "0.333333333333333"),
        (0.00012, "0.00012"),
        (1e-5, "1e-05"),
        (1.5e-7, "1.5e-07"),
        (1e15, "1e+15"),
        (np.int64(3), "3"),
        (np.float64(2.0), "2"),
        (math.inf, "Inf"),
        (True, "TRUE"),
        ("a", "a"),
        (None, None),
        (math.nan, None),
        (pd.NA, None),
    ],
)
def test_as_character(value, text):
    assert _r.as_character(value) == text
