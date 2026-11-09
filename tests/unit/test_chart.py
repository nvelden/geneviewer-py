import pandas as pd
import pytest

from geneviewer import GC, geneviewer_palette, geneviewer_theme

GENES = pd.DataFrame(
    {
        "start": [10, 90, 130, 170, 210],
        "end": [40, 120, 160, 200, 240],
        "name": ["Gene 1", "Gene 3", "Gene 4", "Gene 5", "Gene 6"],
        "group": ["A", "B", "B", "A", "C"],
        "cluster": [1, 1, 2, 2, 2],
    }
)


THEMES = ["default", "classic", "dark", "vintage", "infographic", "publication", "bw"]


def chart(**kwargs):
    return GC(GENES, cluster="cluster", group="group", **kwargs)


# GC -------------------------------------------------------------------------


def test_chart_spec():
    spec = chart(height="200px").to_spec()
    assert spec["version"] == 1
    assert spec["fields"] == {"cluster": "cluster", "group": "group"}
    assert spec["height"] == "200px"
    assert spec["data"]["name"] == GENES["name"].tolist()


def test_chart_checks_the_data():
    with pytest.raises(ValueError, match="group column not found"):
        GC(GENES, group="x")
    with pytest.raises(TypeError, match="DataFrame"):
        GC({"start": [1]})


def test_charts_are_immutable():
    base = chart()
    labelled = base.labels("name")
    assert "labels" not in base.to_spec()
    assert labelled.to_spec()["labels"] == {"label": "name"}
    GENES.loc[0, "name"] = "Changed"
    try:
        assert base.to_spec()["data"]["name"][0] == "Gene 1"
    finally:
        GENES.loc[0, "name"] = "Gene 1"


def test_snake_and_camel_case_arguments():
    a = chart().labels("name", item_style=[{"index": 0}], font_size="10px").to_spec()
    b = chart().labels("name", itemStyle=[{"index": 0}], fontSize="10px").to_spec()
    assert a == b
    assert a["labels"] == {
        "label": "name",
        "itemStyle": [{"index": 0}],
        "fontSize": "10px",
    }


def test_to_spec_validates():
    with pytest.raises(ValueError, match="not valid"):
        chart().genes(marker=3).to_spec()
    assert chart().genes(marker=3).to_spec(validate=False)["genes"] == {"marker": 3}


# checks, as in R ------------------------------------------------------------


def test_labels_checks():
    with pytest.raises(ValueError, match="Please define labels"):
        GC(GENES).labels()
    with pytest.raises(ValueError, match="label column not found"):
        chart().labels("x")
    with pytest.raises(ValueError, match="adjust_labels"):
        chart().labels(adjust_labels="no")
    with pytest.raises(ValueError, match="leader_line"):
        chart().labels(leader_line="grey")
    with pytest.warns(UserWarning, match="out of range"):
        chart().labels(cluster=3)
    assert "labels" not in chart().labels(show=False).to_spec()
    assert chart().labels(False).to_spec()["labels"] == {"show": False}


def test_legend_checks():
    with pytest.raises(ValueError, match="Please define a group"):
        GC(GENES).legend()
    assert GC(GENES).legend(False).to_spec()["legend"] == {"show": False}
    with pytest.raises(ValueError, match="align"):
        chart().legend(align="middle")


def test_title_checks():
    with pytest.raises(ValueError, match="show must be"):
        chart().title("a", show="yes")
    assert "title" not in chart().title().to_spec()


def test_color_both_warns():
    with pytest.warns(UserWarning, match="Both color_scheme"):
        out = chart().color(color_scheme="schemeSet2", custom_colors=["red"])
    assert "color" not in out.to_spec()
    spec = (
        chart().color(color_scheme="schemeSet2", custom_colors={"A": "red"}).to_spec()
    )
    assert spec["color"] == {"colorScheme": "schemeSet2", "customColors": {"A": "red"}}


def test_pattern_keeps_none_where_the_default_is_not_none():
    assert chart().pattern(background=None).to_spec()["pattern"] == {"background": None}
    with pytest.raises(ValueError, match="Unknown pattern"):
        chart().pattern(patterns=["dots", "waves"])


def test_scale_renames_and_checks():
    spec = chart().scale(scale_breaks=True, scale_break_padding="5bp").to_spec()
    assert spec["scale"] == {"autoBreaks": True, "breakPadding": "5bp"}
    with pytest.raises(ValueError, match="scale_break_padding"):
        chart().scale(scale_break_padding="5kb")


def test_grid_margins_with_percentage_height():
    with pytest.raises(ValueError, match="pixel values"):
        chart(height="100%").grid(margin={"top": "10%"})
    assert chart().grid(margin={"top": 10})._height == "100%"


def test_toolbar_and_zoom_checks():
    with pytest.raises(ValueError, match="'hover_color' must be a single colour"):
        chart().toolbar(hover_color=1)
    with pytest.raises(ValueError, match="'buttons'"):
        chart().toolbar(buttons=["print"])
    assert chart().toolbar(buttons="save").to_spec()["toolbar"] == {"buttons": ["save"]}
    with pytest.raises(ValueError, match="relative 'range'"):
        chart().zoom(range=[10, 20])
    spec = chart().zoom(range=(20, 10), mode="absolute").to_spec()
    assert spec["zoom"] == {"range": [20, 10], "mode": "absolute"}


def test_events_checks():
    with pytest.raises(ValueError, match="True or False"):
        chart().events(click="yes")
    spec = chart().events(select=True, legend_click=False).to_spec()
    assert spec["events"] == {"select": True, "legendClick": False}


def test_links_checks():
    with pytest.raises(ValueError, match="Invalid measure"):
        chart().links(data=GENES, measure="score")
    with pytest.raises(ValueError, match="'vmin' must be smaller"):
        chart().links(data=GENES, vmin=50, vmax=10)
    with pytest.raises(NotImplementedError):
        chart().links("group")


def test_align_warnings():
    with pytest.warns(UserWarning, match="define cluster"):
        GC(GENES).align("name", "Gene 3")
    with pytest.warns(UserWarning, match="not found"):
        chart().align("name", "Gene 9")
    spec = chart().align("name", "Gene 3", align="center").to_spec()
    assert spec["align"] == {"idColumn": "name", "id": "Gene 3", "align": "center"}


def test_element_style_checks_and_reset():
    with pytest.raises(ValueError, match="row_id"):
        chart().element_style("genes", row_id=9, fill="red")
    with pytest.raises(ValueError, match="cannot be moved"):
        chart().element_style("genes", row_id=1, move_x=3)
    with pytest.raises(ValueError, match="one group in the legend"):
        chart().element_style("legend", group="Z", part="label")
    spec = (
        chart()
        .element_style("genes", row_id=1, fill="blue")
        .element_style("labels", row_id=2, fill="green")
        .element_style("labels", row_id=2, reset=True)
        .element_style("clusterTitle", cluster=2, part="title", text="Two")
        .to_spec()
    )
    assert spec["elementStyle"] == [
        {"component": "genes", "rowID": 1, "fill": "blue"},
        {"component": "clusterTitle", "cluster": 2, "part": "title", "text": "Two"},
    ]


# theme ------------------------------------------------------------------------


def test_theme_spec():
    spec = (
        chart()
        .theme("dark", links={"normal_color": "red", "show_links": False})
        .theme(
            links={"low_color": "white"},
            pattern=False,
            legend={"legendTextOptions": {"fill": "red"}},
            cluster_title=None,
        )
        .to_spec()
    )
    assert spec["theme"] == {
        "base": "dark",
        "links": {"normalColor": "red", "show": False, "lowColor": "white"},
        "pattern": False,
        "legend": {"textStyle": {"fill": "red"}},
        "clusterTitle": None,
    }


def test_theme_checks():
    with pytest.raises(ValueError, match="name must be one of"):
        chart().theme("neon")
    with pytest.raises(ValueError, match="Unknown theme setting"):
        chart().theme(labelz={"fill": "red"})
    with pytest.raises(ValueError, match="text_scale"):
        chart().theme(text_scale=-1)


def test_geneviewer_theme_returns_a_new_copy():
    first = geneviewer_theme("classic")
    first["textColor"] = "navy"
    first["genes"]["stroke"] = "red"
    assert geneviewer_theme("classic") != first


@pytest.mark.parametrize("name", THEMES)
def test_every_theme_as_settings_gives_a_valid_spec(name):
    spec = chart().theme(**geneviewer_theme(name)).to_spec()
    assert "base" not in spec["theme"]
    spec = chart().theme(geneviewer_theme(name), text_scale=1.2).to_spec()
    assert spec["theme"]["textScale"] == 1.2


def test_geneviewer_palette():
    assert geneviewer_palette() == [
        "okabeIto",
        "classic",
        "vivid",
        "dark",
        "vintage",
        "infographic",
        "greys",
        "blues",
    ]
    assert len(geneviewer_palette("okabeIto")) == 8
    with pytest.raises(ValueError, match="name must be one of"):
        geneviewer_palette("rainbow")
