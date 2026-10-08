#!/usr/bin/env python

# Test the learn.plotting module

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns
from Bio import AlignIO
from matplotlib.container import BarContainer, ErrorbarContainer
from scipy import stats

from teemi.learn import plotting
from teemi.learn.plotting import (
    bar_plot,
    bar_plot_w_hue,
    carpet_barplot,
    color_range_dict,
    correlation_plot,
    grouped_bar_plot,
    horisontal_bar_plot,
    plot_ml_learning_curve,
    plot_phylo_tree,
    plot_stacked_barplot_with_labels,
)


@pytest.fixture(autouse=True)
def shown_figures(monkeypatch):
    """Replace plt.show with a recorder, keep rcParams changes made by the
    plotting functions local to the test and close all figures afterwards."""
    shown = []
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: shown.append(plt.gcf()))
    with mpl.rc_context():
        yield shown
    plt.close("all")


def hex_color(color):
    return mpl.colors.to_hex(color)


def dashed_lines(ax, value, axis):
    """Dashed reference lines drawn at ``value`` (axhline -> axis='y')."""
    found = []
    for line in ax.lines:
        data = line.get_ydata() if axis == "y" else line.get_xdata()
        if line.get_linestyle() == "--" and list(data) == [value, value]:
            found.append(line)
    return found


def spines_visible(ax):
    return {name: spine.get_visible() for name, spine in ax.spines.items()}


NO_SPINES = {"left": False, "right": False, "top": False, "bottom": False}


# ---------------------------------------------------------------- carpet plot
def test_carpet_barplot(shown_figures, tmp_path):
    crosstab = pd.DataFrame(
        {"pG8H": [0.6, 0.3, 0.1], "pCPR": [0.4, 0.7, 0.9]},
        index=pd.Index([1, 2, 3], name="position"),
    )
    color_dict = {"pG8H": "#ff0000", "pCPR": "#0000ff"}

    carpet_barplot(
        crosstab,
        color_dict,
        path=str(tmp_path / "carpet"),
        xlabel="Position",
        ylabel="Proportion",
        size_height=4,
        size_length=6,
        bar_width=0.5,
    )

    assert (tmp_path / "carpet.pdf").is_file()
    assert len(shown_figures) == 1
    fig = shown_figures[0]
    ax = fig.axes[0]

    # One stacked bar container per column, coloured from the color dict
    bars = [c for c in ax.containers if isinstance(c, BarContainer)]
    assert [c.get_label() for c in bars] == ["pG8H", "pCPR"]
    for container in bars:
        for patch in container:
            assert hex_color(patch.get_facecolor()) == color_dict[container.get_label()]
            assert patch.get_width() == pytest.approx(0.5)
    assert [p.get_height() for p in bars[0]] == pytest.approx([0.6, 0.3, 0.1])
    # second layer is stacked on top of the first
    assert [p.get_y() for p in bars[1]] == pytest.approx([0.6, 0.3, 0.1])

    assert ax.get_xlabel() == "Position"
    assert ax.get_ylabel() == "Proportion"
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["pG8H", "pCPR"]
    assert len(ax.get_yticks()) == 0
    assert spines_visible(ax) == NO_SPINES
    assert tuple(fig.get_size_inches()) == (6, 4)


def test_carpet_barplot_does_not_save_without_path(monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)
    crosstab = pd.DataFrame({"a": [1, 2], "b": [2, 1]})

    carpet_barplot(crosstab, {"a": "red", "b": "blue"})
    carpet_barplot(crosstab, {"a": "red", "b": "blue"}, save_pdf=False, path="carpet")

    assert list(tmp_path.iterdir()) == []
    assert len(shown_figures) == 2


# ------------------------------------------------------------- learning curve
def test_plot_ml_learning_curve(shown_figures, tmp_path):
    x = [10, 20, 30]
    y_training = [5.0, 4.0, 3.0]
    y_cv = [8.0, 6.0, 5.0]
    training_sd = [0.5, 0.4, 0.3]
    cv_sd = [1.0, 0.8, 0.5]

    plot_ml_learning_curve(
        x,
        y_training,
        y_cv,
        training_sd,
        cv_sd,
        path=str(tmp_path / "curve"),
        size_height=4,
        size_length=7,
        title="Learning curve",
        y_axis_range=[0, 12],
    )

    assert (tmp_path / "curve.pdf").is_file()
    fig = shown_figures[0]
    ax = fig.axes[0]

    cv_line, training_line = ax.lines
    assert list(cv_line.get_xdata()) == x
    assert list(cv_line.get_ydata()) == y_cv
    assert hex_color(cv_line.get_color()) == "#986e42"
    assert list(training_line.get_ydata()) == y_training
    assert hex_color(training_line.get_color()) == "#0000ff"

    # The two standard-deviation bands span mean -/+ sd
    cv_band, training_band = ax.collections
    cv_vertices = cv_band.get_paths()[0].vertices
    assert cv_vertices[:, 1].min() == pytest.approx(min(np.subtract(y_cv, cv_sd)))
    assert cv_vertices[:, 1].max() == pytest.approx(max(np.add(y_cv, cv_sd)))
    assert hex_color(cv_band.get_facecolor()[0]) == "#ffe7b5"
    training_vertices = training_band.get_paths()[0].vertices
    assert training_vertices[:, 1].min() == pytest.approx(2.7)
    assert training_vertices[:, 1].max() == pytest.approx(5.5)
    assert hex_color(training_band.get_facecolor()[0]) == "#add8e6"

    # Legend entries are attached to the matching artists
    legend = ax.get_legend()
    assert [t.get_text() for t in legend.get_texts()] == [
        "Cross-validation mean MAE",
        "Cross-validation standard-deviation",
        "Model performance MAE",
        "Model standard-deviation",
    ]
    handles = legend.legend_handles
    assert hex_color(handles[0].get_color()) == "#986e42"
    assert hex_color(handles[1].get_facecolor()) == "#ffe7b5"
    assert hex_color(handles[2].get_color()) == "#0000ff"
    assert hex_color(handles[3].get_facecolor()) == "#add8e6"

    assert ax.get_title() == "Learning curve"
    assert ax.get_xlabel() == "Length of the partitioned data"
    assert ax.get_ylabel() == "MAE"
    assert ax.get_ylim() == (0, 12)
    assert spines_visible(ax) == NO_SPINES
    assert tuple(fig.get_size_inches()) == (7, 4)


def test_plot_ml_learning_curve_does_not_save_without_path(monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)

    plot_ml_learning_curve([1, 2], [1, 1], [2, 2], [0, 0], [0, 0])
    plot_ml_learning_curve([1, 2], [1, 1], [2, 2], [0, 0], [0, 0], save_pdf=False, path="c")

    assert list(tmp_path.iterdir()) == []
    assert len(shown_figures) == 2


# ------------------------------------------------------------------- bar plot
def test_bar_plot_with_all_options(shown_figures, tmp_path):
    bar_plot(
        ["G8H", "CPR", "WT"],
        [90, 110, 100],
        error_bar=[5, 10, 2],
        color="#deebf7",
        path=str(tmp_path / "bars"),
        title="Titer",
        x_label="Strain",
        y_label="Relative titer (%)",
        size_height=5,
        size_length=3,
    )

    assert (tmp_path / "bars.pdf").is_file()
    fig = shown_figures[0]
    ax = fig.axes[0]

    assert [p.get_height() for p in ax.patches] == [90, 110, 100]
    for patch in ax.patches:
        assert hex_color(patch.get_facecolor()) == "#deebf7"
        assert hex_color(patch.get_edgecolor()) == "#000000"
    assert [t.get_text() for t in ax.get_xticklabels()] == ["G8H", "CPR", "WT"]

    errorbars = [c for c in ax.containers if isinstance(c, ErrorbarContainer)]
    assert len(errorbars) == 1
    _, _, (vertical_lines,) = errorbars[0].lines
    segments = vertical_lines.get_segments()
    assert [seg[:, 1].tolist() for seg in segments] == [[85, 95], [100, 120], [98, 102]]

    assert len(dashed_lines(ax, 100, axis="y")) == 1
    assert ax.get_title() == "Titer"
    assert ax.get_xlabel() == "Strain"
    assert ax.get_ylabel() == "Relative titer (%)"
    assert spines_visible(ax) == NO_SPINES
    assert tuple(fig.get_size_inches()) == (3, 5)


def test_bar_plot_minimal(monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)

    bar_plot(["a", "b"], [1, 2], horisontal_line=False, save_pdf=False, path="bars")

    ax = shown_figures[0].axes[0]
    assert [p.get_height() for p in ax.patches] == [1, 2]
    assert not any(isinstance(c, ErrorbarContainer) for c in ax.containers)
    assert len(ax.lines) == 0
    assert ax.get_title() == ""
    assert ax.get_xlabel() == ""
    assert ax.get_ylabel() == ""
    assert list(tmp_path.iterdir()) == []


# -------------------------------------------------------- horizontal bar plot
def test_horisontal_bar_plot_with_all_options(shown_figures, tmp_path):
    horisontal_bar_plot(
        ["G8H", "CPR"],
        [80, 120],
        path=str(tmp_path / "hbars"),
        color="#fee6ce",
        title="Titer",
        x_label="Relative titer (%)",
        y_label="Strain",
        size_height=6,
        size_length=6,
        legend=True,
    )

    assert (tmp_path / "hbars.pdf").is_file()
    fig = shown_figures[0]
    ax = fig.axes[0]

    assert [p.get_width() for p in ax.patches] == [80, 120]
    for patch in ax.patches:
        assert hex_color(patch.get_facecolor()) == "#fee6ce"
    # every bar is annotated with its value
    assert [t.get_text() for t in ax.texts] == ["80", "120"]

    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["dbtl_2", "dbtl_1"]
    assert len(dashed_lines(ax, 100, axis="x")) == 1
    assert ax.get_title() == "Titer"
    assert ax.get_xlabel() == "Relative titer (%)"
    assert ax.get_ylabel() == "Strain"
    assert spines_visible(ax) == NO_SPINES
    assert tuple(fig.get_size_inches()) == (6, 6)


def test_horisontal_bar_plot_minimal(monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)

    horisontal_bar_plot(["a", "b"], [1, 2], vertical_line=False, save_pdf=False, path="h")

    ax = shown_figures[0].axes[0]
    assert [p.get_width() for p in ax.patches] == [1, 2]
    assert ax.get_legend() is None
    assert dashed_lines(ax, 100, axis="x") == []
    assert ax.get_title() == ""
    assert list(tmp_path.iterdir()) == []


# ----------------------------------------------------------- correlation plot
@pytest.fixture
def correlation_df():
    return pd.DataFrame(
        {
            "predicted": [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80],
            "measured": [0.15, 0.18, 0.35, 0.38, 0.52, 0.58, 0.75, 0.79],
        }
    )


def test_correlation_plot(correlation_df, shown_figures, tmp_path):
    r, p = stats.pearsonr(correlation_df["predicted"], correlation_df["measured"])

    correlation_plot(
        correlation_df.copy(),
        "predicted",
        "measured",
        path=str(tmp_path / "corr"),
        title="Prediction vs measurement",
        size_height=5,
        size_length=5,
        x_axis_range=[0, 0.9],
        y_axis_range=[0.1, 0.8],
    )

    assert (tmp_path / "corr.pdf").is_file()
    fig = shown_figures[0]
    ax = fig.axes[0]

    # the data points are drawn as a black scatter
    points = ax.collections[0]
    np.testing.assert_allclose(points.get_offsets(), correlation_df.to_numpy())
    assert hex_color(points.get_facecolor()[0]) == "#000000"

    # the fitted regression line follows the least squares fit
    regression_line = ax.lines[0]
    slope, intercept = np.polyfit(correlation_df["predicted"], correlation_df["measured"], 1)
    x_line = np.asarray(regression_line.get_xdata())
    np.testing.assert_allclose(regression_line.get_ydata(), slope * x_line + intercept)
    assert hex_color(regression_line.get_color()) == "#000000"

    # the statistics are reported in the figure title
    suptitle = fig._suptitle.get_text()
    assert suptitle.startswith("R-squared = ")
    assert f"P-value = {p:.3E}" in suptitle

    assert ax.get_title() == "Prediction vs measurement"
    assert ax.get_xlim() == (0, 0.9)
    assert ax.get_ylim() == (0.1, 0.8)
    assert tuple(fig.get_size_inches()) == (5, 5)


def test_correlation_plot_does_not_save_without_path(correlation_df, monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)

    correlation_plot(correlation_df.copy(), "predicted", "measured")
    correlation_plot(correlation_df.copy(), "predicted", "measured", save_pdf=False, path="c")

    assert list(tmp_path.iterdir()) == []
    assert len(shown_figures) == 2


# ----------------------------------------------------------- bar plot w. hue
@pytest.fixture
def hue_df():
    return pd.DataFrame(
        {
            "strain": ["s1", "s1", "s2", "s2", "s1", "s1", "s2", "s2"],
            "titer": [90, 110, 40, 60, 140, 160, 70, 90],
            "category": ["dbtl_1"] * 4 + ["dbtl_2"] * 4,
        }
    )


def test_bar_plot_w_hue(hue_df, shown_figures, tmp_path):
    bar_plot_w_hue(
        hue_df,
        "strain",
        "titer",
        path=str(tmp_path / "hue"),
        palette=["#1f77b4", "#ff7f0e"],
        title="Titer",
        x_label="Strain",
        y_label="Titer (%)",
        size_height=4,
        size_length=6,
    )

    assert (tmp_path / "hue.pdf").is_file()
    fig = shown_figures[0]
    ax = fig.axes[0]

    # one bar container per hue level, bar height is the group mean
    bars = [c for c in ax.containers if isinstance(c, BarContainer)]
    assert len(bars) == 2
    assert [p.get_height() for p in bars[0]] == pytest.approx([100, 50])
    assert [p.get_height() for p in bars[1]] == pytest.approx([150, 80])
    # seaborn draws the palette colours slightly desaturated
    assert {hex_color(p.get_facecolor()) for p in bars[0]} == {
        hex_color(sns.desaturate("#1f77b4", 0.75))
    }
    assert {hex_color(p.get_facecolor()) for p in bars[1]} == {
        hex_color(sns.desaturate("#ff7f0e", 0.75))
    }
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["dbtl_1", "dbtl_2"]

    assert len(dashed_lines(ax, 100, axis="y")) == 1
    assert ax.get_title() == "Titer"
    assert ax.get_xlabel() == "Strain"
    assert ax.get_ylabel() == "Titer (%)"
    assert not ax.spines["top"].get_visible()
    assert not ax.spines["right"].get_visible()
    assert tuple(fig.get_size_inches()) == (6, 4)


def test_bar_plot_w_hue_without_line_and_pdf(hue_df, monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)

    bar_plot_w_hue(hue_df, "strain", "titer", horisontal_line=False, save_pdf=False, path="h")

    ax = shown_figures[0].axes[0]
    assert dashed_lines(ax, 100, axis="y") == []
    assert list(tmp_path.iterdir()) == []


# --------------------------------------------------------------- color ranges
def test_color_range_dict():
    colors = color_range_dict()

    assert list(colors) == ["yellow", "orange", "blue", "green"]
    assert {name: len(values) for name, values in colors.items()} == {
        "yellow": 19,
        "orange": 22,
        "blue": 24,
        "green": 18,
    }
    for values in colors.values():
        assert len(set(values)) == len(values)
        for value in values:
            assert value.startswith("#") and len(value) == 7
            assert mpl.colors.is_color_like(value)
    # every call returns a fresh dictionary
    assert color_range_dict() is not colors


# ---------------------------------------------------------------- phylo tree
@pytest.fixture
def alignment(tmp_path):
    alignment_file = tmp_path / "homologs.fasta"
    alignment_file.write_text(
        ">G8H_A\nMKTAYIAKQRQISFVK\n"
        ">G8H_B\nMKTAYIAKQRQISFVR\n"
        ">G8H_C\nMKSAYLAEQRQLSFAK\n"
        ">G8H_D\nGGSWYLAEHHNLTWAE\n"
    )
    return AlignIO.read(alignment_file, "fasta")


def test_plot_phylo_tree(alignment, monkeypatch, shown_figures, tmp_path):
    drawn_trees = []
    original_draw = plotting.Phylo.draw

    def spy_draw(tree, *args, **kwargs):
        drawn_trees.append(tree)
        return original_draw(tree, *args, **kwargs)

    monkeypatch.setattr(plotting.Phylo, "draw", spy_draw)

    plot_phylo_tree(alignment, path=str(tmp_path / "tree"), height=6, wideness=6)

    assert (tmp_path / "tree.pdf").is_file()

    # UPGMA on identity distances groups the two near identical sequences
    (tree,) = drawn_trees
    assert sorted(leaf.name for leaf in tree.get_terminals()) == [
        "G8H_A",
        "G8H_B",
        "G8H_C",
        "G8H_D",
    ]
    assert len(tree.common_ancestor("G8H_A", "G8H_B").get_terminals()) == 2

    fig = shown_figures[0]
    assert fig.dpi == 300
    assert tuple(fig.get_size_inches()) == (6, 6)
    ax = fig.axes[0]
    assert not ax.axison
    labels = {t.get_text().strip() for t in ax.texts}
    assert {"G8H_A", "G8H_B", "G8H_C", "G8H_D"} <= labels


def test_plot_phylo_tree_does_not_save_without_path(alignment, monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)
    before = sorted(tmp_path.iterdir())

    plot_phylo_tree(alignment)
    plot_phylo_tree(alignment, save_pdf=False, path="tree")

    assert sorted(tmp_path.iterdir()) == before
    assert len(shown_figures) == 2


# ------------------------------------------------- stacked barplot w. labels
@pytest.fixture
def occurrences_df():
    return pd.DataFrame(
        {"pG8H": [25.0, 62.5], "pCPR": [75.0, 37.5]},
        index=["Promoter 1", "Promoter 2"],
    )


def test_plot_stacked_barplot_with_labels(occurrences_df, tmp_path):
    plot_stacked_barplot_with_labels(
        occurrences_df,
        ["#e6550d", "#3182bd"],
        title="Occurrences",
        path=str(tmp_path) + "/",
        size_length=8,
        size_heigth=4,
    )

    assert (tmp_path / "Occurences of each part sampled.pdf").is_file()
    fig = plt.gcf()
    ax = fig.axes[0]

    bars = [c for c in ax.containers if isinstance(c, BarContainer)]
    assert [c.get_label() for c in bars] == ["pG8H", "pCPR"]
    assert {hex_color(p.get_facecolor()) for p in bars[0]} == {"#e6550d"}
    assert {hex_color(p.get_facecolor()) for p in bars[1]} == {"#3182bd"}
    assert {hex_color(p.get_edgecolor()) for c in bars for p in c} == {"#000000"}

    # each box is labelled with its column name and percentage
    assert [t.get_text() for t in ax.texts] == [
        "pG8H \n25.0 %",
        "pG8H \n62.5 %",
        "pCPR \n75.0 %",
        "pCPR \n37.5 %",
    ]
    assert ax.get_title() == "Occurrences"
    assert ax.get_legend().get_texts() == []
    assert spines_visible(ax) == NO_SPINES
    assert tuple(fig.get_size_inches()) == (8, 4)


def test_plot_stacked_barplot_with_labels_does_not_save_without_path(occurrences_df, monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    plot_stacked_barplot_with_labels(occurrences_df, ["red", "blue"])

    assert list(tmp_path.iterdir()) == []


# ----------------------------------------------------------- grouped barplot
def test_grouped_bar_plot(shown_figures, tmp_path):
    grouped_bar_plot(
        ["G8H\nWT", "G8H\nmut", "CPR\nWT", "CPR\nmut"],
        [100, 80, 100, 120],
        ["white", "black", "white", "black"],
        ["Wild type", "Mutant"],
        title="Titer",
        y_label="Relative titer (%)",
        x_label="Strain",
        path=str(tmp_path) + "/",
        size_height=6,
        size_length=6,
    )

    assert (tmp_path / "grouped_bar_plot.pdf").is_file()
    fig = shown_figures[0]
    ax = fig.axes[0]

    assert [p.get_height() for p in ax.patches] == [100, 80, 100, 120]
    assert [hex_color(p.get_facecolor()) for p in ax.patches] == [
        "#ffffff",
        "#000000",
        "#ffffff",
        "#000000",
    ]

    legend = ax.get_legend()
    assert [t.get_text() for t in legend.get_texts()] == ["Wild type", "Mutant"]
    assert [hex_color(h.get_facecolor()) for h in legend.legend_handles] == [
        "#ffffff",
        "#000000",
    ]
    assert not legend.get_frame_on()

    assert len(dashed_lines(ax, 100, axis="y")) == 1
    assert ax.get_title() == "Titer"
    assert ax.get_xlabel() == "Strain"
    assert ax.get_ylabel() == "Relative titer (%)"
    assert spines_visible(ax) == NO_SPINES
    assert tuple(fig.get_size_inches()) == (6, 6)


def test_grouped_bar_plot_without_line_and_pdf(monkeypatch, tmp_path, shown_figures):
    monkeypatch.chdir(tmp_path)

    grouped_bar_plot(["a", "b"], [1, 2], ["white", "black"], ["x", "y"], axhline=False)

    ax = shown_figures[0].axes[0]
    assert dashed_lines(ax, 100, axis="y") == []
    assert list(tmp_path.iterdir()) == []
