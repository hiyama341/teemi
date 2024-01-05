#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

""" Module used for producing/reproducing plots"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy import stats
import pandas as pd
import matplotlib.patches as mpatches


# for phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo
import pylab
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


def carpet_barplot(
    pd_dataframe_cross_tab_prop,
    colorDict: dict,
    save_pdf=True,
    path="",
    xlabel="",
    ylabel="",
    size_height: int = 10,
    size_length: int = 20,
    bar_width=1.0,
    ylim=[0, 25],
) -> None:
    """Plotting stacked barplots from a pandas dataframe cross tab df

    Parameters
    ----------
    pd_dataframe_cross_tab_prop : pd.DataFrames
        pd.crosstab
    colorDict : Dict
        key as coloumn name, value hex color code
    save_pdf : bool
    path : string

    Returns
    -------
    stacked barplot
    """

    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Plot
    ax = pd_dataframe_cross_tab_prop.plot(
        kind="bar", stacked=True, color=colorDict, width=bar_width
    )

    plt.legend(loc="upper right", ncol=2)
    plt.xlabel(xlabel, size=25, fontname="Helvetica", fontweight="bold")
    plt.ylabel(ylabel, size=25, fontname="Helvetica", fontweight="bold")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    plt.xscale("linear")

    # make axis integers
    plt.locator_params(axis="both", integer=True, tight=True)

    # remove y axis labes
    plt.yticks([])
    plt.xticks(rotation=0, fontweight="bold", size=20)  # changing x scale by own

    ## size matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_length, size_height)

    # tight layout
    plt.tight_layout()

    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    if save_pdf and path != "":
        ## save pdf
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def plot_ml_learning_curve(
    x_partitioned_data: list,
    y_training: list,
    y_cv: list,
    training_sd,
    cv_sd: list,
    save_pdf=True,
    path="",
    size_height: int = 10,
    size_length: int = 10,
    title="",
    linewidth: int = 1.5,
    y_axis_range: list = [0, 25],
) -> None:
    """Plotting a learning curve from partitioned dataframes.

    Parameters
    ----------
    x_partitioned_data : list
        x coordinates
    y_training : list
        y coordinates for trained models
    y_cv : list
        y coordinates for cross-validated models
    training_sd : list
        standard-deviation for the models
    cv_sd : list
        standard-deviation for cross-validated models

    save_pdf : bool
    path : str

    Returns
    -------
    Learning curve
    """
    # making sd into numpy array
    training_sd = np.array(training_sd)
    cv_sd = np.array(cv_sd)

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # CV_mae
    plt.plot(x_partitioned_data, y_cv, color="#986e42", linewidth=linewidth)
    plt.fill_between(x_partitioned_data, y_cv - cv_sd, y_cv + cv_sd, color="#ffe7b5")

    # Model plotted
    plt.plot(x_partitioned_data, y_training, color="Blue", linewidth=linewidth)
    plt.fill_between(
        x_partitioned_data,
        y_training - training_sd,
        y_training + training_sd,
        color="#ADD8E6",
    )

    # Remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Add labels and titel
    ax.set_xlabel(
        "Length of the partitioned data",
        size=20,
        fontname="Helvetica",
        fontweight="bold",
    )
    ax.set_ylabel("MAE", size=20, fontname="Helvetica")
    ax.set_title(title, size=30, fontname="Helvetica")
    ax.legend(
        [
            "Cross-validation mean MAE",
            "Cross-validation standard-deviation",
            "Model performance MAE",
            "Model standard-deviation",
        ],
        loc="upper right",
        shadow=True,
        fontsize="small",
    )

    # Set color
    ax.set_facecolor("white")

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_length, size_height)

    # change y-axis range
    plt.ylim(y_axis_range)

    if save_pdf and path != "":
        ## save pdf
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    # show
    plt.show()


def bar_plot(
    x: list,
    y: list,
    error_bar: list = None,
    horisontal_line=True,
    save_pdf=True,
    color="white",
    path="",
    title=None,
    x_label=None,
    y_label=None,
    size_height: int = 25,
    size_length: int = 15,
) -> None:
    """Plotting a bar_plot .

    Parameters
    ----------
    x : list
        x coordinates
    y : list
        y coordinates
    error_bar : list  (Optional)
        lits of errorbars matching the x coordinates
    horisontal_line : bool
    save_pdf : bool
    path : str
    color : str
        can be matplotlib color names or hex color codes
    title : str
    x_label : str
    y_label : str

    Returns
    -------
    bar_plot"""

    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # Plot
    plt.bar(x, y, edgecolor="black", color=color)  # white

    # Errorbar
    if error_bar != None:
        plt.errorbar(
            x,
            y,
            yerr=error_bar,
            fmt="o",
            color="black",
        )  # ms = 2)

    # add horisontal line
    if horisontal_line:
        plt.axhline(y=100, color="black", linestyle="dashed")

    # Title and labels
    if title is not None:
        ax.set_title(title, size=30, fontname="Helvetica")
    if x_label is not None:
        ax.set_xlabel(x_label, size=20, fontname="Helvetica")
    if y_label is not None:
        ax.set_ylabel(y_label, size=20, fontname="Helvetica")

    # Set color
    ax.set_facecolor("white")

    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_length, size_height)

    if save_pdf and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def horisontal_bar_plot(
    x: list,
    y: list,
    vertical_line: bool = True,
    save_pdf: bool = True,
    path: str = "",
    color="white",
    title=None,
    x_label=None,
    y_label=None,
    size_height: int = 10,
    size_length: int = 8,
    legend=False,
) -> None:
    """Plotting a horisontal_bar_plot .

    Parameters
    ----------
    x : list
        x coordinates
    y : list
        y coordinates
    vertical_line : bool
    save_pdf : bool
    path : str
    color : str
    can be matplotlib color names or hex color codes
    title : str
    x_label : str
    y_label : str

    Returns
    -------
    horisontal_bar_plot"""
    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # Plot
    plt.barh(x, y, edgecolor="black", color=color)  # white # '#deebf7'

    # Change x labels rotation
    ax.tick_params(rotation=90)

    # LEGEND
    if legend:
        (lineOne,) = ax.plot([], label="dbtl_2")
        (lineTwo,) = ax.plot([], label="dbtl_1")
        ax.legend(handles=[lineOne, lineTwo], loc="best")

    # remove gridlines
    ax.grid(False)

    # Title and labels
    if title is not None:
        ax.set_title(title, size=30, fontname="Helvetica")
    if x_label is not None:
        ax.set_xlabel(x_label, size=20, fontname="Helvetica")
    if y_label is not None:
        ax.set_ylabel(y_label, size=20, fontname="Helvetica")

    # Background
    ax.set_facecolor("white")

    # add horisontal line
    if vertical_line:
        plt.axvline(x=100, color="black", linestyle="dashed")

    # rotate sticks
    ax.tick_params(rotation=0)

    # ## adding the labels on the bar
    for c in ax.containers:
        ax.bar_label(c, padding=10)
    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_height, size_length)

    if save_pdf and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def correlation_plot(
    dataframe,
    x: str,
    y: str,
    save_pdf: bool = True,
    path: str = "",
    title: str = "",
    size_height: int = 10,
    size_length: int = 10,
    y_axis_range: list = [0, 1],
    x_axis_range: list = [0, 1],
) -> None:
    """Plotting a correlation_plot.

    Parameters
    ----------
    dataframe : pd.DataFrame
    x : str
        x coordinates
    y : str
        y coordinates
    save_pdf : bool
    path : str
    size_height : int
    size_length : int
    x_axis_range : list
    y_axis_range : list

    Returns
    -------
    correlation_plot"""

    # set seaborn plotting aesthetics as default
    sns.set_context("paper", font_scale=2.0, rc={"lines.linewidth": 1.5})
    dataframe["color"] = "black"
    g = sns.lmplot(
        data=dataframe,
        x=x,
        y=y,
        hue="color",
        palette=["#000000"],
        fit_reg=True,
        height=10,
        line_kws={"color": "black"},
        ci=False,
        legend=False,
    )
    r, p = stats.pearsonr(dataframe[x], dataframe[y])
    # add white background
    ax = plt.gca()
    ax.set_facecolor("white")
    # add R and P values
    plt.suptitle(
        f"R-squared = {r:.3f} \n  P-value = {p:.3E}",
        y=0.8,
        x=0.3,
        size=15,
        fontname="Helvetica",
        fontweight="bold",
    )

    # add title
    ax.set_title(title)

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_height, size_length)

    # change y-axis range
    plt.xlim(x_axis_range)
    plt.ylim(y_axis_range)

    if save_pdf and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def bar_plot_w_hue(
    dataframe,
    x: str,
    y: str,
    save_pdf=True,
    path="",
    hue: str = "category",
    palette="dark",
    title="",
    x_label="",
    y_label="",
    horisontal_line: bool = True,
    size_height: int = 10,
    size_length: int = 10,
) -> None:
    """Plotting a correlation_plot.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Dataframe with categories in seperate column from x and y
    x : str
        x coordinates
    y : str
        y coordinates
    save_pdf : bool
    path : str
        path to folder with name of the file
    title : str
    x_label : str
    y_label : str

    Returns
    -------
    bar_plot_w_hue"""
    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    ax = sns.barplot(x=x, y=y, hue=hue, data=dataframe, palette=palette)
    # Remove spines
    sns.despine()
    # add labels
    ax = plt.gca()
    ax.set_xlabel(x_label, size=20, fontname="Helvetica", fontweight="bold")
    ax.set_ylabel(y_label, size=20, fontname="Helvetica", fontweight="bold")
    ax.set_title(title, size=30, fontname="Helvetica", fontweight="bold")

    # white background
    ax.set_facecolor("white")
    plt.xscale("linear")

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_length, size_height)

    if horisontal_line:
        # normalized line
        ax.axhline(100, color="black", linestyle="dashed")

    if save_pdf == True and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300)

    plt.show()


def color_range_dict() -> dict:
    """Returns a dictionary of color ranges.

    Returns
    -------
    dict
        A dictionary of color ranges containing keys 'yellow', 'orange', 'blue' and 'green'
    """
    return {
        "yellow": [
            "#a59a00",
            "#aa9e00",
            "#b0a300",
            "#b5a800",
            "#bbad00",
            "#c0b200",
            "#c5b700",
            "#cbbc00",
            "#d0c100",
            "#d6c600",
            "#dbcb03",
            "#e1d011",
            "#e6d51b",
            "#ecda24",
            "#f2df2b",
            "#f7e432",
            "#fde938",
            "#ffef40",
            "#fff647",
        ],
        "orange": [
            "#cd6511",
            "#d26916",
            "#d76d1b",
            "#dc7120",
            "#e17624",
            "#e67a28",
            "#eb7e2d",
            "#f08231",
            "#f58635",
            "#fa8a39",
            "#fe8f3e",
            "#ff9544",
            "#ff9c4b",
            "#ffa351",
            "#ffaa58",
            "#ffb15e",
            "#ffb764",
            "#ffbd6a",
            "#ffc470",
            "#ffca76",
            "#ffd07c",
            "#ffd682",
        ],
        "blue": [
            "#2d89bc",
            "#348dc0",
            "#3b92c5",
            "#4197ca",
            "#479bcf",
            "#4da0d4",
            "#52a5d9",
            "#58aade",
            "#5daee3",
            "#63b3e8",
            "#68b8ee",
            "#6ebdf3",
            "#73c2f8",
            "#78c7fd",
            "#7eccff",
            "#84d1ff",
            "#8ad7ff",
            "#8fdcff",
            "#95e2ff",
            "#9be7ff",
            "#a1edff",
            "#a6f3ff",
            "#acf8ff",
            "#b2feff",
        ],
        "green": [
            "#24b161",
            "#2bb565",
            "#32ba69",
            "#38bf6d",
            "#3ec371",
            "#44c876",
            "#4acd7a",
            "#4fd27e",
            "#54d683",
            "#5adb87",
            "#5fe08b",
            "#64e590",
            "#69ea94",
            "#6eee99",
            "#73f39d",
            "#78f8a2",
            "#7efda7",
            "#91ffb9",
        ],
    }


def plot_phylo_tree(alignment_file, save_pdf=True, path="", height=10, wideness=8):
    """Plotting a phylogenetic tree.

    Parameters
    ----------
    alignment_file : Bio.Align.MultipleSeqAlignment
        A multiple sequence alignemt made into a Bio.Align.MultipleSeqAlignment obejct
        Dataframe with categories in seperate column from x and y
    save_pdf : bool
    path : str
        path to folder with name of the file

    Returns
    -------
    plot_phylo_tree"""

    # calculate distances in protein identety:
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment_file)

    # Construct tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Plot tree
    # set the size of the figure
    plt.rc("font", family="Helvetica")
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    fig = plt.figure(figsize=(height, wideness), dpi=300)

    plt.rcParams.update({"font.size": 10})
    axes = fig.add_subplot(1, 1, 1)

    fig1 = plt.gcf()
    Phylo.draw(tree, axes=axes, do_show=False, branch_labels=None)
    pylab.axis("off")
    # homologs
    plt.rcParams.update({"font.size": 10})

    # Saving the plot
    if save_pdf and path != "":
        ## save pdf
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def plot_stacked_barplot_with_labels(
    df: pd.DataFrame,
    colors: list,
    title="",
    y_label="",
    x_label="",
    path="",
    size_length: int = 20,
    size_heigth: int = 10,
):
    """Plots a stacked barplot from a dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the data to be plotted.
    colors : list
        A list of colors for the different bars in the plot.
    title : str, optional
        The title of the plot. Default is an empty string.
    y_label : str, optional
        The label for the y-axis of the plot. Default is an empty string.
    x_label : str, optional
        The label for the x-axis of the plot. Default is an empty string.
    path : str, optional
        The path to the directory where the plot will be saved. Default is an empty string.
    """
    # Initialize
    plt.rc("font", family="Helvetica")
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    ax = df.plot(
        kind="bar",
        stacked=True,
        figsize=(size_length, size_heigth),
        color=colors,
        edgecolor="Black",
    )  # , cmap="coolwarm"
    # Add Title and Labels
    plt.title(title, fontsize=40)
    plt.xlabel(y_label, fontsize=25, weight="bold")
    plt.ylabel(x_label, fontsize=25, weight="bold")
    ax.tick_params(axis="both", which="major", labelsize=25)

    # removes the borders around the plot
    sns.despine(bottom=True, left=True)
    ax.legend([], [], frameon=False)  # around the legend

    # adding laves to each box
    for c in ax.containers:
        # this one writes label and percent
        labels_for_bars = [f"{c.get_label()} \n{round(v.get_height(),2)} %" for v in c]

        # remove the labels parameter if it's not needed for customized labels
        ax.bar_label(
            c,
            labels=labels_for_bars,
            label_type="center",
            fmt="str",
            size=20,
            weight="bold",
        )

    if path != "":
        name = "Occurences of each part sampled"
        plt.savefig(path + name + ".pdf", format="pdf", dpi=300)


def grouped_bar_plot(
    x: list,
    y: list,
    colors: list,
    category_labels: list,
    title="",
    y_label="",
    x_label="",
    path="",
    axhline=True,
    size_height: int = 10,
    size_length: int = 10,
):
    """
    Create a grouped bar plot from input data and save the plot in a pdf format

    Parameters
    ----------
    x : list
        A list of x-coordinates of the bars in the plot.
    y : list
        A list of y-coordinates of the bars in the plot.
    colors : list
        A list of color codes for the bars in the plot.
    category_labels : list
        A list of labels for the categories represented by the bars in the plot.
    title : str, optional
        The title of the plot. Default is an empty string.
    y_label : str, optional
        The label for the y-axis of the plot. Default is an empty string.
    x_label : str, optional
        The label for the x-axis of the plot. Default is an empty string.
    path : str, optional
        The path to the directory where the plot will be saved. Default is an empty string.
    axhline : bool, optional
        The flag to show an horisontal line
    """

    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # Plot
    plt.bar(x, y, edgecolor="black", color=colors)  # white - orange = #fee6ce

    # Add labels and titel
    ax.set_ylabel(y_label, size=20, fontname="Helvetica")
    ax.set_xlabel(x_label, size=30, fontname="Helvetica")
    ax.set_title(title, size=30, fontname="Helvetica")

    if axhline:
        # add horisontal line
        plt.axhline(y=100, color="black", linestyle="dashed")

    # Set color
    ax.set_facecolor("white")

    white_patch = mpatches.Patch(color="white", label=category_labels[0])
    black_patch = mpatches.Patch(color="black", label=category_labels[1])
    ax.legend(handles=[white_patch, black_patch], fontsize=20, frameon=False)

    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # SIze matters
    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_height, size_length)

    if path != "":
        name = "grouped_bar_plot"
        plt.savefig(path + name + ".pdf", format="pdf", dpi=300)

    plt.show()
