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

""" This part of the lab module is used for simulating and calculating PCR reactions."""

# standard libraries
import textwrap as _textwrap
import math
import csv
import json

# Extra
import pandas as pd
from pydna._pretty import pretty_str as _pretty_str
import requests


def primer_tm_neb(primer, conc=0.5, prodcode="q5-0"):
    """Calculates a single primers melting temp from NEB.

    Parameters
    ----------
    primer1 : str
    conc : float
    prodcode : str
        find product codes on nebswebsite: https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    tm : int
        primer melting temperature

    """

    url = "https://tmapi.neb.com/tm/batch"
    seqpairs = [[primer]]

    input = {"seqpairs": seqpairs, "conc": conc, "prodcode": prodcode}
    headers = {"content-type": "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)

    r = json.loads(res.content)

    if r["success"]:
        for row in r["data"]:
            return row["tm1"]
    else:
        print("request failed")
        print(r["error"][0])


def primer_ta_neb(primer1, primer2, conc=0.5, prodcode="q5-0"):
    """Calculates primer pair melting temp TA,  from NEB.

    Parameters
    ----------
    primer1 : str
        first primer to be used for finding the optimal ta
    primer2 : str
        second primer to be used for finding the optimal ta
    conc : float
    prodcode : str
        find product codes on nebswebsite: https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    ta : int
        primer pair annealing temp

    """

    url = "https://tmapi.neb.com/tm/batch"
    seqpairs = [[primer1, primer2]]

    input = {"seqpairs": seqpairs, "conc": conc, "prodcode": prodcode}
    headers = {"content-type": "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)

    r = json.loads(res.content)

    if r["success"]:
        for row in r["data"]:
            return row["ta"]

    else:
        print("request failed")
        print(r["error"][0])


def grouper(iterable, max_diff):
    """Groups objects into distinct groups based on differences"""
    prev = None
    group = []
    for item in iterable:
        if not prev or item - prev <= max_diff:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group


def calculate_volumes(vol_p_reac=0, no_of_reactions=1, standard_reagents=[], standard_volumes=[]):
    """
    Makes a reaction scheme for PCR master mixes.

    Parameters
    ----------
    vol_p_reac : int, optional
        Volume per reaction. Default is 0.
    no_of_reactions : int, optional
        Number of reactions. Default is 1.
    standard_reagents : list
        List of standard reagents.
    standard_volumes : list
        List of volumes for standard reagents.

    Returns
    -------
    volumes_df : pd.DataFrame
        DataFrame containing volume information.

    Examples
    --------
    .. code-block:: python

        calculate_volumes(
            vol_p_reac = 10,
            no_of_reactions = 6,
            standard_reagents = ["DNA","Buffer, Cutsmart","H20","Enz, USER"],
            standard_volumes = [1,1,7,1]
        )


    .. code-block:: none

        The following reaction scheme will be made:

            vol_p_reac    vol_p_x_reac
        DNA            1.0           6.0
        Buffer, Cutsmart 1.0         6.0
        H20            7.0           42.0
        Enz, USER      1.0           6.0
        Total          10.0          60.0

    """

    standard_total_volume = sum(standard_volumes)
    volumes_p_x = [val / standard_total_volume * vol_p_reac for val in standard_volumes]
    volumes_p_x_p_y_reactions = [val * no_of_reactions for val in volumes_p_x]

    volumes_p_x_plus_total = volumes_p_x + [sum(volumes_p_x)]
    volumes_p_x_p_y_reactions_plus_total = volumes_p_x_p_y_reactions + [
        sum(volumes_p_x_p_y_reactions)
    ]
    reagents_plus_total = standard_reagents + ["Total"]

    volumes_df = pd.DataFrame(
        data={
            "vol_p_reac": volumes_p_x_plus_total,
            "vol_p_"
            + str(no_of_reactions)
            + "_reac": volumes_p_x_p_y_reactions_plus_total,
        },
        index=reagents_plus_total,
    )
    return volumes_df


def calculate_processing_speed(amplicon):
    """Determines process speed based on the
    which polymerase is used.

    Parameters
    ----------
    amplicon : pydna.amplicon

    Returns
    -------
    amplicon : pydna.amplicon
        Adds annotations to the amplicon object dependent
        on which polymerase was used

    Notes
    -----
    The amplicon needs to have the following dict incorporated:
    amplicon.annotations["polymerase"]
    """

    if "proc_speed" in amplicon.forward_primer.annotations:
        print("proc_speed already set")
        return amplicon

    # proc_speed units are seconds/kb
    if amplicon.annotations["polymerase"] == "OneTaq Hot Start":
        proc_speed = 60
    elif amplicon.annotations["polymerase"] == "Q5 Hot Start":
        proc_speed = 30
    elif amplicon.annotations["polymerase"] == "Phusion":
        proc_speed = 30

    amplicon.annotations["proc_speed"] = proc_speed

    return amplicon


def calculate_elongation_time(amplicon):
    """Determines elongation time for an amplicon
    and add the elongation time to the amplicon annotations.

    Parameters
    ----------
    amplicon : pydna.amplicon

    Returns
    -------
    Adds the elongation time to the amplicon annotations

    Notes
    -----
    The amplicon needs to have a dict called proc_speed shown as follows:
    amplicon.annotations["proc_speed"]
    This dict within the annotations can be made with the function proc_speed.

    """

    if "elongation_time" in amplicon.forward_primer.annotations:
        print("elongation_time already set")
        return amplicon

    # elongation_time units are seconds
    elongation_time = amplicon.annotations["proc_speed"] * len(amplicon) / 1000
    amplicon.annotations["elongation_time"] = math.ceil(elongation_time)

    return amplicon


def calculate_required_thermal_cyclers(
    amplicons: list, polymerase: str, elong_time_max_diff=15
):
    """Determines the number of thermalcyclers that is needed
    based on elongation time differences

    Parameters
    ----------
    amplicons : list
        of pydna.amplicon objects
    polymerase : str

    Returns
    -------
    pd.DataFrame
        dataframe of grouped amplicons

    """

    amp_names = [amplicon.name for amplicon in amplicons]
    elong_times = [amplicon.annotations["elongation_time"] for amplicon in amplicons]
    tas = [amplicon.annotations["ta " + polymerase] for amplicon in amplicons]
    order = list(range(0, len(amplicons)))

    list_of_tuples = list(zip(amp_names, tas, elong_times, order))

    list_of_tuples.sort()

    groups = dict(enumerate(grouper(elong_times, elong_time_max_diff), 1))

    list_of_lists = [list(elem) for elem in list_of_tuples]

    for gNo, gTimes in groups.items():
        # print(gNo, gTimes)
        for idx, lst in enumerate(list_of_lists):
            if lst[2] in gTimes:
                list_of_lists[idx][2] = max(gTimes)

    thermal_cyclers = pd.DataFrame(
        list_of_lists, columns=["amplicons", "tas", "elong_times", "order"]
    )
    thermal_cyclers = thermal_cyclers.sort_values(["order"])
    thermal_cyclers = (
        thermal_cyclers.groupby(["tas", "elong_times"])["amplicons"]
        .apply(", ".join)
        .reset_index()
    )

    return thermal_cyclers


def pcr_locations(amplicons: list):
    """Obtain annotation information for amplicons.

    Parameters
    ----------
    amplicons : list
        List of amplicon objects from `pydna.amplicon()`

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame with amplicon locations
    """
    # initialization
    product_loc = []
    product_names = []
    template_loc = []
    fw_loc = []
    rv_loc = []

    for i in range(0, len(amplicons)):
        product_names.append(amplicons[i].name)

        # Test if batches is present
        if (
            "batches" in amplicons[i].template.annotations.keys()
            and len(amplicons[i].template.annotations["batches"]) != 0
        ):
            product_loc.append(
                amplicons[i].template.annotations["batches"][0]["location"]
            )
            template_loc.append(
                amplicons[i].template.annotations["batches"][0]["location"]
            )
        elif (
            "batches" in amplicons[i].annotations.keys()
            and len(amplicons[i].annotations["batches"]) != 0
        ):
            product_loc.append(amplicons[i].annotations["batches"][0]["location"])
            template_loc.append(amplicons[i].annotations["batches"][0]["location"])
        else:
            product_loc.append("Empty")
            template_loc.append("Empty")
            print(
                "No batches were found for "
                + str(amplicons[i].name)
                + ". Please check the object."
            )

        # Save primer locations
        if (
            "batches" in amplicons[i].forward_primer.annotations.keys()
            and len(amplicons[i].forward_primer.annotations["batches"]) != 0
        ):
            fw_loc.append(
                amplicons[i].forward_primer.annotations["batches"][0]["location"]
            )
        else:
            fw_loc.append("Empty")
            print(str(amplicons[i].name) + ": Foward primer location was not found")
        if (
            "batches" in amplicons[i].reverse_primer.annotations.keys()
            and len(amplicons[i].reverse_primer.annotations["batches"]) != 0
        ):
            rv_loc.append(
                amplicons[i].reverse_primer.annotations["batches"][0]["location"]
            )
        else:
            rv_loc.append("Empty")
            print(str(amplicons[i].name) + ": Reverse primer location was not found")

    # Save information as dataframe
    df_pcr = pd.DataFrame(
        list(zip(product_loc, product_names, template_loc, fw_loc, rv_loc)),
        columns=["location", "name", "template", "fw", "rv"],
    )

    return df_pcr


def nanophotometer_concentrations(path=""):
    """Reads a CSV file with nanophotometer concentraions
    and returns the concentrations in a list.

    Parameters
    ----------
    path : str
        path to file

    Returns
    -------
    concentrations : list
        list of concentrations from the file as floats
    """
    concentrations = []
    with open(path, encoding="Latin1") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        next(reader)[4]
        for row in reader:
            conc = float(row[4].replace(",", "."))
            concentrations.append(conc)

    return concentrations


def amplicon_by_name(name: str, amplicons_lst: list):
    """Returns amplicon with specified name.

    Parameters
    ----------
    name : str
    amplicons_lst : list

    Returns
    -------
    amplicon : pydna.amplicon

    """
    for amplicon in amplicons_lst:
        if amplicon.name == name:
            return amplicon


def Q5_NEB_PCR_program(amplicon):
    """Simple PCR program designed to give a quick visual representation
    of a PCR reaction.

    Parameters
    ----------
    amplicon : pydna.amplicon
        pydna amplicon object

    Returns
    -------
    str
        schematic representation of a Q5 program
    """
    # Determine elongation time and process speed.
    amplicon = calculate_elongation_time(amplicon)
    amplicon = calculate_processing_speed(amplicon)

    # ta
    amplicon.annotations["ta Q5 Hot Start"] = primer_ta_neb(
        str(amplicon.forward_primer.seq), str(amplicon.reverse_primer.seq)
    )

    # tm forward and reverse
    amplicon.forward_primer.annotations["tm Q5 Hot Start"] = primer_tm_neb(
        str(amplicon.forward_primer.seq)
    )
    amplicon.reverse_primer.annotations["tm Q5 Hot Start"] = primer_tm_neb(
        str(amplicon.reverse_primer.seq)
    )

    r"""Returns a string containing a text representation of a suggested
    PCR program using Taq or similar polymerase.
    ::
     |98°C|98°C               |    |tmf:59.5
     |____|_____          72°C|72°C|tmr:59.7
     |30s |10s  \ 59.1°C _____|____|30s/kb
     |    |      \______/ 0:32|5min|GC 51%
     |    |       30s         |    |1051bp
    """

    formated = _textwrap.dedent(
        r"""
                            |98°C|98°C               |    |tmf:{tmf:.1f}
                            |____|_____          72°C|72°C|tmr:{tmr:.1f}
                            |30 s|10s  \ {ta:.1f}°C _____|____|{rate}s/kb
                            |    |      \______/{0:2}:{1:2}|2min|GC {GC_prod}%
                            |    |       20s         |    |{size}bp
                            """[
            1:-1
        ].format(
            rate=amplicon.annotations["proc_speed"],
            size=len(amplicon.seq),
            ta=amplicon.annotations["ta Q5 Hot Start"],
            tmf=amplicon.forward_primer.annotations["tm Q5 Hot Start"],
            tmr=amplicon.reverse_primer.annotations["tm Q5 Hot Start"],
            GC_prod=round(amplicon.gc() * 100, 2),
            *map(int, divmod(amplicon.annotations["elongation_time"], 60)),
        )
    )

    return _pretty_str(formated)


def set_plate_locations(amplicons: list):
    """Makes a dataframe from amplicons.

    Parameters
    ----------

    amplicons : list
        list of pydna.amplicon objects

    Returns
    -------
    pd.DataFrame
        with overview of plate locations"""

    plate_locations = []
    for amplicon in amplicons:
        plate_locations.append(
            [
                amplicon.name,
                amplicon.annotations["batches"][0]["location"],
                amplicon.annotations["template_name"],
                amplicon.template.annotations["batches"][0]["location"],
                amplicon.forward_primer.id,
                amplicon.forward_primer.annotations["batches"][0]["location"],
                amplicon.reverse_primer.id,
                amplicon.reverse_primer.annotations["batches"][0]["location"],
            ]
        )
    amplicon_df = pd.DataFrame(
        plate_locations,
        columns=[
            "name",
            "location",
            "template_name",
            "template_location",
            "fw_name",
            "fw_location",
            "rv_name",
            "rv_location",
        ],
    )
    amplicon_df = amplicon_df.set_index("name")

    return amplicon_df


def update_amplicon_annotations(
    amplicon_names: list,
    amplicons: list,
    locations: list,
    concentrations: list,
    volumes: list,
) -> None:
    """Updates the annotations of amplicons in the amplicon list.

    Parameters
    ----------
    amplicon_names : list
        List of amplicon names.
    locations : list
        List of locations for each amplicon.
    concentrations : list
        List of concentrations for each amplicon.
    volumes : list
        List of volumes for each amplicon.

    Returns
    -------
    None
    """
    for i in range(len(amplicon_names)):
        amplicon_by_name(amplicon_names[i], amplicons).annotations["batches"][0][
            "location"
        ] = locations[i]
        amplicon_by_name(amplicon_names[i], amplicons).annotations["batches"][0][
            "concentration"
        ] = concentrations[i]
        amplicon_by_name(amplicon_names[i], amplicons).annotations["batches"][0][
            "volume"
        ] = volumes[i]


def get_amplicons_by_row(row, amplicon_df, amplicons):
   """Returns a list of amplicons in a given gel row.

   Parameters
   ----------
   row : str
       Name of the gel row.
   amplicon_df : pandas DataFrame
       DataFrame with amplicon information, including the column 'prow' indicating the gel row.
   amplicons : list of Amplicon
       List of Amplicon objects.

   Returns
   -------
   list of Amplicon
       List of Amplicon objects in the given gel row.
   """
   row_names = amplicon_df[amplicon_df['prow']==row][['name']]['name'].tolist()

   row_amplicons = []
   for name in row_names:
       for amplicon in amplicons:
           if amplicon.name == name:
               row_amplicons.append([amplicon])

   return(row_amplicons)


def get_amplicons_by_column(col, amplicon_df, amplicons):
   """
   Returns a list of amplicons in a given gel column.

   Parameters
   ----------
   col : str
       Name of the gel column.
   amplicon_df : pandas DataFrame
       DataFrame with amplicon information, including the column 'pcol' indicating the gel column.
   amplicons : list of Amplicon
       List of Amplicon objects.

   Returns
   -------
   list of Amplicon
       List of Amplicon objects in the given gel column.
   """
   col_names = amplicon_df[amplicon_df['pcol']==col][['name']]['name'].tolist()

   col_amplicons = []
   for name in col_names:
       for amplicon in amplicons:
           if amplicon.name == name:
               col_amplicons.append([amplicon])

   return(col_amplicons)
