#!/usr/bin/env python
# MIT License
# Copyright (c) 2022, Technical University of Denmark (DTU)
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


""" This part of the lab module is used for making transformations"""
# standard libraries
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
from pydna.dseqrecord import Dseqrecord

# for plotting with pyplot
import matplotlib.pyplot as plt

def ng_to_nmol(ng: float, bp: float):
    """Calculates nanogram to nanomol for transformation mixes.

    To do a transformation it is important to have the right ratio
    of plasmid to insert. In other words this is done by calculating
    the nanomolar ratios and this tool can do that

    Parameters
    ----------

    ng : float
        eg. nanogram
    param: float
        eg. number of basepairs. Can also be an int

    Returns
    -------
    ng_to_nmol : float

    Note
    ----
    It calculates the nmol in the following way:
    nmol = ng/(bp*650)
    """
    if ng > 0 and bp > 0:
        return ng / (bp * 650)
    else:
        return "non-valid_input"


def ODtime(initialOD: float, time: float, td: float = 0.5):
    """Calculates the OD based on doupling time.

    Parameters
    ----------
    initialOD : float
        in OD
    time : float
        in hours
    td : float
        doupling time i.e. td in h^-1

    Returns
    -------
    OD : float
        the OD after a certain time()
    """
    if initialOD >= 0 and time >= 0:
        return round(initialOD * 2 ** (time * td), 3)
    else:
        return "non-valid_input"


def time_to_inculate(
    initialOD=0.0025, td=0.4, verbose=False, transformation_time: int = 12
):

    """Calculates when a starter culture is ready to be inoculated
     with a transformation mixture.

    Parameters
    ----------
    initialOD : float
    td : float
        is doubling time
    transformation_time : int
        The time you want to transform
    verbose : Bool
        Provides extra information

    Returns
    -------
    A plot of cell growth at different td

    Notes
    -----
    This is used to calculate when the cells should be used for transformation.
    For example:
    OD_1 = 1 * 10^7 cells / ml
    For a succesfull S.cerevisiae transformation between 1 to 2 × 10^7 cells/ml should be used
    Normal doupling time is between 3-6 hours


    """
    if verbose:
        print("GOAL: to get enough cells in exponential phase for transformation")
    print("Assumed that: ")
    print(
        "- transformation time "
        + str(transformation_time)
        + " (reached OD=1 the day after)"
    )

    times = list(range(0, 37))
    ods_025_3 = [ODtime(initialOD, time, td=0.3) for time in times]
    ods_025_input = [ODtime(initialOD, time, td=td) for time in times]
    ods_025_5 = [ODtime(initialOD, time, td=0.5) for time in times]
    fig = plt.figure()
    ax = plt.axes()

    # ax.set_xlim(lims)
    ax.set_ylim([0, 2.2])
    ax.plot(times, [2] * len(times), "r-", label="end of exponential phase")
    ax.plot(times, [1] * len(times), "k-", label="target")
    ax.plot(times, ods_025_3, label="iOD=" + str(initialOD) + ", td=0.3")
    ax.plot(times, ods_025_input, label="iOD=" + str(initialOD) + ", td=" + str(td))
    ax.plot(times, ods_025_5, label="iOD=" + str(initialOD) + ", td=0.5")

    plt.xlabel("time, h^-1")
    plt.ylabel("OD")
    plt.legend()
    plt.show()

    def inoculation_time(times, ods):
        def find_closest(A, target):
            # A must be sorted
            idx = A.searchsorted(target)
            idx = np.clip(idx, 1, len(A) - 1)
            left = A[idx - 1]
            right = A[idx]
            idx -= target - left < right - target
            return idx

        # In how many hours will the cells have reached OD 1?
        hours_to_OD1 = times[find_closest(np.array(ods), 1)]
        print("Hours to OD = 1: \t" + str(hours_to_OD1) + " hours")

        ### When do u need to innoculate?
        when_to_inoculate = transformation_time - hours_to_OD1

        if when_to_inoculate < 0:
            print("Transformation time has been set to ", transformation_time)
            print(
                "Time of inoculation: \t"
                + str(when_to_inoculate + 24)
                + " (the day before)"
            )

        else:
            print("Transformation time has been set to ", transformation_time)
            print(
                "Time of inoculation: \t" + str(when_to_inoculate) + "(the day before)"
            )

        # If i innoculate now?

        print(
            "\nIf you innoculate now, the cells will have reached OD= 1 by:  ",
            datetime.now() + timedelta(hours=hours_to_OD1),
        )

    inoculation_time(times, ods_025_input)

    if verbose:
        print()
        print(
            "How to hit initialOD = 0.0025 (e.g. from colony)? Guess. Inoculate 9/10 + 1/10 'normal' colony per ~10 ml"
        )
        print("How much volume? ~2 ml per transformation")


def transformation_mix(
    reaction_names, reaction_participants, wanted_amounts, water_dna_p_reac, media=""
):

    """Makes a pandas dataframe of the parts (their location)
     that needs to be put into the transformation mixes.

    Parameters
    ----------
    reaction_names : list
        list of reaction names
    reaction_participants : list
        list of pydna.Dseqrecord objects of Bio.seqrecord objects
    wanted_concentrations : dict
        dict of the names of the reactants with their calculated nmol
    water_dna_p_reac : int
        the amount of water wanted for the reaction
    media : list
        list of names of the media used. e.g. ['LB_AMP']

    Returns
    -------
    pandas.DataFrame
        with a transformation scheme showing which parts should be
        mixed for each reaction including positive and negative
        controls.

    Examples
    --------
    # 1. Mention which reacion names you have
    reaction_names = ["insert", "n.ctr", "n.ctr", "n.ctr", "p. ctr"]

    # 2. Add reaction reaction_participants
    reaction_participants = [[vector, gRNA1_pcr_prod,gRNA2_pcr_prod], #the insert we want
                     [vector],                                        #negative control
                     [gRNA1_pcr_prod],                                #negative control
                     [gRNA2_pcr_prod],                                #negative control
                     [LEU_plasmid]]                                   #positive control

    # 2. Calculate nmol:
    nmol_vector = ng_to_nmol(ng = 15, bp = len(vector))
    nmol_gRNA = ng_to_nmol(ng = 30, bp = len(gRNA1_pcr_prod))
    nmol_pctr = ng_to_nmol(ng = 10, bp = len(LEU_plasmid))

    # 3. Add the concentrations
    wanted_concentrations = {'p0056\\(pESC-LEU-ccdB-USER)' : nmol_vector,
                     'ATF1'                        : nmol_gRNA,
                     'CroCPR'                      : nmol_gRNA,
                     'LEU_plasmid'                 : nmol_pctr}

    # 4. what media the transformants are plated on (5 transformations here)
    media = ['LB_AMP'] * 5

    # 5. initate the function
    transformation_mix(reaction_names, reaction_participants, wanted_amounts =
    (wanted_concentrations, water_dna_p_reac = 7, media = media)

    Return:
                #these are freezer locations
            name	l4_I06	l4_I07	l4_I08	p1_F06	water	plate on
    0	insert	0.1	   0.6 	   0.6	    NaN	    5.7	    LB_AMP
    1	n.ctr	0.1 	NaN    	NaN    	NaN    	6.9    	LB_AMP
    2	n.ctr	NaN 	0.6    	NaN    	NaN    	6.4    	LB_AMP
    3	n.ctr	NaN 	NaN    	0.6    	NaN    	6.4    	LB_AMP
    4	p. ctr	NaN    	NaN    	NaN    	0.1    	6.9    	LB_AMP
    """

    df_comb = pd.DataFrame()

    for name, parts in zip(reaction_names, reaction_participants):
        names = [part.name for part in parts]
        locations = [part.annotations["batches"][0]["location"] for part in parts]
        concentrations = [
            part.annotations["batches"][0]["concentration"] for part in parts
        ]  # ng/ul
        part_names = [part.name for part in parts]  # ng/ul
        sizes = [len(part) for part in parts]  # in bp

        part_mass = [
            round(wanted_amounts.get(pname, "") * int(size) * 650, 1)
            for pname, size in zip(part_names, sizes)
        ]  # in ng = nmol * bp * 650 ng/(nmol * bp)

        part_volume = [
            round(mass / con, 1) for mass, con in zip(part_mass, concentrations)
        ]  # in µl

        di = dict(zip(names, part_volume))

        df = pd.DataFrame(data=di, index=[name])  # ,index  = reagents_plus_total

        df_comb = pd.concat([df_comb, df], sort=False)

    df_comb["water"] = water_dna_p_reac - df_comb.sum(axis=1)

    df_comb = df_comb.reset_index()

    df_comb = df_comb.rename(columns={"index": "name"})

    if media != "":
        df_comb["plate on"] = media

    return df_comb


def wanted_mass(wanted_moles, size):
    """Calculates the mass needed from the specified amount 
    of moles and size. 

    Parameters
    ----------
    wanted_moles : int
        wanted moles in nmol
    size : int
        size in bp

    Returns
    -------
    w_mass_rounded : int
        in ng. Mass wanted for the reaction.
    """
    w_mass = wanted_moles * size * 650
    w_mass_rounded = round(w_mass, 1)
    return w_mass_rounded


def wanted_volume(wanted_mass, actual_concentration):
    """Calculates the wanted volume from the mass and 
    concentration.

    Parameters
    ----------
    wanted_mass : int
        wanted mass in ng

    actual_concentration : int
            actual_concentration in ng/ul

    Returns
    -------
    wanted_volume_rounded : int
        return in ul
    """
    wanted_volume = wanted_mass / actual_concentration
    wanted_volume_rounded = round(wanted_volume, 1)
    return wanted_volume_rounded


def transformation_partitipants(reaction_participants, amnt =0.0005 , sgRNA_plasmid_name= None, sgRNA_plasmid_conc= None):
    """Returns a dict with the µl amounts needed in a transformation reaction.

    Parameters
    ----------
    reaction_participants : list of list of Dseqrecord
        List of lists of Dseqrecord objects representing the reaction participants.
    amnt : float, optional
        Amount in µl of the reagents other than `sgRNA_plasmid_name`. Default is 0.0005.
    sgRNA_plasmid_name : str, optional
        Name of the sgRNA plasmid. If not provided, `amnt` is used for all reaction participants.
    sgRNA_plasmid_conc : float, optional
        Concentration in µl of the sgRNA plasmid. If not provided, `amnt` is used for all reaction participants.

    Returns
    -------
    dict
        Dict with the µl amounts needed for the transformation reaction, 
        with keys being the names of the reaction participants and values being the corresponding µl amounts.
    """
    ...
    # Initialize two lists 
    wanted_amounts = [[] for i in range(len(reaction_participants))]
    names_matrix = [[] for i in range(len(reaction_participants))]
    for reac_no, reac in enumerate(reaction_participants):
        for parti_no, parti in enumerate(reac):
            if sgRNA_plasmid_name == None and sgRNA_plasmid_conc == None: 
                wanted_amounts[reac_no].append(amnt)
                names_matrix[reac_no].append(parti.name)
            else:
                names_matrix[reac_no].append(parti.name)
                if parti.name == sgRNA_plasmid_name:
                    wanted_amounts[reac_no].append(sgRNA_plasmid_conc)
                    names_matrix[reac_no].append(parti.name)
                else: 
                    wanted_amounts[reac_no].append(amnt)


    # making the reaction participants into Dseqrecords and changing names        
    new_dict_with_wanted_amounts = dict()
    for i in range(len(reaction_participants)):
        for j in range(len(reaction_participants[i])): 
            reaction_participants[i][j] = Dseqrecord(reaction_participants[i][j])
            new_dict_with_wanted_amounts[reaction_participants[i][j].name] = wanted_amounts[i][j]
    
    
    return new_dict_with_wanted_amounts 


def calculate_volume_and_total_concentration(amplicons, amplicon_parts_amounts_total, n= 1):
    """Calculates the volume and total concentration of 
    a list of DNA parts.
    
    Parameters
    ----------
    amplicons : list
        A list of amplicon objects
    amplicon_parts_amounts_total : dict
        A dictionary of amplicon names and their respective total amounts
    n : int (optional)
        Gives the option of multiplying the volume is needed. Optional set to 1. 

    Returns
    -------
    volumes : list
        List of volumes of each amplicon
    ngs : list
        List of ngs of each amplicon
    total_conc : float
        Total concentration of all amplicons
    """
    print('name, volume, concentration, location')

    volumes = []
    ngs = []
    for amp in amplicons:
        w_moles = amplicon_parts_amounts_total[amp.name]
        w_mass = wanted_mass(wanted_moles=w_moles, size=len(amp))
        act_conc = amp.annotations['batches'][0]['concentration']
        w_volume = wanted_volume(w_mass, act_conc)*n
        volumes.append(w_volume)
        ngs.append(w_volume * act_conc)
        print(amp.name, w_volume, act_conc, '\t', amp.annotations['batches'][0]['location'])

    #Count total concentrtaion expected   
    total_vol = sum(volumes)
    total_ngs = sum(ngs)
    total_conc = total_ngs/total_vol
    print('total volume: ', sum(volumes))
    print()
    print('total ngs: ', sum(ngs))
    print('total conc: ', total_conc)
    
    return volumes, ngs, total_conc
    

def pool_parts(amplicons:list, part_names:list,part_amounts:list,  pool_names:list, pool_lengths)->dict: 
    """Pools amplicon parts and returns a dictionary of pooled volumes.

    Parameters
    ----------
    amplicons : list
        List of amplicon objects.
    part_names : list
        List of part names.
    part_amounts : list
        List of amounts of each part.
    pool_names : list
        List of pool names.
    pool_lengths : list
        List of pool lengths.

    Returns
    -------
    pooled_volumes : dict
        Dictionary containing the pooled volumes for each amplicon part.
    """
    # intialize
    pooled_volumes = {}
    
    #Iterate through the parts that are avilable
    for amplicon in amplicons:
        if amplicon.template.name in part_names:
            
            # calculate volume needed 
            ind1 = part_names.index(amplicon.template.name)
            amount = part_amounts[ind1]
            
            ind2 = pool_names.index(amplicon.template.name)
            vol = (pool_lengths[ind2]*650*amount)/amplicon.annotations['batches'][0]['concentration']
            
            # add it to the dictionary
            if amplicon.template.name in pooled_volumes:
                pooled_volumes[amplicon.template.name][amplicon.name] = {'volume_to_mix':round(vol,1),'location':amplicon.annotations['batches'][0]['location'], 'concentration':amplicon.annotations['batches'][0]['concentration']}
            else:
                pooled_volumes[amplicon.template.name] = {amplicon.name: {'volume_to_mix':round(vol,1),'location':amplicon.annotations['batches'][0]['location'], 'concentration':amplicon.annotations['batches'][0]['concentration']}}

    return pooled_volumes


def print_pooled_parts(pooled_volumes: dict) -> None:
    """Prints the pooled parts and calculated concentrations.

    Parameters
    ----------
    pooled_volumes : dict
        Dictionary containing the pooled volumes for each amplicon part.

    Returns
    -------
    None
    """
    print("To be pooled together")
    con_per_part = {}
    for key in pooled_volumes:
        print(key)
        total_vol = 0
        total_con = 0
        total_ng = 0
        for ke in pooled_volumes[key]:
            print(ke, pooled_volumes[key][ke])
            total_vol += pooled_volumes[key][ke]['volume_to_mix']
            total_con += pooled_volumes[key][ke]['concentration']
            total_ng += pooled_volumes[key][ke]['concentration'] * pooled_volumes[key][ke]['volume_to_mix']
        print("vol", round(total_vol, 1))
        print("calculated con", total_ng / total_vol, '\n')
        con_per_part[key] = round(total_ng / total_vol)
        total_con = 0
        total_vol = 0
        total_ng = 0
