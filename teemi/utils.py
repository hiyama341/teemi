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

""" A script that provides utility functions"""

# Typing
from typing import Dict, Any, List, Tuple
import Bio


def mean(lst: list):
    """Get mean from a list.

    Parameters
    ----------
    lst : list of floats.

    Return
    ------
    mean : float
    """
    if not lst:
        raise ValueError("List cannot be empty.")
    mean = sum(lst) / len(lst)
    return mean


def counting_occurences(data_with_occurences: dict):
    """
    Count the occurences of each key in the input dict and returns the percentage of each key in the total values

    Parameters
    ----------
    data_with_occurences : dict
        The dictionary containing the data for counting occurences.
    Returns
    -------
    tuple
        A tuple containing two lists, the first one is the occurence percentage of each key, the second one is the list of keys.
    """
    columns = []
    data = []
    for key, value in data_with_occurences.items():

        values = data_with_occurences.values()
        total = sum(values)

        data.append((value / total) * 100)
        columns.append(key)

    return data, columns


def unnest_dict(dictionary: Dict[Any, Any], key_unnest_dict: str) -> Dict[Any, Any]:
    """Unnest a dictionary by merging the values of a nested dictionary
    with the original dictionary, and renaming the key "label" to "name"

    Parameters:
    -----------
    dictionary : Dict[Any, Any]
        The input dictionary
    key_unnest_dict : str
        The key of the nested dictionary to unnest

    Returns:
    --------
    dictionary : Dict[Any, Any]
        The unnested dictionary
    """
    nested_dict = dictionary.pop(key_unnest_dict)

    unnested_dictionary = {**dictionary, **nested_dict}

    if "label" in unnested_dictionary.keys():
        unnested_dictionary["name"] = unnested_dictionary.pop("label")

    return unnested_dictionary


def nest_dict(
    dictionary: Dict[Any, Any],
    key_for_nested_dict: str,
    first_order_keys: List[str] = None,
) -> Dict[Any, Any]:
    """ "
    Nest a dictionary by moving the values of specified keys to a nested dictionary.

    Parameters:
    -----------
    dictionary : Dict[Any, Any]
        The input dictionary
    key_for_nested_dict : str
        The key to use for the nested dictionary
    first_order_keys : List[str], optional
        List of keys to keep in the first level of the dictionary, by default None

    Returns:
    --------
    dictionary : Dict[Any, Any]
        The nested dictionary

    """
    if first_order_keys is None:
        first_order_keys = []
    else:
        first_order_keys = first_order_keys

    first_order_dict, second_order_dict = split_based_on_keys(
        dictionary, first_order_keys
    )

    # Create qualifer dict
    dict_to_be_nested = {}
    dict_to_be_nested[key_for_nested_dict] = second_order_dict

    # Merge dictionaries by nesting qualifer dict
    nested_dictionary = {**first_order_dict, **dict_to_be_nested}

    return nested_dictionary


def start_end_to_location(dictionary: Dict[str, Any], length: int) -> Dict[str, Any]:
    """Start and End Key Value pair to Compound Location Key Value pair.

    Parameters:
    -----------
    dictionary : Dict[str, Any]
        The input dictionary containing "start" and "end" keys
    length : int
        The length of the sequence

    Returns:
    --------
    dictionary : Dict[str, Any]
        The dictionary with "location" key added and start, end removed.
    """
    start = dictionary.pop("start")
    end = dictionary.pop("end")
    start_pos = Bio.SeqFeature.ExactPosition(start)
    end_pos = Bio.SeqFeature.ExactPosition(end)
    if start_pos < end_pos:
        location = Bio.SeqFeature.FeatureLocation(start_pos, end_pos)
    else:
        f1 = Bio.SeqFeature.FeatureLocation(start_pos, length - 1)
        f2 = Bio.SeqFeature.FeatureLocation(0, end_pos)
        location = Bio.SeqFeature.CompoundLocation([f1, f2])

    dictionary["location"] = location
    return dictionary


def split_based_on_keys(
    dictionary: Dict[Any, Any], key_list: List[str]
) -> Tuple[Dict[Any, Any], Dict[Any, Any]]:
    """Split a dictionary into two based on a list of keys.
    Parameters:
    -----------
    dictionary : Dict[Any, Any]
        The input dictionary
    key_list : List[str]
        The list of keys to split the dictionary on

    Returns:
    --------
    first_dict : Dict[Any, Any]
        The dictionary containing the keys specified in key_list
    other_dict : Dict[Any, Any]
        The dictionary containing the keys not specified in key_list
    """
    # Split keys
    first_keys = key_list
    other_keys = [k for k in dictionary.keys() if k not in first_keys]

    # Create dicts
    first_dict = {k: dictionary[k] for k in first_keys}
    other_dict = {k: dictionary[k] for k in other_keys}

    return (first_dict, other_dict)


def rename_dict_keys(dictionary: Dict, trans_dictionary: Dict) -> Dict:
    """rename the keys of a dictionary using another dictionary
    Parameters:
    -----------
    dictionary : Dict[K, V]
        The input dictionary
    trans_dictionary : Dict[K, K]
        The dictionary containing the keys to be replaced as keys and the new keys as values

    Returns:
    --------
    dictionary : Dict[K, V]
        The dictionary with the keys renamed.
    """

    keys = dictionary.keys()
    values = dictionary.values()

    new_keys = [trans_dictionary.get(k, k) for k in keys]
    return dict(zip(new_keys, values))


def location_to_start_end_strand(dictionary: Dict[str, Any]) -> Dict[str, Any]:
    """Convert "location" as Bio.SeqFeature.CompoundLocation
    to "start", "end", "strand" key-value pairs.

    Parameters:
    -----------
    dictionary : Dict[str, Any]
        The input dictionary containing "location" key

    Returns:
    --------
    dictionary : Dict[str, Any]
        The dictionary with "start", "end", "strand" key added and location removed.
    """

    Compoundloc = dictionary.pop("location")

    dictionary["start"] = Compoundloc.start.real
    dictionary["end"] = Compoundloc.end.real
    dictionary["strand"] = Compoundloc.strand

    return dictionary
