#!/usr/bin/env python

# Test utils module

from teemi.utils import *

def test_mean():
    assert mean([1, 2, 3]) == 2.0
    # Test for a normal case
    assert mean([1, 2, 3, 4, 5]) == 3
    # Test for a case with minus
    assert mean([-1, -2, -3]) == -2.0
    # Test for a list containing only one element
    assert mean([5]) == 5
    # test floats
    assert mean([1.5, 2.5, 3.5]) == 2.5

def test_counting_occurrences():
    data = {'A': 3, 'B': 2, 'C': 1}
    result = counting_occurences(data)
    assert result == ([50.0, 33.33333333333333, 16.666666666666664], ['A', 'B', 'C'])

def test_counting_occurrences_empty_dict():
    data = {}
    result = counting_occurences(data)
    assert result == ([], [])

def test_counting_occurrences_big_numbers():
    data = {'A': 300, 'B': 200, 'C': 100}
    result = counting_occurences(data)
    assert result == ([50.0, 33.33333333333333, 16.666666666666664], ['A', 'B', 'C'])

def test_counting_occurrences_negative_numbers():
    data = {'A': -3, 'B': -2, 'C': -1}
    result = counting_occurences(data)
    assert result == ([50.0, 33.33333333333333, 16.666666666666664], ['A', 'B', 'C'])

def test_unnest_dict():
    input_dict = {'A': 1, 'B': 2, 'nested': {'C': 3, 'label': 'D'}}
    result = unnest_dict(input_dict, 'nested')
    assert result == {'A': 1, 'B': 2, 'C': 3, 'name': 'D'}


def test_unnest_dict_no_label_key():
    input_dict = {'A': 1, 'B': 2, 'nested': {'C': 3}}
    result = unnest_dict(input_dict, 'nested')
    assert result == {'A': 1, 'B': 2, 'C': 3}

def test_start_end_to_location():
    input_dict = {'start': 10, 'end': 20, 'name': 'test'}
    result = start_end_to_location(input_dict, 100)
    assert "location" in result.keys()

def test_start_end_to_location_start_greater_end():
    input_dict = {'start': 20, 'end': 10, 'name': 'test'}
    result = start_end_to_location(input_dict, 100)
    assert "location" in result.keys()


def test_split_based_on_keys():
    # Test for a normal case
    dictionary = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    key_list = ['a', 'c']
    first_dict, other_dict = split_based_on_keys(dictionary, key_list)
    assert first_dict == {'a': 1, 'c': 3}
    assert other_dict == {'b': 2, 'd': 4}
    
    # Test for a case where key_list is empty
    dictionary = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    key_list = []
    first_dict, other_dict = split_based_on_keys(dictionary, key_list)
    assert first_dict == {}
    assert other_dict == {'a': 1, 'b': 2, 'c': 3, 'd': 4}

def test_rename_dict_keys():
    # Test for a normal case
    dictionary = {'a': 1, 'b': 2, 'c': 3}
    trans_dictionary = {'a': 'A', 'b': 'B'}
    renamed_dict = rename_dict_keys(dictionary, trans_dictionary)
    assert renamed_dict == {'A': 1, 'B': 2, 'c': 3}
        
    # Test for a case where trans_dictionary is empty
    dictionary = {'a': 1, 'b': 2, 'c': 3}
    trans_dictionary = {}
    renamed_dict = rename_dict_keys(dictionary, trans_dictionary)
    assert renamed_dict == {'a': 1, 'b': 2, 'c': 3}
    
def test_location_to_start_end_strand():
    # Test for a normal case
    from Bio.SeqFeature import FeatureLocation
    dictionary = {"location": FeatureLocation(5, 10, 1)}
    transformed_dict = location_to_start_end_strand(dictionary)
    assert transformed_dict == {"start": 5, "end": 10, "strand": 1}
