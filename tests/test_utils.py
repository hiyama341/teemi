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
    data = {"A": 3, "B": 2, "C": 1}
    result = counting_occurences(data)
    assert result == ([50.0, 33.33333333333333, 16.666666666666664], ["A", "B", "C"])


def test_counting_occurrences_empty_dict():
    data = {}
    result = counting_occurences(data)
    assert result == ([], [])


def test_counting_occurrences_big_numbers():
    data = {"A": 300, "B": 200, "C": 100}
    result = counting_occurences(data)
    assert result == ([50.0, 33.33333333333333, 16.666666666666664], ["A", "B", "C"])


def test_counting_occurrences_negative_numbers():
    data = {"A": -3, "B": -2, "C": -1}
    result = counting_occurences(data)
    assert result == ([50.0, 33.33333333333333, 16.666666666666664], ["A", "B", "C"])


def test_unnest_dict():
    input_dict = {"A": 1, "B": 2, "nested": {"C": 3, "label": "D"}}
    result = unnest_dict(input_dict, "nested")
    assert result == {"A": 1, "B": 2, "C": 3, "name": "D"}


def test_nest_dict():
    input_dict = {
        "key1": "value1",
        "key2": "value2",
        "key3": "value3",
        "key4": "value4",
    }
    key_for_nested_dict = "nested_dict_key"
    first_order_keys = ["key1", "key2"]
    expected_output_dict = {
        "key1": "value1",
        "key2": "value2",
        "nested_dict_key": {"key3": "value3", "key4": "value4"},
    }
    output_dict = nest_dict(input_dict, key_for_nested_dict, first_order_keys)
    assert output_dict == expected_output_dict


def test_unnest_dict_no_label_key():
    input_dict = {"A": 1, "B": 2, "nested": {"C": 3}}
    result = unnest_dict(input_dict, "nested")
    assert result == {"A": 1, "B": 2, "C": 3}


def test_start_end_to_location():
    input_dict = {"start": 10, "end": 20, "name": "test"}
    result = start_end_to_location(input_dict, 100)
    assert "location" in result.keys()


def test_start_end_to_location_start_greater_end():
    input_dict = {"start": 20, "end": 10, "name": "test"}
    result = start_end_to_location(input_dict, 100)
    assert "location" in result.keys()


def test_split_based_on_keys():
    # Test for a normal case
    dictionary = {"a": 1, "b": 2, "c": 3, "d": 4}
    key_list = ["a", "c"]
    first_dict, other_dict = split_based_on_keys(dictionary, key_list)
    assert first_dict == {"a": 1, "c": 3}
    assert other_dict == {"b": 2, "d": 4}

    # Test for a case where key_list is empty
    dictionary = {"a": 1, "b": 2, "c": 3, "d": 4}
    key_list = []
    first_dict, other_dict = split_based_on_keys(dictionary, key_list)
    assert first_dict == {}
    assert other_dict == {"a": 1, "b": 2, "c": 3, "d": 4}


def test_rename_dict_keys():
    # Test for a normal case
    dictionary = {"a": 1, "b": 2, "c": 3}
    trans_dictionary = {"a": "A", "b": "B"}
    renamed_dict = rename_dict_keys(dictionary, trans_dictionary)
    assert renamed_dict == {"A": 1, "B": 2, "c": 3}

    # Test for a case where trans_dictionary is empty
    dictionary = {"a": 1, "b": 2, "c": 3}
    trans_dictionary = {}
    renamed_dict = rename_dict_keys(dictionary, trans_dictionary)
    assert renamed_dict == {"a": 1, "b": 2, "c": 3}


def test_location_to_start_end_strand():
    # Test for a normal case
    from Bio.SeqFeature import FeatureLocation

    dictionary = {"location": FeatureLocation(5, 10, 1)}
    transformed_dict = location_to_start_end_strand(dictionary)
    assert transformed_dict == {"start": 5, "end": 10, "strand": 1}


def test_multiply_list():
    input_list = [1, 2, 3, 4]
    expected_output = 24
    output = multiply_list(input_list)
    assert output == expected_output


def test_remove_tuple_duplicates():
    input_list = [(1, 2), (3, 4), (1, 2), (5, 6), (3, 4)]
    expected_output = [(1, 2), (3, 4), (5, 6)]
    output = remove_tuple_duplicates(input_list)
    assert output == expected_output


################################
# For testing of remove_tuple_duplicates_with_name_attributes
class Record:
    def __init__(self, name, age):
        self.name = name
        self.age = age


def test_remove_duplicates_with_name_attribute():
    input_list = [
        Record("Alice", 30),
        Record("Bob", 35),
        Record("Alice", 40),
        Record("Charlie", 45),
        Record("Bob", 50),
    ]
    expected_output = [Record("Alice", 30), Record("Bob", 35), Record("Charlie", 45)]
    output = remove_duplicates_with_name_attribute(input_list)
    assert len(output) == len(expected_output)
    assert all([output[i].name == expected_output[i].name for i in range(len(output))])
    assert all([output[i].age == expected_output[i].age for i in range(len(output))])
