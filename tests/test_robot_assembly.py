#!/usr/bin/env python

# Test RobotAssembly
import sys
import pytest
import pandas as pd

from teemi.build.robot_assembly import RobotAssembly, well_keys_96, make_virtual_plates_fromDF, picklist_from_plates

# Read in test data
test_data = pd.read_excel('../teemi/tests/files_for_testing/Random_PCR_list.xlsx')

# Add test primers and templates
fwd_primers = ['F0069','F0078', 'F0083','F0088','F0093','F0098','F00103','F00108']
rev_primers = ['R00174','R00175','R00176','R00177']
templates = ['pRPL15B']

# Initiate the Robotassembly class
test_Robotassemly = RobotAssembly(test_data,fwd_primers ,rev_primers, templates)


def test_Robotassemly_correct_input():
    assert type(test_Robotassemly.pandas_PCR) == pd.core.frame.DataFrame
    assert type(test_Robotassemly.forward_primers) == list
    assert type(test_Robotassemly.reverse_primers )  == list
    assert type(test_Robotassemly.templates ) == list


def test_Robotassemly_primers():
    assert test_Robotassemly.forward_primers ==['F0069', 'F0078', 'F0083', 'F0088', 'F0093', 'F0098', 'F00103', 'F00108']
    assert test_Robotassemly.reverse_primers == ['R00174', 'R00175', 'R00176', 'R00177']
    assert test_Robotassemly.templates == ['pRPL15B']



# Helper functions
def test_well_keys96(): 

    well_keys_row = well_keys_96()

    assert type(well_keys_row) == list
    assert well_keys_row == ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12']

    well_keys_column = well_keys_96(row = False )

    assert type(well_keys_column) == list
    assert well_keys_column == ['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3', 'A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4', 'A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5', 'A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6', 'A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7', 'A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8', 'A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9', 'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10', 'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11', 'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12']


def test_make_virtual_plates_fromDF(): 
    f, r, t, pcr_mix = make_virtual_plates_fromDF(fwd_primers, rev_primers, templates, test_data)

    assert f.index_to_wellname(1) == 'A1'
    assert r.index_to_wellname(1) == 'A1'
    assert f.num_wells == 96
    assert t.num_wells == 96
    assert f.columns == {1: ['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1'], 2: ['A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2'], 3: ['A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3'], 4: ['A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4'], 5: ['A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5'], 6: ['A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6'], 7: ['A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7'], 8: ['A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8'], 9: ['A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9'], 10: ['A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10'], 11: ['A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11'], 12: ['A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12']}
    assert f.name == '1'
    assert pcr_mix.name == '4'
    assert f.wells["A1"].volume == 5.0
    assert f.wells["A1"].is_empty == False



def test_picklist_from_plates():
    f, r, t, pcr_mix = make_virtual_plates_fromDF(fwd_primers, rev_primers, templates, test_data)

    picklist = picklist_from_plates(f, r, t, pcr_mix, test_data)

    assert picklist.to_plain_string()[:38] == 'Transfer 1.00E+00L from 1 A1 into 5 A1'
    assert picklist.total_transferred_volume() == 640
