#!/usr/bin/env python

# Test csv_database
import pandas as pd
import numpy as np

# Import the csv_database modules
from teemi.lims.csv_database import *

def test_get_unique_id(): 
    unique_id = get_unique_id(path = '../teemi/tests/files_for_testing/csv_database_tests')
    assert unique_id == 10041


def test_get_database(): 
    ds_dna_box = get_database('ds_dna_box', path= '../teemi/tests/files_for_testing/csv_database_tests/')
    assert ds_dna_box.iloc[0]['ID'] == 10011.0


def test_get_plate(): 
    plasmid_plates = get_database('plasmid_plates', path= '../teemi/tests/files_for_testing/csv_database_tests/')
    box1 = get_plate(0,plasmid_plates )
    assert box1.iloc[0]['ID'] == 10000
    

def test_get_box(): 
    ds_dna_box = get_database('ds_dna_box', path= '../teemi/tests/files_for_testing/csv_database_tests/')
    box1 = get_box(0,ds_dna_box )
    assert box1.iloc[0]['name'] == 'UP_XI-2'


def test_add_unique_ids():
    test_templates = []
    for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    
    add_unique_ids(test_templates, path= '../teemi/tests/files_for_testing/csv_database_tests')

    assert test_templates[0].id == str(10041)
    assert test_templates[1].id == str(10042)
    assert test_templates[2].id == str(10043)
    assert test_templates[3].id == str(10044)


def test_add_annotations(): 
    test_templates = []
    for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    test_templates

    add_annotations(test_templates, concentration = 100)
    assert test_templates[0].annotations['batches'][0]['concentration'] == 100


def test_get_dna_from_plate_name():
    my_dna = get_dna_from_plate_name('pRS416.gb','plasmid_plates',database_path =  '../teemi/tests/files_for_testing/csv_database_tests/')
    assert len(my_dna.seq) == 4898


def test_get_dna_from_box_name():
    my_dna = get_dna_from_box_name('AanCPR_tCYC1','ds_dna_box', database_path= '../teemi/tests/files_for_testing/csv_database_tests/')
    assert len(my_dna.seq) == 2298


def test_change_row(): 
    #databse
    plasmid_plates = get_database('plasmid_plates', path= '../teemi/tests/files_for_testing/csv_database_tests/')

    # new insert
    test_templates = []
    for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    # adding annotations    
    biopython_object = [test_templates[0]]
    biopython_object = add_annotations(biopython_object )[0]
    biopython_object.id = 999999
    # changing row
    change_row(0, plasmid_plates,  biopython_object)
    print(biopython_object)
    print(plasmid_plates)
    assert plasmid_plates.iloc[0]['ID'] == 999999


def test_delete_row_df(): 
    #database
    plasmid_plates = get_database('plasmid_plates', path= '../teemi/tests/files_for_testing/csv_database_tests/')
    delete_row_df(0, plasmid_plates)
    
    assert str(plasmid_plates.iloc[0]['ID']) == 'nan'


def test_add_sequences_to_dataframe(): 
    test_templates = []

    for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    add_annotations(test_templates)
    add_unique_ids(test_templates, path= '../teemi/tests/files_for_testing/csv_database_tests')
    plasmid_plates = get_database('plasmid_plates', path= '../teemi/tests/files_for_testing/csv_database_tests/')

    add_sequences_to_dataframe(test_templates,plasmid_plates, index= 20 )

    assert str(plasmid_plates.iloc[20]['name']) == 'pCCW12'
    assert str(plasmid_plates.iloc[21]['name']) == 'pTPI1'
    assert str(plasmid_plates.iloc[22]['name']) == 'pCYC1'
    assert str(plasmid_plates.iloc[23]['name']) == 'pENO2'


def test_update_database(): 
    # making a test to test if pandas writes a csv correctly
    pass


# ---------------------------------------------------------------------------
# Additional coverage tests (self-contained, file writes go to tmp_path)
# ---------------------------------------------------------------------------
import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from teemi.lims import csv_database as csv_db

_CSV_DB_DIR = (
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "files_for_testing", "csv_database_tests")
    + os.sep
)


def _make_annotated_records(ids, names, seqs):
    records = [
        SeqRecord(Seq(seq), id=str(id_), name=name, description=name + " description")
        for id_, name, seq in zip(ids, names, seqs)
    ]
    return csv_db.add_annotations(
        records,
        concentration=42.5,
        reference="test-ref",
        volume=12,
        comments="added by test",
        location="freezer_1",
    )


def test_get_unique_id_returns_10000_when_database_has_no_ids(tmp_path):
    pd.DataFrame({"ID": [np.nan, np.nan], "name": [np.nan, np.nan]}).to_csv(
        tmp_path / "empty_plates.csv", index=False
    )
    # non-csv files in the database folder must be ignored
    (tmp_path / "notes.txt").write_text("ID\n99999\n")

    assert csv_db.get_unique_id(path=str(tmp_path)) == 10000


def test_get_unique_id_ignores_non_csv_files(tmp_path):
    pd.DataFrame({"ID": [10003.0, np.nan, 10007.0]}).to_csv(tmp_path / "db.csv", index=False)
    (tmp_path / "notes.txt").write_text("ID\n99999\n")

    assert csv_db.get_unique_id(path=str(tmp_path)) == 10008


def test_add_sequences_to_dataframe_fills_first_blank_rows():
    plasmid_plates = csv_db.get_database("plasmid_plates", path=_CSV_DB_DIR)
    first_blank = plasmid_plates.index[plasmid_plates["ID"].isna()][0]
    assert first_blank == 11
    records = _make_annotated_records(
        ids=[20001, 20002], names=["new_part_1", "new_part_2"], seqs=["ATGC", "GGGCCCAA"]
    )

    csv_db.add_sequences_to_dataframe(records, plasmid_plates)  # index=0 -> auto placement

    first, second = plasmid_plates.loc[11], plasmid_plates.loc[12]
    assert first["ID"] == 20001
    assert first["name"] == "new_part_1"
    assert first["description"] == "new_part_1 description"
    assert first["seq"] == "ATGC"
    assert first["size"] == 4
    assert first["concentration"] == 42.5
    assert first["volume"] == 12.0
    assert first["location"] == "freezer_1"
    assert first["comments"] == "added by test"
    assert first["reference"] == "test-ref"
    # plate coordinates of the row are untouched
    assert (first["plate"], first["row"], first["col"]) == (0, "A", 12)
    assert second["ID"] == 20002
    assert second["seq"] == "GGGCCCAA"
    assert second["size"] == 8
    assert (second["row"], second["col"]) == ("B", 1)
    # existing entries are not overwritten and the next row is still blank
    assert plasmid_plates.loc[10, "ID"] == 10010
    assert np.isnan(plasmid_plates.loc[13, "ID"])


def test_update_database_writes_csv_that_round_trips(tmp_path):
    plasmid_plates = csv_db.get_database("plasmid_plates", path=_CSV_DB_DIR)
    records = _make_annotated_records(ids=[30000], names=["round_trip"], seqs=["ATGAAATAG"])
    csv_db.add_sequences_to_dataframe(records, plasmid_plates, index=40)

    csv_db.update_database(plasmid_plates, "plasmid_plates_copy", path=str(tmp_path) + os.sep)

    assert sorted(os.listdir(tmp_path)) == ["plasmid_plates_copy.csv"]
    reloaded = csv_db.get_database("plasmid_plates_copy", path=str(tmp_path) + os.sep)
    assert list(reloaded.columns) == list(plasmid_plates.columns)
    assert len(reloaded) == len(plasmid_plates)
    assert reloaded.loc[40, "ID"] == 30000
    assert reloaded.loc[40, "name"] == "round_trip"
    assert reloaded.loc[40, "seq"] == "ATGAAATAG"
    assert reloaded.loc[0, "ID"] == 10000
    # the new ID is picked up by the unique id generator
    assert csv_db.get_unique_id(path=str(tmp_path)) == 30001


def _write_genbank(path, seq, name):
    record = SeqRecord(Seq(seq), id=name, name=name, description="test genbank record")
    record.annotations["molecule_type"] = "DNA"
    SeqIO.write(record, str(path), "genbank")


def test_get_dna_from_plate_name_genbank(tmp_path):
    pd.DataFrame(
        {
            "ID": [10500.0, 10501.0],
            "name": ["other_plasmid", "my_plasmid"],
            "plate": [1, 1],
            "row": ["B", "B"],
            "col": [2, 3],
            "seq": ["AAAA", "CCCC"],
            "concentration": [10.0, 55.5],
            "volume": [5.0, 20.0],
            "location": ["fridge", "freezer"],
            "description": ["other", "mine"],
        }
    ).to_csv(tmp_path / "plates.csv", index=False)
    genbank_dir = tmp_path / "genbank"
    genbank_dir.mkdir()
    _write_genbank(genbank_dir / "10501.gb", "ATGCATGCATGC", "my_plasmid")

    record = csv_db.get_dna_from_plate_name(
        "my_plasmid",
        "plates",
        database_path=str(tmp_path) + os.sep,
        genbank_files_path=str(genbank_dir) + os.sep,
        genbank=True,
    )

    # the sequence comes from the genbank file, not from the csv "seq" column
    assert str(record.seq) == "ATGCATGCATGC"
    assert record.name == "my_plasmid"
    assert record.annotations == {
        "plate": 1,
        "row": "B",
        "col": 3,
        "batches": [{"location": "freezer_1_B3", "volume": 20.0, "concentration": 55.5}],
    }


def test_get_dna_from_box_name_genbank(tmp_path):
    pd.DataFrame(
        {
            "ID": [10600.0, 10601.0],
            "name": ["my_fragment", "other_fragment"],
            "box": [2, 2],
            "row": ["C", "C"],
            "col": [7, 8],
            "seq": ["AAAA", "CCCC"],
            "concentration": [80.0, 1.0],
            "volume": [15.0, 1.0],
            "location": ["ds_dna_box_rack", "elsewhere"],
            "description": ["mine", "other"],
        }
    ).to_csv(tmp_path / "boxes.csv", index=False)
    _write_genbank(tmp_path / "10600.gb", "GGGTTTAAACCC", "my_fragment")

    record = csv_db.get_dna_from_box_name(
        "my_fragment",
        "boxes",
        database_path=str(tmp_path) + os.sep,
        genbank_files_path=str(tmp_path) + os.sep,
        genbank=True,
    )

    assert str(record.seq) == "GGGTTTAAACCC"
    assert record.name == "my_fragment"
    assert record.annotations == {
        "box": 2,
        "row": "C",
        "col": 7,
        "batches": [
            {"location": "ds_dna_box_rack_2_C7", "volume": 15.0, "concentration": 80.0}
        ],
    }


def test_get_database_missing_file_prints_hint_and_returns_none(tmp_path, capsys):
    result = csv_db.get_database("does_not_exist", path=str(tmp_path) + os.sep)

    assert result is None
    assert "Couldnt find that databse" in capsys.readouterr().out
