import os

import pytest
from teemi.build.containers_wells_picklists import Plate96, Transfer, PickList

##### From picklists

source = Plate96(name="Source")
destination = Plate96(name="Destination")
source_well = source.wells["A1"]
destination_well = destination.wells["B2"]
volume = 25 * 10 ** (-6)
transfer_1 = Transfer(source_well, destination_well, volume)
picklist = PickList()


def test_add_transfer():
    picklist.add_transfer(transfer=transfer_1)
    assert isinstance(picklist.transfers_list[0], Transfer)


def test_to_plain_string():
    assert (
        picklist.to_plain_string()
        == "Transfer 2.50E-05L from Source A1 into Destination B2"
    )


def test_to_plain_textfile(tmpdir):
    path = os.path.join(str(tmpdir), "test.txt")
    picklist.to_plain_textfile(filename=path)
    assert os.path.exists(path)


def test_simulate():
    with pytest.raises(ValueError):
        picklist.simulate(inplace=False)


def test_restricted_to():
    new_picklist = picklist.restricted_to(
        source_well=destination_well, destination_well=destination_well
    )
    assert len(new_picklist.transfers_list) == 0

    new_picklist_2 = picklist.restricted_to(
        source_well=source_well, destination_well=destination_well
    )
    assert len(new_picklist_2.transfers_list) == 0


def test_sorted_by():
    assert isinstance(PickList().sorted_by(), PickList)


def test_total_transferred_volume():
    assert picklist.total_transferred_volume() == 25 * 10 ** (-6)


def test_enforce_maximum_dispense_volume():
    new_picklist = picklist.enforce_maximum_dispense_volume(5 * 10 ** (-6))
    assert len(new_picklist.transfers_list) == 5


def test_merge_picklists():
    new_picklist = picklist.merge_picklists([picklist, picklist])
    assert len(new_picklist.transfers_list) == 2


### Fron test_plate

from teemi.build.containers_wells_picklists import Plate96, Transfer, PickList, Well


def condition(well):
    return well.volume > 20 * 10 ** (-6)


def test_find_unique_well_by_condition():
    with pytest.raises(Exception):
        Plate96().find_unique_well_by_condition(condition)


def test_find_unique_well_containing():
    with pytest.raises(Exception):
        Plate96().find_unique_well_containing("testquery")


def test_list_well_data_fields():
    with pytest.raises(KeyError):
        Plate96().list_well_data_fields()


def test_return_column():
    assert isinstance(Plate96().return_column(5)[0], Well)
    assert len(Plate96().return_column(5)) == 8


def test_list_wells_in_column():
    assert isinstance(Plate96().list_wells_in_column(5)[0], Well)


def test_return_row():
    assert isinstance(Plate96().return_row("A")[0], Well)
    assert isinstance(Plate96().return_row(1)[0], Well)
    assert len(Plate96().return_row("A")) == 12


def test_list_wells_in_row():
    assert isinstance(Plate96().list_wells_in_row(5)[0], Well)


def test_list_filtered_wells():
    def condition(well):
        return well.volume > 50

    assert Plate96().list_filtered_wells(condition) == []


def test_wells_grouped_by():
    assert len(Plate96().wells_grouped_by()[0][1]) == 96


def test_get_well_at_index():
    well = Plate96().get_well_at_index(5)
    assert well.name == "A5"


wellname_data = [
    ("A5", "row", 5),
    ("A5", "column", 33),
    ("C6", "row", 30),
    ("C6", "column", 43),
]
inverted_wellname_data = [[s[-1], s[1], s[0]] for s in wellname_data]


@pytest.mark.parametrize("wellname, direction, expected", wellname_data)
def test_wellname_to_index(wellname, direction, expected):
    assert Plate96().wellname_to_index(wellname, direction) == expected


@pytest.mark.parametrize("index, direction, expected", inverted_wellname_data)
def test_index_to_wellname(index, direction, expected):
    assert Plate96().index_to_wellname(index, direction) == expected


def test_iter_wells():
    result = Plate96().iter_wells()
    assert isinstance(next(result), Well)


def test___repr__():
    assert Plate96().__repr__() == "Plate96(None)"


# From test_well
from teemi.build.containers_wells_picklists import TransferError


plate = Plate96()
well = plate.get_well_at_index(1)


def test_volume():
    assert well.volume == 0


def test_iterate_sources_tree():
    result = well.iterate_sources_tree()
    assert isinstance(next(result), Well)


def test_add_content():
    plate = Plate96()
    well = plate.get_well_at_index(1)
    components_quantities = {"Compound_1": 5}
    volume = 20 * 10 ** (-6)  # 20 uL
    well.add_content(components_quantities, volume=volume)
    assert well.content.quantities == {"Compound_1": 5}

    well2 = plate.get_well_at_index(2)
    well2.add_content(components_quantities, volume=20, unit_volume="uL")
    assert well2.content.concentration() == 250000.00000000003


def test_subtract_content():
    components_quantities = {"Compound_1": 5}
    volume = 30 * 10 ** (-6)  # 30 uL
    with pytest.raises(TransferError):
        well.subtract_content(components_quantities, volume)


def test_empty_completely():
    well.empty_completely()
    assert well.content.volume == 0


def test___repr__():
    assert well.__repr__() == "(None-A1)"


def test_pretty_summary():
    result = well.pretty_summary()
    expected = "(None-A1)\n  Volume: 0\n  Content: \n  Metadata: "
    assert result == expected


def test_to_dict():
    result = well.to_dict()
    expected = {
        "name": "A1",
        "content": {"volume": 0, "quantities": {}},
        "row": 1,
        "column": 1,
    }
    assert result == expected


def test_index_in_plate():
    result = well.index_in_plate()
    expected = 1
    assert result == expected


other_well = plate.get_well_at_index(2)


def test_is_after():
    assert well.is_after(other_well) is False
    assert other_well.is_after(well) is True


def test___lt__():
    assert True


# From wellcontent
from teemi.build.containers_wells_picklists import WellContent


wellcontent = WellContent(
    quantities={"Compound_1": 5, "Compound_2": 10}, volume=25
)  # 30 L [sic]


def test_concentration():
    assert WellContent().concentration() == 0
    assert WellContent(quantities={"Compound_1": 5}).concentration() == 0

    assert wellcontent.concentration() == 0.2
    assert wellcontent.concentration("Compound_1") == 0.2
    assert wellcontent.concentration("Compound_2") == 0.4
    assert wellcontent.concentration("Compound_3") == 0  # not in wellcontent


def test_to_dict():
    result = wellcontent.to_dict()
    expected = {"volume": 25, "quantities": {"Compound_1": 5, "Compound_2": 10}}
    assert result == expected


def test_make_empty():
    wellcontent = WellContent(quantities={"Compound_1": 5, "Compound_2": 10}, volume=25)
    wellcontent.make_empty()
    assert wellcontent.volume == 0
    assert wellcontent.quantities == {}


def test_components_as_string():
    assert wellcontent.components_as_string() == "Compound_1 Compound_2"


# From test_transfer


def test_TransferError():
    with pytest.raises(ValueError):
        raise TransferError()


source = Plate96(name="Source")
destination = Plate96(name="Destination")
source_well = source.wells["A1"]
destination_well = destination.wells["B2"]
volume = 25 * 10 ** (-6)
transfer = Transfer(source_well, destination_well, volume)


def test_to_plain_string():
    assert (
        transfer.to_plain_string()
        == "Transfer 2.50E-05L from Source A1 into Destination B2"
    )


def test_to_short_string():
    assert (
        transfer.to_short_string()
        == "Transfer 2.50E-05L (Source-A1) -> (Destination-B2)"
    )


def test_with_new_volume():
    new_volume = 50 * 10 ** (-7)
    new_transfer = transfer.with_new_volume(new_volume)
    assert new_transfer.volume == new_volume


def test_apply():
    with pytest.raises(ValueError):
        transfer.apply()

    source_2 = Plate96(name="Source_2")
    source_2.wells["A1"].add_content({"Compound_1": 1}, volume=5 * 10 ** (-6))
    destination_2 = Plate96(name="Destination_2")
    transfer_2 = Transfer(source_2.wells["A1"], destination_2.wells["B2"], volume)

    with pytest.raises(ValueError):
        transfer_2.apply()

    source_2.wells["A1"].add_content({"Compound_1": 1}, volume=25 * 10 ** (-6))
    destination_2.wells["B2"].capacity = 3 * 10 ** (-6)
    with pytest.raises(ValueError):
        transfer_2.apply()

    destination_2.wells["B2"].capacity = 50 * 10 ** (-6)
    transfer_2.apply()
    assert destination_2.wells["B2"].volume == volume


def test___repr__():
    assert (
        transfer.__repr__() == "Transfer 2.50E-05L from Source A1 into Destination B2"
    )


# From test_helper_functions
from teemi.build.containers_wells_picklists import (
    compute_rows_columns,
    rowname_to_number,
    number_to_rowname,
    wellname_to_coordinates,
    wellname_to_index,
    coordinates_to_wellname,
)


def invert_sublists(l):
    return [sl[::-1] for sl in l]


@pytest.mark.parametrize(
    "num_wells, expected",
    [(48, (6, 8)), (96, (8, 12)), (384, (16, 24)), (1536, (32, 48))],
)
def test_compute_rows_columns(num_wells, expected):
    assert compute_rows_columns(num_wells) == expected


rowname_data = [("A", 1), ("E", 5), ("AA", 27), ("AE", 31)]


@pytest.mark.parametrize("rowname, expected", rowname_data)
def test_rowname_to_number(rowname, expected):
    assert rowname_to_number(rowname) == expected


@pytest.mark.parametrize("number, expected", invert_sublists(rowname_data))
def test_number_to_rowname(number, expected):
    assert number_to_rowname(number) == expected


coordinates_data = [
    ("A1", (1, 1)),
    ("C2", (3, 2)),
    ("C04", (3, 4)),
    ("H11", (8, 11)),
    ("AA7", (27, 7)),
    ("AC07", (29, 7)),
]


@pytest.mark.parametrize("wellname, expected", coordinates_data)
def test_wellname_to_coordinates(wellname, expected):
    assert wellname_to_coordinates(wellname) == expected


coord_to_name_data = [
    ((1, 1), "A1"),
    ((3, 2), "C2"),
    ((3, 4), "C4"),
    ((8, 11), "H11"),
    ((27, 7), "AA7"),
    ((29, 7), "AC7"),
]


@pytest.mark.parametrize("coords, expected", coord_to_name_data)
def test_coordinates_to_wellname(coords, expected):
    assert coordinates_to_wellname(coords) == expected


wellname_data = [
    ("A5", 96, "row", 5),
    ("A5", 96, "column", 33),
    ("C6", 96, "row", 30),
    ("C6", 96, "column", 43),
    ("C6", 384, "row", 54),
    ("C6", 384, "column", 83),
]
inverted_wellname_data = [[s[-1], s[1], s[2], s[0]] for s in wellname_data]


@pytest.mark.parametrize("wellname, nwells, direction, expected", wellname_data)
def test_wellname_to_index(wellname, nwells, direction, expected):
    assert wellname_to_index(wellname, nwells, direction) == expected


# ---------------------------------------------------------------------------
# Additional coverage tests (imported from the non-deprecated legacy module)
# ---------------------------------------------------------------------------
import numpy as np

from teemi.legacy.build import containers_wells_picklists as cwp


def _filled_plate(name="Src"):
    plate = cwp.Plate96(name=name)
    plate["A1"].add_content({"dna": 10}, volume=20, unit_volume="uL")
    plate["A2"].add_content({"primer": 4}, volume=10, unit_volume="uL")
    return plate


def test_legacy_iterate_sources_tree_recurses_through_parent_wells():
    plate = cwp.Plate96(name="P")
    grandparent, parent, child = plate["A1"], plate["A2"], plate["A3"]
    parent.sources = [grandparent]
    child.sources = [parent, "external_reagent"]

    assert list(child.iterate_sources_tree()) == [
        grandparent,
        parent,
        "external_reagent",
        child,
    ]


def test_legacy_iterate_sources_tree_yields_applied_transfers():
    source = _filled_plate()
    destination = cwp.Plate96(name="Dst")
    transfer = cwp.Transfer(source["A1"], destination["B2"], 5e-6)
    transfer.apply()

    assert list(destination["B2"].iterate_sources_tree()) == [transfer, destination["B2"]]


def test_legacy_add_content_over_capacity_raises():
    well = cwp.Plate96(name="P")["A1"]
    well.capacity = 50e-6
    well.add_content({"dna": 1}, volume=40, unit_volume="uL")

    with pytest.raises(cwp.TransferError, match="brings volume over capacity"):
        well.add_content({"dna": 1}, volume=20, unit_volume="uL")

    assert well.volume == pytest.approx(40e-6)
    assert well.content.quantities == {"dna": 1}


def test_legacy_subtract_content_removes_fully_subtracted_component():
    well = cwp.Plate96(name="P")["A1"]
    well.add_content({"X": 5, "Y": 2}, volume=10e-6)

    well.subtract_content({"X": 5, "Y": 1}, volume=4e-6)

    assert well.content.quantities == {"Y": 1}
    assert well.volume == pytest.approx(6e-6)


def test_legacy_well_coordinates_and_ordering():
    plate = cwp.Plate96(name="P")
    assert plate["C5"].coordinates == (3, 5)
    assert plate["H12"].coordinates == (8, 12)
    # wells are ordered by their string representation "(plate-wellname)"
    assert plate["A2"] < plate["B1"]
    assert not plate["B1"] < plate["A2"]
    assert sorted([plate["B1"], plate["A3"], plate["A10"]]) == [
        plate["A10"],
        plate["A3"],
        plate["B1"],
    ]


def test_legacy_find_unique_well_by_condition():
    plate = _filled_plate()

    def has_content(well):
        return not well.is_empty

    with pytest.raises(cwp.NoUniqueWell, match="several wells"):
        plate.find_unique_well_by_condition(has_content)
    with pytest.raises(cwp.NoUniqueWell, match="No wells found"):
        plate.find_unique_well_by_condition(lambda well: well.volume > 1)
    assert plate.find_unique_well_by_condition(lambda w: w.name == "A2") is plate["A2"]
    assert plate.find_unique_well_containing("primer") is plate["A2"]


def test_legacy_list_wells_in_row_by_letter():
    wells = cwp.Plate2x4(name="P").list_wells_in_row("B")
    assert [well.name for well in wells] == ["B1", "B2", "B3", "B4"]


def test_legacy_wells_grouped_by_data_field_sorted_and_ignoring_none():
    plate = cwp.Plate2x4(
        name="P",
        wells_data={"A1": {"sample": "b"}, "A2": {"sample": "a"}, "B1": {"sample": "b"}},
    )

    grouped = plate.wells_grouped_by(data_field="sample", sort_keys=True, ignore_none=True)
    assert [(key, [w.name for w in wells]) for key, wells in grouped] == [
        ("a", ["A2"]),
        ("b", ["A1", "B1"]),
    ]

    grouped_with_none = plate.wells_grouped_by(data_field="sample")
    assert [key for key, _ in grouped_with_none] == ["b", "a", None]
    assert len(grouped_with_none[2][1]) == 5


def test_legacy_iter_wells_by_column():
    names = [well.name for well in cwp.Plate2x4().iter_wells(direction="column")]
    assert names == ["A1", "B1", "A2", "B2", "A3", "B3", "A4", "B4"]


def test_legacy_to_pandas_dataframe_with_fields():
    plate = cwp.Plate2x4(name="P")
    plate["A2"].add_content({"dna": 1}, volume=3)

    dataframe = plate.to_pandas_dataframe(fields=["name", "row"], direction="column")

    assert list(dataframe.columns) == ["name", "row"]
    assert list(dataframe["name"]) == ["A1", "B1", "A2", "B2", "A3", "B3", "A4", "B4"]
    assert list(dataframe["row"]) == [1, 2, 1, 2, 1, 2, 1, 2]


def test_legacy_plate_repr():
    assert repr(cwp.Plate96(name="Source")) == "Plate96(Source)"
    assert repr(cwp.Plate2x4()) == "Plate2x4(None)"


def test_legacy_add_transfer_from_parameters():
    source, destination = _filled_plate(), cwp.Plate96(name="Dst")
    picklist = cwp.PickList()

    picklist.add_transfer(
        source_well=source["A1"],
        destination_well=destination["C3"],
        volume=2e-6,
        data={"tip": "new"},
    )

    (transfer,) = picklist.transfers_list
    assert isinstance(transfer, cwp.Transfer)
    assert transfer.source_well is source["A1"]
    assert transfer.destination_well is destination["C3"]
    assert transfer.volume == 2e-6
    assert transfer.data == {"tip": "new"}


def _two_transfer_picklist():
    source, destination = _filled_plate(), cwp.Plate96(name="Dst")
    picklist = cwp.PickList()
    picklist.add_transfer(source_well=source["A1"], destination_well=destination["B1"], volume=5e-6)
    picklist.add_transfer(source_well=source["A2"], destination_well=destination["B1"], volume=5e-6)
    return picklist, source, destination


def test_legacy_simulate_not_inplace_returns_copies_and_leaves_plates_untouched():
    picklist, source, destination = _two_transfer_picklist()

    new_plates = picklist.simulate(inplace=False)

    assert set(new_plates) == {source, destination}
    new_destination_well = new_plates[destination]["B1"]
    assert new_destination_well.volume == pytest.approx(10e-6)
    assert new_destination_well.content.quantities == pytest.approx({"dna": 2.5, "primer": 2})
    assert new_plates[source]["A1"].volume == pytest.approx(15e-6)
    # the original plates are not modified
    assert destination["B1"].is_empty
    assert source["A1"].volume == pytest.approx(20e-6)


def test_legacy_simulate_inplace_modifies_plates():
    picklist, source, destination = _two_transfer_picklist()

    assert picklist.simulate(inplace=True) is None

    assert destination["B1"].volume == pytest.approx(10e-6)
    assert destination["B1"].content.quantities == pytest.approx({"dna": 2.5, "primer": 2})
    assert source["A2"].volume == pytest.approx(5e-6)
    assert source["A2"].content.quantities == pytest.approx({"primer": 2})


def test_legacy_sorted_by_callable():
    plate = cwp.Plate96(name="P")
    picklist = cwp.PickList(
        [
            cwp.Transfer(plate["A1"], plate["B1"], 3),
            cwp.Transfer(plate["A2"], plate["B1"], 1),
            cwp.Transfer(plate["A3"], plate["B1"], 2),
        ]
    )

    sorted_picklist = picklist.sorted_by(lambda transfer: transfer.volume)

    assert [t.volume for t in sorted_picklist.transfers_list] == [1, 2, 3]
    assert sorted_picklist.data == {"parent": picklist}


@pytest.mark.xfail(
    raises=KeyError,
    strict=True,
    reason=(
        "Bug: in PickList.sorted_by the nested `def sorting_method` shadows the "
        "string argument, so transfer.__dict__ is indexed with the function itself"
    ),
)
def test_legacy_sorted_by_attribute_name():
    plate = cwp.Plate96(name="P")
    picklist = cwp.PickList(
        [
            cwp.Transfer(plate["B1"], plate["H1"], 1),
            cwp.Transfer(plate["A1"], plate["H1"], 1),
        ]
    )

    sorted_picklist = picklist.sorted_by("source_well")

    assert [t.source_well.name for t in sorted_picklist.transfers_list] == ["A1", "B1"]


def test_legacy_enforce_maximum_dispense_volume_with_remainder():
    plate = cwp.Plate96(name="P")
    picklist = cwp.PickList([cwp.Transfer(plate["A1"], plate["B1"], 25, data={"k": 1})])

    split = picklist.enforce_maximum_dispense_volume(10)

    assert [t.volume for t in split.transfers_list] == [10, 10, 5]
    assert all(t.source_well is plate["A1"] for t in split.transfers_list)
    assert all(t.destination_well is plate["B1"] for t in split.transfers_list)
    assert all(t.data == {"k": 1} for t in split.transfers_list)
    assert split.total_transferred_volume() == 25


def test_legacy_flowbot_instructions():
    source, destination = cwp.Plate96(name="1"), cwp.Plate96(name="5")
    transfer_1 = cwp.Transfer(source["A1"], destination["B2"], 20)
    transfer_2 = cwp.Transfer(source["C3"], destination["D4"], 50.7)

    assert transfer_1.to_flowbot_instructions() == "1:A1, 5:B2, 20 "
    picklist = cwp.PickList([transfer_1, transfer_2])
    assert picklist.to_flowbot_instructions_string() == "1:A1, 5:B2, 20 \n1:C3, 5:D4, 50.7 "


def test_legacy_rowname_to_number_rejects_invalid_name():
    with pytest.raises(ValueError):
        cwp.rowname_to_number("a")


@pytest.mark.xfail(
    raises=AssertionError,
    strict=True,
    reason=(
        "Bug: rowname_to_number catches IndexError, but str.index raises "
        "ValueError('substring not found'), so the friendly message is never used"
    ),
)
def test_legacy_rowname_to_number_invalid_name_message():
    with pytest.raises(ValueError, match="is not a valid row name"):
        cwp.rowname_to_number("a")


def test_legacy_wellname_to_index_invalid_direction():
    with pytest.raises(ValueError, match="`direction` must be in"):
        cwp.wellname_to_index("A1", 96, direction="diagonal")


def test_legacy_index_to_row_column():
    assert cwp.index_to_row_column(14, 96, direction="row") == (2, 2)
    assert cwp.index_to_row_column(14, 96, direction="column") == (6, 2)
    with pytest.raises(ValueError, match="`direction` must be in"):
        cwp.index_to_row_column(1, 96, direction="diagonal")


def test_legacy_replace_nans_in_dict_nested():
    dct = {"a": np.nan, "b": {"c": np.nan, "d": 1.5}, "e": "text"}

    cwp.replace_nans_in_dict(dct, replace_by="missing")

    assert dct == {"a": "missing", "b": {"c": "missing", "d": 1.5}, "e": "text"}


def test_legacy_plate_to_dict_replaces_nans_in_well_data():
    plate = cwp.Plate2x4(name="P", wells_data={"A1": {"conc": np.nan, "id": 7}})

    dct = plate.to_dict()
    assert dct["wells"]["A1"]["conc"] == "null"
    assert dct["wells"]["A1"]["id"] == 7

    raw = plate.to_dict(replace_nans_by=None)
    assert raw["wells"]["A1"]["conc"] is np.nan
