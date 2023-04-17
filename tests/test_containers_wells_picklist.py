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
