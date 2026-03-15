# MIT License
#
# Copyright (c) 2019 Global BioFoundry Alliance
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
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# This script was writting by the Global BioFoundry Alliance and it is used
# as a foundation building easy to use functions as the RobotAssembly class.
# The original code can be found at:
# https://github.com/Global-Biofoundries-Alliance/SynBioPython

from copy import deepcopy
import math
import re
import numpy as np
from collections import OrderedDict
import pandas
from copy import deepcopy

# from Well.py
"""This module contains a generic class for a well."""


class Well:
    """Generic class for a well.

    :param plate: The plate on which the well is located
    :param row: The well's row (a number, starting from 0)
    :param column: The well's column (a number, starting from 0)
    :param name: The well's name, for instance "A1"
    :param data: A dictionary storing data on the well, used in algorithms and reports.
    """

    capacity = None
    dead_volume_per_transfer_class = None

    def __init__(self, plate, row, column, name, data=None):
        self.plate = plate
        self.row = row
        self.column = column
        self.name = name
        self.data = data or {}
        self.sources = []
        self.content = WellContent()

    @property
    def volume(self):
        """Return volume."""
        return self.content.volume

    def iterate_sources_tree(self):
        """Iterate through the tree of sources."""
        for source in self.sources:
            if isinstance(source, Well):
                for parent in source.iterate_sources_tree():
                    yield parent
            else:
                yield source
        yield self

    def add_content(self, components_quantities, volume=None, unit_volume="L"):
        """Add content to well.

        :param components_quantities: Dictionary of components and quantities
          (default: gram). Example `{"Compound_1": 5}`.
        :param volume: Volume (default: liter).
        :param unit_volume: Unit of volume (default: liter). Options: liter (L),
            milliliter (mL), microliter (uL), nanoliter (nL).
        """
        volume = volume * unit_factors[unit_volume]
        if volume > 0:
            final_volume = self.content.volume + volume
            if (self.capacity is not None) and (final_volume > self.capacity):
                raise TransferError(
                    "Transfer of %.2e L to %s brings volume over capacity."
                    % (volume, self)
                )
            self.content.volume = final_volume
        for component, quantity in components_quantities.items():
            if component not in self.content.quantities:
                self.content.quantities[component] = 0
            self.content.quantities[component] += quantity

    def subtract_content(self, components_quantities, volume=0):
        """Subtract content from well."""
        if volume > 0:
            if volume > self.volume:
                raise TransferError(
                    (
                        "Subtraction of %.2e L from %s is impossible."
                        " Current volume: %.2e L"
                    )
                    % (volume, self, self.volume)
                )
            self.content.volume -= volume
        for component, quantity in components_quantities.items():
            if self.content.quantities[component] == quantity:
                self.content.quantities.pop(component)
            else:
                self.content.quantities[component] -= quantity

    def empty_completely(self):
        """Empty the well."""
        self.content.quantities = {}
        self.content.volume = 0

    @property
    def coordinates(self):
        """Return (well.row, well.column)."""
        return (self.row, self.column)

    @property
    def is_empty(self):
        """Return true if the well's volume is 0."""
        return self.volume == 0

    def __repr__(self):
        return "(%s-%s)" % (self.plate.name, self.name)

    def pretty_summary(self):
        """Return a summary string of the well."""
        data = "\n    ".join(
            [""] + [("%s: %s" % (key, value)) for key, value in self.data.items()]
        )
        content = "\n    ".join(
            [""]
            + [
                ("%s: %s" % (key, value))
                for key, value in self.content.quantities.items()
            ]
        )
        return (
            "{self}\n"
            "  Volume: {self.volume}\n"
            "  Content: {content}\n"
            "  Metadata: {data}"
        ).format(self=self, content=content, data=data)

    def to_dict(self):
        """Convert well to dict"""
        return dict(
            [
                ["name", self.name],
                ["content", self.content.to_dict()],
                ["row", self.row],
                ["column", self.column],
            ]
            + list(self.data.items())
        )

    def index_in_plate(self, direction="row"):
        """Return the index of the well in the plate."""
        return self.plate.wellname_to_index(self.name, direction=direction)

    def is_after(self, other, direction="row"):
        """Return whether this well is located strictly after the other well.

        Example: iterate over all free wells after the last non-free well:

        >>> direction = 'row'
        >>> last_occupied_well = plate.last_nonempty_well(direction=direction)
        >>> free_wells = (w for w in plate.iter_wells(direction=direction)
        >>>               if w.is_after(last_occupied_well))
        >>> for well in free_wells: ...
        """
        well_index = self.index_in_plate(direction=direction)
        other_index = other.index_in_plate(direction=direction)
        return well_index > other_index

    def __lt__(self, other):
        return str(self) < str(other)


# from Plate.py
class NoUniqueWell(Exception):
    """NoUniqueWell exception class."""


class Plate:
    """Base class for all plates.

    See the builtin_containers for usage classes, such as generic microplate
    classes (Plate96, Plate384, etc).

    :param name: Name or ID of the Plate as it will appear in strings and reports
    :param wells_data: A dict {"A1": {data}, "A2": ...}.
        The format of the data is left free
    :param plate_data: plate data
    """

    well_class = Well

    def __init__(self, name=None, wells_data=None, plate_data=None):

        self.name = name
        self.data = plate_data or {}
        self.wells_data = wells_data or {}
        self.num_wells = self.num_rows * self.num_columns
        self.wells = {}
        self.columns = {column: [] for column in range(1, self.num_columns + 1)}
        self.rows = {number_to_rowname(row): [] for row in range(1, self.num_rows + 1)}
        for row in range(1, self.num_rows + 1):
            for column in range(1, self.num_columns + 1):
                wellname = coordinates_to_wellname((row, column))
                data = self.wells_data.get(wellname, {})
                well = self.well_class(
                    plate=self,
                    row=row,
                    column=column,
                    name=wellname,
                    data=data,
                )
                self.wells[wellname] = well
                self.columns[column] += [wellname]
                self.rows[number_to_rowname(row)] += [wellname]

    def __getitem__(self, k):
        """Return e.g. well A1's dict when calling `myplate['A1']`."""
        return self.wells[k]

    def find_unique_well_by_condition(self, condition):
        """Return the unique well of the plate satisfying the condition.

        The ``condition`` method should have a signature of Well=>True/False.

        Raises a NoUniqueWell error if 0 or several wells satisfy the condition.
        """
        wells = [well for name, well in self.wells.items() if condition(well)]
        if len(wells) > 1:
            raise NoUniqueWell("Query returned several wells: %s" % wells)
        if len(wells) == 0:
            raise NoUniqueWell("No wells found matching the condition")
        return wells[0]

    def find_unique_well_containing(self, query):
        """Return the unique well whose content contains the query."""

        def condition(well):
            return query in well.content.quantities.keys()

        return self.find_unique_well_by_condition(condition)

    def list_well_data_fields(self):
        """Return all fields used in well data in the plate."""
        return sorted(list(set(field for well in self for field in well.data.keys())))

    def return_column(self, column_number):
        """Return the list of all wells of the plate in the given column."""
        return [self.wells[wellname] for wellname in self.columns[column_number]]

    def list_wells_in_column(self, column_number):
        """Return the list of all wells of the plate in the given column.

        Examples:

        >>> for well in plate.list_wells_in_column(5):
        >>>      print(well.name)
        """
        return [well for well in self.iter_wells() if well.column == column_number]

    def return_row(self, row):
        """Return the list of all wells of the plate in the given row.

        The `row` can be either a row number (1,2,3) or row letter(s) (A,B,C).
        """
        if isinstance(row, int):
            row = number_to_rowname(row)
        return [self.wells[wellname] for wellname in self.rows[row]]

    def list_wells_in_row(self, row):
        """Return the list of all wells of the plate in the given row.

        The `row` can be either a row number (1,2,3) or row letter(s) (A,B,C).

        Examples:

        >>> for well in plate.list_wells_in_row("H"):
        >>>      print(well.name)

        """
        if isinstance(row, str):
            row = rowname_to_number(row)
        return [well for well in self.iter_wells() if well.row == row]

    def list_filtered_wells(self, well_filter):
        """List filtered wells.

        Examples:

        >>> def condition(well):
        >>>     return well.volume > 50
        >>> for well in myplate.list_filtered_wells(condition):
        >>>     print(well.name)
        """
        return list(filter(well_filter, self.wells.values()))

    def wells_grouped_by(
        self,
        data_field=None,
        key=None,
        sort_keys=False,
        ignore_none=False,
        direction_of_occurence="row",
    ):
        """Return wells grouped by key."""
        if key is None:

            def key(well):
                return well.data.get(data_field, None)

        dct = OrderedDict()
        for well in self.iter_wells(direction=direction_of_occurence):
            well_key = key(well)
            if well_key not in dct:
                dct[well_key] = [well]
            else:
                dct[well_key].append(well)
        if ignore_none:
            dct.pop(None, None)
        keys = dct.keys()
        if sort_keys:
            keys = sorted(keys)
        return [(k, dct[k]) for k in keys]

    def get_well_at_index(self, index, direction="row"):
        """Return the well at the corresponding index.

        Examples:

        >>> plate.get_well_at_index(1)  # well A1
        >>> plate.get_well_at_index(2)  # well A2
        >>> plate.get_well_at_index(2, direction="column")  # well B1
        """
        return self[self.index_to_wellname(index, direction=direction)]

    def index_to_wellname(self, index, direction="row"):
        """Return the name of the well at the corresponding index.

        Examples:

        >>> plate.index_to_wellname(1)  # "A1"
        >>> plate.get_well_at_index(2)  # "A2"
        >>> plate.get_well_at_index(2, direction="column")  # "B1"
        """
        return index_to_wellname(index, self.num_wells, direction=direction)

    def wellname_to_index(self, wellname, direction="row"):
        """Return the index of the well in the plate.

        Examples:
        >>> plate.wellname_to_index("A1")  # 1
        >>> plate.wellname_to_index("A2")  # 2
        >>> plate.wellname_to_index("A1", direction="column")  # 9 (8x12 plate)
        """
        return wellname_to_index(wellname, self.num_wells, direction=direction)

    def wells_sorted_by(self, sortkey):
        """Return wells sorted by sortkey"""
        return (e for e in sorted(self.wells.values(), key=sortkey))

    def iter_wells(self, direction="row"):
        """Iter through the wells either by row or by column.

        Examples:

        >>> for well in plate.iter_wells():
        >>>     print (well.name)
        """
        if direction == "row":
            return self.wells_sorted_by(lambda w: (w.row, w.column))
        else:
            return self.wells_sorted_by(lambda w: (w.column, w.row))

    def to_dict(self, replace_nans_by="null"):
        """Convert plate to dict."""
        dct = {
            "data": self.data,
            "wells": {well.name: well.to_dict() for well in self.wells.values()},
        }
        if replace_nans_by is not None:
            replace_nans_in_dict(dct, replace_by=replace_nans_by)
        return dct

    def to_pandas_dataframe(self, fields=None, direction="row"):
        """Return a dataframe with the info on each well."""
        dataframe = pandas.DataFrame.from_records(self.to_dict()["wells"]).T
        by = ["row", "column"] if direction == "row" else ["column", "row"]
        dataframe = dataframe.sort_values(by=by)
        if fields is not None:
            dataframe = dataframe[fields]
        return dataframe

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.name)


# From PickList.py
"""Classes to represent picklists and liquid transfers in general."""


class PickList:
    """Representation of a list of well-to-well transfers.

    :param transfers_list: A list of Transfer objects that will be part of the same
        dispensing operation, in the order in which they are meant to be simulated.
    :param data: A dict with information on the picklist.
    """

    def __init__(self, transfers_list=(), data=None):

        self.transfers_list = list(transfers_list)
        self.data = {} if data is None else data

    def add_transfer(
        self,
        source_well=None,
        destination_well=None,
        volume=None,
        data=None,
        transfer=None,
    ):
        """Add a transfer to the picklist's tranfers list.

        You can either provide a ``Transfer`` object with the ``transfer``
        parameter, or the parameters.
        """
        if transfer is None:
            transfer = Transfer(
                source_well=source_well,
                destination_well=destination_well,
                volume=volume,
                data=data,
            )
        self.transfers_list.append(transfer)

    def to_plain_string(self):
        """Return the list of transfers in human-readable format."""
        return "\n".join(transfer.to_plain_string() for transfer in self.transfers_list)

    def to_plain_textfile(self, filename):
        """Write the picklist in a file in a human reable format."""
        with open(filename, "w+") as f:
            f.write(self.to_plain_string())

    def simulate(self, content_field="content", inplace=True):
        """Simulate the execution of the picklist."""

        if not inplace:
            all_plates = set(
                plate
                for transfer in self.transfers_list
                for plate in [
                    transfer.source_well.plate,
                    transfer.destination_well.plate,
                ]
            )
            new_plates = {plate: deepcopy(plate) for plate in all_plates}

            new_transfer_list = []
            for transfer in self.transfers_list:
                new_source_plate = new_plates[transfer.source_well.plate]
                new_dest_plate = new_plates[transfer.destination_well.plate]
                new_source_well = new_source_plate.wells[transfer.source_well.name]
                new_dest_well = new_dest_plate.wells[transfer.destination_well.name]
                new_transfer_list.append(
                    Transfer(
                        volume=transfer.volume,
                        source_well=new_source_well,
                        destination_well=new_dest_well,
                    )
                )

            new_picklist = PickList(transfers_list=new_transfer_list)
            new_picklist.simulate(
                content_field=content_field,
                inplace=True,
            )
            return new_plates

        else:
            for transfer in self.transfers_list:
                transfer.apply()
            return None

    def restricted_to(
        self, transfer_filter=None, source_well=None, destination_well=None
    ):
        """Return a version of the picklist restricted to transfers with the
        right source/destination well.

        You can provide ``source_well`` and ``destination_well`` or
        alternatively just a function ``transfer_filter`` with signature
        (transfer)=>True/False that will be used to filter out transfers
        (for which it returns false).
        """
        if transfer_filter is None:

            def transfer_filter(tr):
                source_well_is_ok = (source_well is None) or (
                    source_well == tr.source_well
                )
                dest_well_is_ok = (destination_well is None) or (
                    destination_well == tr.destination_well
                )
                return source_well_is_ok and dest_well_is_ok

        transfers = [tr for tr in self.transfers_list if transfer_filter(tr)]
        return PickList(transfers, data={"parent": self})

    def sorted_by(self, sorting_method="source_well"):
        """Return a new version of the picklist sorted by some parameter.

        The ``sorting_method`` is either the name of an attribute of the
        transfers, such as "source_well", or a function f(transfer) -> value.
        """
        if not hasattr(sorting_method, "__call__"):

            def sorting_method(transfer):
                return transfer.__dict__[sorting_method]

        return PickList(
            sorted(self.transfers_list, key=sorting_method),
            data={"parent": self},
        )

    def total_transferred_volume(self):
        """Return the sum of all volumes from all transfers."""
        return sum([transfer.volume for transfer in self.transfers_list])

    def enforce_maximum_dispense_volume(self, max_dispense_volume):
        """Return a new picklist were every too-large dispense is broken down
        into smaller dispenses."""
        transfers = []
        for trf in self.transfers_list:
            n_additional_dispense = int(trf.volume / max_dispense_volume)
            rest = trf.volume - n_additional_dispense * max_dispense_volume
            for _ in range(n_additional_dispense):
                transfers.append(trf.with_new_volume(max_dispense_volume))
            if rest > 0:
                transfers.append(trf.with_new_volume(rest))
        return PickList(transfers_list=transfers)

    def to_flowbot_instructions_string(self):
        # Made to accomodate flowbot instructions - not part of synbiopython
        """Return the list of transfers in Flowbot format."""
        return "\n".join(
            transfer.to_flowbot_instructions() for transfer in self.transfers_list
        )

    def __add__(self, other):
        return PickList(self.transfers_list + other.transfers_list)

    @staticmethod
    def merge_picklists(picklists_list):
        """Merge the list of picklists into a single picklist.

        The transfers in the final picklist are the concatenation of the
        transfers in the different picklists, in the order in which they appear
        in the list.
        """
        return sum(picklists_list, PickList([]))


# From Transfer.py
class TransferError(ValueError):
    pass


class Transfer:
    """Class representing a transfer from a source well to a destination well.

    :param source_well: A Well object from which to transfer.
    :param destination_well: A Well object to which to transfer.
    :param volume: Volume to be transferred, expressed in liters.
    :param data: A dict containing any useful information about the transfer.
        This information can be used later e.g. as parameters for the transfer
        when exporting a picklist.
    """

    def __init__(self, source_well, destination_well, volume, data=None):

        self.volume = volume
        self.source_well = source_well
        self.destination_well = destination_well
        self.data = data

    def to_plain_string(self):
        """Return "Transfer {volume}L from {source_well} into {dest_well}"."""
        return (
            "Transfer {self.volume:.02E}L from {self.source_well.plate.name} "
            "{self.source_well.name} into "
            "{self.destination_well.plate.name} "
            "{self.destination_well.name}"
        ).format(self=self)

    def to_short_string(self):
        """Return "Transfer {volume}L {source_well} -> {dest_well}"."""
        return (
            "{self.__class__.__name__} {self.volume:.02E}L {self.source_well} -> {self.destination_well}"
        ).format(self=self)

    def with_new_volume(self, new_volume):
        """Return a version of the transfer with a new volume."""
        return self.__class__(
            source_well=self.source_well,
            destination_well=self.destination_well,
            volume=new_volume,
            data=self.data,
        )

    def apply(self):
        # error_prefix = "%s error:" % self.to_short_string()

        if self.source_well.is_empty:
            raise TransferError("Source well is empty!")

        #  pre-check in both source and destination wells that transfers
        #  are valid
        if self.volume > self.source_well.volume:
            raise TransferError(
                ("Subtraction of %.2e L from %s impossible." " Current volume: %.2e L")
                % (self.volume, self, self.source_well.volume)
            )
        final_destination_volume = self.destination_well.volume + self.volume
        if (self.destination_well.capacity is not None) and (
            final_destination_volume > self.destination_well.capacity
        ):
            raise TransferError(
                "Transfer of %.2e L from %s to %s brings volume over capacity."
                % (self.volume, self, self.destination_well)
            )

        #  If you arrive here, it means that the transfer is valid, do it.
        factor = float(self.volume) / self.source_well.volume

        quantities_transferred = {
            component: quantity * factor
            for component, quantity in self.source_well.content.quantities.items()
        }
        self.destination_well.add_content(quantities_transferred, volume=self.volume)
        self.source_well.subtract_content(quantities_transferred, volume=self.volume)
        if self not in self.destination_well.sources:
            self.destination_well.sources.append(self)

    def to_flowbot_instructions(self):
        """
        Return Flowbot instructions.

        Example:
        
        .. code-block:: none

            source, destination, volume
            4:A3, 4:A6, 20
            3:A1, 7, 50.7
            2:A, 2:B-F, 100
        """
        # Made to accommodate flowbot instructions - not part of synbiopython
        return (
            "{self.source_well.plate.name}:"
            "{self.source_well.name},"
            " {self.destination_well.plate.name}:"
            "{self.destination_well.name}, {self.volume} "
        ).format(self=self)


    def __repr__(self):
        """Return  "Transfer {volume}L from {source_well} into {dest_well}"."""
        return self.to_plain_string()


# From WellContent.py
class WellContent:
    """Class to represent the volume and quantities of a well.

    Having the well content represented as a separate object makes it possible
    to have several wells share the same content, e.g. in throughs.
    """

    def __init__(self, quantities=None, volume=0):
        if quantities is None:
            quantities = {}
        self.volume = volume
        self.quantities = quantities

    def concentration(self, component=None, default=0):
        """Return concentration of component."""
        if self.quantities == {}:
            return default
        if self.volume == 0:
            return default
        if component is None:
            component = list(self.quantities.keys())[0]
        if component not in self.quantities:
            return default
        return 1.0 * self.quantities[component] / self.volume

    def to_dict(self):
        """Return a dict {volume: 0.0001, quantities: {...:...}}."""
        return {"volume": self.volume, "quantities": self.quantities}

    def make_empty(self):
        """Empty the well."""
        self.volume = 0
        self.quantities = {}

    def components_as_string(self, separator=" "):
        """Return a string representation of what's in the well mix."""
        return separator.join(sorted(self.quantities.keys()))


# From builtin_containers.py
"""Classes to represent plates"""


class Plate96(Plate):
    """Base class for standard 96-well plates"""

    num_rows = 8
    num_columns = 12


class Plate2x4(Plate):
    """Class for 8-well (2 x 4) plates such as colony plating plates"""

    num_rows = 2
    num_columns = 4


# From helper_functions.py
def compute_rows_columns(num_wells):
    """Convert 96->(8,12), 384->(16,24), etc."""
    a = math.sqrt(num_wells / 6)
    n_rows = int(round(2 * a))
    n_columns = int(round(3 * a))
    return n_rows, n_columns


def rowname_to_number(name):
    "Convert A->1 Z->26 AA->27 etc."
    if len(name) == 2:
        return 26 * rowname_to_number(name[0]) + rowname_to_number(name[1])
    try:
        return "ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(name) + 1
    except IndexError:
        raise ValueError(name + " is not a valid row name.")


def number_to_rowname(number):
    "Convert 1->A 26->Z 27->AA etc."
    if number > 26:
        return number_to_rowname(int(number / 26)) + number_to_rowname(number % 26)
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[number - 1]


def wellname_to_coordinates(wellname):
    """Convert A1->(1,1), H11->(8, 11), etc."""
    rowname, colname = re.match("([a-zA-Z]+)([0-9]+)", wellname).groups()
    return rowname_to_number(rowname), int(colname)


def coordinates_to_wellname(coords):
    """Convert (1,1)->A1, (4,3)->D3, (12, 12)->H12, etc."""
    row, column = coords
    return number_to_rowname(row) + str(column)


def wellname_to_index(wellname, num_wells, direction="row"):
    """Convert e.g. A1..H12 into 1..96
    direction is either row for A1 A2 A3... or column for A1 B1 C1 D1 etc.

    :param wellname: the name of the well
    :param num_wells: number of wells on the plate
    :type num_wells: int
    :param direction: the direction of counting. Either "row" or "column".
    :type direction: str
    """
    n_rows, n_columns = compute_rows_columns(num_wells)
    row, column = wellname_to_coordinates(wellname)
    if direction == "row":
        return column + n_columns * (row - 1)
    if direction == "column":
        return row + n_rows * (column - 1)
    raise ValueError("`direction` must be in (row, column)")


def index_to_row_column(index, num_wells, direction="row"):
    n_rows, n_columns = compute_rows_columns(num_wells)
    if direction == "row":
        row = 1 + int((index - 1) / n_columns)
        column = 1 + ((index - 1) % n_columns)
    elif direction == "column":
        row, column = 1 + ((index - 1) % n_rows), 1 + int((index - 1) / n_rows)
    else:
        raise ValueError("`direction` must be in (row, column)")
    return row, column


def index_to_wellname(index, num_wells, direction="row"):
    """Convert e.g. 1..96 into A1..H12

    :param index: the index of the well
    :type index: int
    :param num_wells: number of wells on the plate
    :type num_wells: int
    :param direction: the direction of counting. Either "row" or "column".
    :type direction: str
    """
    row, column = index_to_row_column(index, num_wells, direction)
    return coordinates_to_wellname((row, column))


# From tools.py
unit_factors = {
    prefix + unit: factor
    for unit in "glL"
    for prefix, factor in [("", 1), ("m", 1e-3), ("u", 1e-6), ("n", 1e-9)]
}

volume_values_and_units = sorted(
    [(value, unit) for (unit, value) in unit_factors.items() if unit.endswith("L")]
)


def replace_nans_in_dict(dictionary, replace_by="null"):
    """Replace NaNs in a dictionary with a string.

    :param dictionary: the dictionary
    :type dictionary: dict
    :param replace_by: replacement
    :type replace_by: str
    """
    for key, value in dictionary.items():
        if isinstance(value, dict):
            replace_nans_in_dict(value, replace_by=replace_by)
        elif value is np.nan:
            dictionary[key] = replace_by
