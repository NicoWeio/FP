import numpy as np
import pint
from uncertainties import ufloat

stringFormatter = str


class Column:
    def __init__(self, data):
        self.data = data
        self.formatter = stringFormatter  # TODO

    def get_cell(self, row_index):
        cell_data = self.data[row_index]
        if cell_data is not None:
            # format the cell data
            s = self.formatter(cell_data)
            # pad the cell data with spaces to the width of the column
            s = '{0: >{width}}'.format(s, width=self.max_width)
            return s
        else:
            return "–"

    @property
    def max_width(self):
        # TODO: This doesn't check placeholder widths
        return max([len(self.formatter(x)) for x in self.data])


def generate_table(*columns):
    """
    Generate a table from a list of columns.
    """
    output_columns = []
    columns = [Column(data=c) for c in columns]
    max_rows = max([len(c.data) for c in columns])
    for row_index in range(max_rows):
        cells = [c.get_cell(row_index) for c in columns]
        output_columns += [" & ".join(cells) + r" \\"]

    return '\n'.join(output_columns)


print(generate_table(
    [1, 2, 3, 4],
    [5, 6, 7, 8],
    [9, 10, 11, 12],
)
)
