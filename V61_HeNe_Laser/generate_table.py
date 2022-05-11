import numpy as np
import pint
from uncertainties import ufloat


class StringFormatter():
    def __call__(self, value):
        return str(value)


class PintFormatter():
    def __call__(self, value):
        magnitude = value.m
        return f'{magnitude:.2f}'


class Column:
    def __init__(self, data, formatter=None):
        self.data = data
        self.formatter = formatter or StringFormatter()

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

    @property
    def units(self):
        first_cell = self.data[0]
        if isinstance(first_cell, pint.Quantity):
            return first_cell.units
        else:
            return None

    name = '?'  # TODO


def generate_table_content(*columns):
    """
    Generate a table from a list of columns.
    """
    output_columns = []
    # columns = [Column(data=c) for c in columns]
    max_rows = max([len(c.data) for c in columns])
    for row_index in range(max_rows):
        cells = [c.get_cell(row_index) for c in columns]
        output_columns += [" & ".join(cells) + r" \\"]

    return '\n'.join(output_columns)


def generate_table(*columns):
    formatters = [PintFormatter() for c in columns]
    columns = [Column(data=c, formatter=f) for c, f in zip(columns, formatters)]
    # TODO: Assuming all-pint for now
    coltypes = ['S' for c in columns]  # TODO
    units_siunitx = [f'{c.units:Lx}' for c in columns]
    names = [f'{c.name}' for c in columns]

    output = []
    output += [r"\begin{tabular}{" f"{' '.join(coltypes)}" "}"]
    output += [r"\toprule"]
    foo = [
        f'${name}' r' \mathbin{/} ' f'{unit}$'
        for name, unit in zip(names, units_siunitx)
    ]
    output += [" & \n".join(foo) + r" \\"]

    output += [r"\midrule"]
    output += [generate_table_content(*columns)]
    output += [r"\bottomrule"]
    output += [r"\end{tabular}"]

    return '\n'.join(output)


# print(generate_table(
#     [1, 2, 3, 4],
#     [5, 6, 7, 8],
#     [9, 10, 11, 12],
# ))

ureg = pint.UnitRegistry()
r = np.arange(0, 10) * ureg('m')
I = np.arange(0, 10) * ureg('dBm')

tab = generate_table(r, I)
print(tab)
