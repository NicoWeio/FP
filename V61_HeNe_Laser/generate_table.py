import numpy as np
import pint
from uncertainties import ufloat


class StringFormatter():
    def __call__(self, value):
        return str(value)


class PintFormatter():
    def __init__(self, unit):
        self.unit = unit

    def __call__(self, value):
        magnitude = value.to(self.unit).m
        return f'{magnitude:.2f}'


class Column:
    # name = '?'

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

    # @property
    # def units(self):
    #     first_cell = self.data[0]
    #     if isinstance(first_cell, pint.Quantity):
    #         return first_cell.units
    #     else:
    #         return None


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


def generate_table_pint(filename, *column_tuples):
    names, units, raw_columns = zip(*column_tuples, strict=True)

    # verify input
    assert all(isinstance(n, str) for n in names)
    assert all(isinstance(c, pint.Quantity) for c in raw_columns)
    for u in units:
        if isinstance(u, pint.Quantity) and u.magnitude == 1:
            raise TypeError(f'"{u}" is not a pint.Unit. Hint: Use ureg.… instead of ureg("…")')
        assert isinstance(u, pint.Unit), f'"{u}" is not a pint.Unit'
    assert all(isinstance(u, pint.Unit) for u in units)

    # generate output

    formatters = [PintFormatter(u) for c, u in zip(raw_columns, units)]
    columns = [Column(data=c, formatter=f) for c, f in zip(raw_columns, formatters)]

    coltypes = ['S' for c in columns]  # TODO
    units_siunitx = [f'{u:Lx}' for u in units]

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

    final_output = '\n'.join(output)

    # write output
    with open(filename, 'w') as f:
        f.write(final_output)


# ureg = pint.UnitRegistry()
# r = np.arange(0, 10) * ureg('m')
# I = np.arange(0, 10) * ureg('dBm')

# tab = generate_table_pint(('r', ureg.cm, r), ('I', ureg.dBm, I))
# print(tab)
