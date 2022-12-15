import numpy as np
import pint
from uncertainties import ufloat


class StringFormatter():
    def __call__(self, value):
        return str(value)


class PintFormatter:
    def __init__(self, unit, decimals=1):
        self.unit = unit
        self.decimals = decimals

    def __call__(self, value):
        magnitude = value.to(self.unit).m
        # return f'{magnitude:.2f}'
        return '{0:.{decimals}f}'.format(magnitude, decimals=self.decimals)


class Column:
    # name = '?'

    def __init__(self, input_tuple):
        self.name: str = input_tuple[0]
        self.unit: pint.Unit | type = input_tuple[1]
        self.data = input_tuple[2]
        if len(input_tuple) == 4:
            assert isinstance(input_tuple[3], int)
            self.formatter = PintFormatter(self.unit, decimals=input_tuple[3])
        else:
            if self.unit == str:
                self.formatter = StringFormatter()
            else:
                self.formatter = PintFormatter(self.unit)
                # verify input [TODO: Move to separate function]
                if isinstance(self.unit, pint.Quantity) and self.unit.magnitude == 1:
                    raise TypeError(f'"{self.unit}" is not a pint.Unit. Hint: Use ureg.… instead of ureg("…")')
                assert isinstance(self.unit, pint.Unit), f'"{self.unit}" is not a pint.Unit'


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
    def header(self):
        # if self.unit == pint.Unit('kiloelectron_volt'):
        if str(self.unit) == 'kiloelectron_volt':
            # TODO: Fix upstream or generalize to all eV units
            return (f'${self.name}' r' \mathbin{/} ' r'\si[]{\kilo\electronvolt}$')
        return (
            (f'${self.name}' r' \mathbin{/} ' f'{self.unit:Lx}$')
            if isinstance(self.unit, pint.Unit) and str(self.unit) != 'dimensionless'
            else f'${self.name}$'
        )

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
    columns = [Column(ct) for ct in column_tuples]

    coltypes = ['S' for c in columns]  # TODO

    output = []
    output += [r"\begin{tabular}{" f"{' '.join(coltypes)}" "}"]
    output += [r"\toprule"]
    foo = ['{' + c.header + '}' for c in columns]
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
