import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from numpy.linalg import inv
import pandas as pd
import pint
import rich
from rich.console import Console
import tools
ureg = pint.UnitRegistry()
console = Console()


def A_row_from_indices(indices):
    indices = indices.strip()  # Zwischenlösung
    baserow = np.zeros(9)  # 9 = Anz. Würfel
    if '/' in indices:
        assert '|' not in indices
        # → diagonal
        # baserow[np.array(map(int,indices.split('/')))] = np.sqrt(2)
        for i in map(int, indices.split('/')):
            baserow[i-1] = np.sqrt(2)
    elif '|' in indices:
        assert '/' not in indices
        # → parallel
        # baserow[np.array(map(int,indices.split('|')))] = 1
        for i in map(int, indices.split('|')):
            baserow[i-1] = 1
    return baserow


def A_from_indices(all_indices):
    return np.row_stack(list(map(A_row_from_indices, all_indices)))


def calc_I(N, T):
    """
    Zählrate aus Anzahl N und Zeit T.
    Input: ohne Einheiten, Output: Liste mit Einheiten.
    """
    I = unp.uarray(N / T, np.sqrt(N)/T) / ureg.s
    # N *= ureg.dimensionless  # TODO
    # T *= ureg.s
    return I

# TODO: Nullmessung


def analyze(filename, I_0_getter):
    dat = pd.read_csv(filename)
    # dat.rename(str.strip)
    # dat.rename_axis()

    A = A_from_indices(dat['indices'])
    # print(A.round(1))

    N = dat['N']
    T = dat['T']

    I = calc_I(N, T)
    N *= ureg.dimensionless  # TODO
    T *= ureg.s

    I_0 = I_0_getter(dat['indices'])
    y = unp.log((I_0 / I).to('dimensionless').m)
    μ_mat = np.linalg.inv(A.T @ A) @ A.T @ y
    μ = μ_mat.mean()

    # print(f'{μ_mat=}')
    print(f"μ = {μ:.3f}")
    return µ


def get_data(filename):
    dat = pd.read_csv(filename, comment='#')

    A = A_from_indices(dat['indices'])

    N = dat['N']
    T = dat['T']

    I = calc_I(N, T)
    # NOTE: Ohne pint-pandas wird implizit zu einheitenlosen Arrays konvertiert!

    return {
        'N': N * ureg.dimensionless,
        'T': T * ureg.s,
        'I': I,
    }


def get_closest_material(µ):
    μ_LIT = {
        'Al': 0.211,
        'Pb': 1.419,
        'Fe': 0.606,
        'Messing': 0.638,
        'Delrin': 0.121,
    }

    diff_tuples = [(abs(µ - µ_lit_single) / µ_lit_single, name) for name, µ_lit_single in µ_LIT.items()]
    diff_tuples.sort()
    return diff_tuples[0]


dat_Nullmessung = get_data('dat2/Nullmessung.csv')
assert dat_Nullmessung['I'].check('1/[time]')
I_0_parallel, I_0_hauptdiag, I_0_nebendiag = dat_Nullmessung['I']


def get_I_0_for_indices(indices):
    if '/' in indices:
        if len(indices.split('/')) == 3:
            return I_0_hauptdiag
        else:
            assert len(indices.split('/')) == 2
            return I_0_nebendiag
    else:
        return I_0_parallel


# analyze('dat2/Nullmessung.csv')
µ = analyze('dat2/Würfel2.csv', get_I_0_for_indices)
# Würfel 2 SOLL: Eisen
mat_abw, mat = get_closest_material(µ)
print(f"Best fit: {mat} mit Abweichung {mat_abw:.1%}")
