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


def d_from_indices(all_indices):
    return tools.pintify(list(map(d_row_from_indices, all_indices)))

def I_0_from_indices(all_indices):
    return tools.pintify(list(map(I_0_row_from_indices, all_indices)))

def d_row_from_indices(indices):
    if '/' in indices:
        if len(indices.split('/')) == 3:
            return 3 * d_Einzelwürfel / np.sqrt(2)
        else:
            assert len(indices.split('/')) == 2
            return 2 * d_Einzelwürfel / np.sqrt(2)
    else:
        return 3 * d_Einzelwürfel


def I_0_row_from_indices(indices):
    if '/' in indices:
        if len(indices.split('/')) == 3:
            return I_0_hauptdiag
        else:
            assert len(indices.split('/')) == 2
            return I_0_nebendiag
    else:
        return I_0_parallel


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


# d_Würfel = ureg('3 cm')
d_Einzelwürfel = ureg('1 cm')


def analyze_homogen(dat, I_0_getter):
    A = A_from_indices(dat['indices'])
    # print(A.round(1))

    I_0 = I_0_getter(dat['indices'])
    d = d_from_indices(dat['indices'])
    μ_vec = (
        unp.log((I_0 / dat['I']).to('dimensionless').m) /
        d
    )

    print(f'{μ_vec=}')
    µ = abs(µ_vec).mean() # TODO: SEM
    print(f'{μ=}')
    return µ


def analyze(dat, I_0_getter):
    A = A_from_indices(dat['indices'])
    # print(A.round(1))

    I_0 = I_0_getter(dat['indices'])
    y = unp.log((I_0 / dat['I']).to('dimensionless').m)
    d = d_from_indices(dat['indices'])
    μ_vec = (
        np.linalg.inv(A.T @ A) @ A.T @ y /
        d
    )

    print(f'{μ_vec=}')
    µ = abs(µ_vec).mean()
    print(f'{μ=}')
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
        'indices': dat['indices'],
    }


μ_LIT = {
    'Al': ureg('0.211 / cm'),
    'Pb': ureg('1.419 / cm'),
    'Fe': ureg('0.606 / cm'),
    'Messing': ureg('0.638 / cm'),
    'Delrin': ureg('0.121 / cm'),
}


def get_closest_material(µ):
    diff_tuples = [(abs(µ - µ_lit_single), name) for name, µ_lit_single in µ_LIT.items()]
    diff_tuples.sort()
    return diff_tuples[0]


dat_Nullmessung = get_data('dat2/Nullmessung.csv')
assert dat_Nullmessung['I'].check('1/[time]')
I_0_parallel, I_0_hauptdiag, I_0_nebendiag = dat_Nullmessung['I']
print(f"Nullmessung: {dat_Nullmessung['I']:.2f}")

WÜRFEL = [
    {
        'num': 2,
        'material': 'Fe',  # Mampfzwerg, SanjoR
    },
]

for würfel in WÜRFEL:
    dat = get_data(f'dat2/Würfel{würfel["num"]}.csv')
    µ = analyze_homogen(dat, I_0_from_indices)
    print(f"μ = {μ:.3f}")
    mat_abw, mat = get_closest_material(µ)
    print(f"Best fit: {mat} mit Abweichung {mat_abw:.2f}")
    if mat == würfel['material']:
        print("→ wie erwartet :)")
    else:
        print(f"→ sollte {würfel['material']} sein :(")
    print(
        "Abweichung best-fit vs. µ_lit des tatsächlichen Materials:\n",
        tools.fmt_compare_to_ref(µ, µ_LIT[würfel['material']])
    )
