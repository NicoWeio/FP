import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from numpy.linalg import inv
import pint
import rich
from rich.console import Console
import tools
ureg = pint.UnitRegistry()
console = Console()


def calc_I(N, T):
    """
    Zählrate aus Anzahl N und Zeit T.
    Input: ohne Einheiten, Output: Liste mit Einheiten.
    """
    I = unp.uarray(N / T, np.sqrt(N)/T) / ureg.s
    # N *= ureg.dimensionless  # TODO
    # T *= ureg.s
    return I


def kleinsteQuadrate(y, W, A):
    temp = np.dot(np.linalg.inv(np.dot(A.T, np.dot(W, A))), A.T)
    a = np.dot(temp, np.dot(W, y))
    a_err = np.linalg.inv(np.dot(A.T, np.dot(W, A)))
    return a, np.sqrt(np.diag(a_err))


s = np.sqrt(2)
A = np.matrix([[0, 0, 0, 0, 0, s, 0, s, 0],
               [0, 0, s, 0, s, 0, s, 0, 0],
               [0, s, 0, s, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 1, 1, 1],
               [0, 0, 0, 1, 1, 1, 0, 0, 0],
               [1, 1, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, s, 0, 0, 0, s, 0],
               [s, 0, 0, 0, s, 0, 0, 0, s],
               [0, s, 0, 0, 0, s, 0, 0, 0],
               [1, 0, 0, 1, 0, 0, 1, 0, 0],
               [0, 1, 0, 0, 1, 0, 0, 1, 0],
               [0, 0, 1, 0, 0, 1, 0, 0, 1]])

μ_LIT = {
    'Al': 0.211,
    'Pb': 1.419,
    'Fe': 0.606,
    'Messing': 0.638,
    'Delrin': 0.121,
}

MATERIALS = [
    {
        'name': 'Alu',
        'μ_lit': 0.202,
    },
    {
        'name': 'Blei',
        'μ_lit': 1.25,
    },
    {
        'name': 'Unb',
    },
    {
        'name': 'Luft',
    },
]

# TODO: █ Leermessung (Luft)
# TODO: Spektrum plotten?


# █ Würfel X → nur Aluminiumhülle
# → das nennen wir I_0 (!)
console.rule("Nullmessung")

# TODO: irgendwie Orientierung richtig einlesen oder anders handhaben – Pandas?
N_0, t_0, orient_0 = np.genfromtxt(f'dat/Nullmessung.csv', delimiter=',', skip_header=1, unpack=True)
# TODO: Nullmessung pro Material!?
I_0 = calc_I(N_0, t_0).mean() # TODO: hier wirklich mitteln?
print(f"{I_0=}")


# █ Würfel YZ

for material in MATERIALS:
    console.rule(material['name'])
    N, t = np.genfromtxt(f'dat/{material["name"]}.csv', delimiter=',', skip_header=1, unpack=True)
    I = calc_I(N, t)
    N *= ureg.dimensionless  # TODO
    t *= ureg.s

    # I_0 = ureg('184/s')  # TODO: Nullmessung pro Material
    y = unp.log((I_0 / I).to('dimensionless').m)
    μ_mat = np.linalg.inv(A.T @ A) @ A.T @ y
    μ = μ_mat.mean()

    # print(f'{μ_mat=}')
    print(f"{μ=}")

    if 'μ_lit' in material:
        # TODO: korrekte Einheit von μ ist 1/cm
        print(tools.fmt_compare_to_ref(material['μ_lit'] * ureg.dimensionless, μ))


    mu_mat = np.linalg.inv(A.T @ A) @ A.T @ y

    print(f'{material["name"]}:')
    # print(f'\t{mu_mat}')

    μ = mu_mat.mean()

    print(f"{μ=}")

    if 'lit' in material:
        print(tools.fmt_compare_to_ref(material['lit'] * ureg.dimensionless, μ))
