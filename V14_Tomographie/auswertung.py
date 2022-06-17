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


def kleinsteQuadrate(y, W, A):
    temp = np.dot(np.linalg.inv(np.dot(A.T, np.dot(W, A))), A.T)
    a = np.dot(temp, np.dot(W, y))
    a_err = np.linalg.inv(np.dot(A.T, np.dot(W, A)))
    return a, np.sqrt(np.diag(a_err))


A = np.matrix([[0, 0, 0, 0, 0, np.sqrt(2), 0, np.sqrt(2), 0],
               [0, 0, np.sqrt(2), 0, np.sqrt(2), 0, np.sqrt(2), 0, 0],
               [0, np.sqrt(2), 0, np.sqrt(2), 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 1, 1, 1],
               [0, 0, 0, 1, 1, 1, 0, 0, 0],
               [1, 1, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, np.sqrt(2), 0, 0, 0, np.sqrt(2), 0],
               [np.sqrt(2), 0, 0, 0, np.sqrt(2), 0, 0, 0, np.sqrt(2)],
               [0, np.sqrt(2), 0, 0, 0, np.sqrt(2), 0, 0, 0],
               [1, 0, 0, 1, 0, 0, 1, 0, 0],
               [0, 1, 0, 0, 1, 0, 0, 1, 0],
               [0, 0, 1, 0, 0, 1, 0, 0, 1]])

A_t = np.transpose(A)  # TODO: besser direkt A.T

MATERIALS = [
    {
        'name': 'Alu',
        'lit': 0.202,
    },
    {
        'name': 'Blei',
        'lit': 1.25,
    },
    {
        'name': 'Unb',
    },
    {
        'name': 'Luft',
    },
]


for material in MATERIALS:
    console.rule(material['name'])
    N, t = np.genfromtxt(f'dat/{material["name"]}.csv', delimiter=',', skip_header=1, unpack=True)
    I = unp.uarray(N / t, np.sqrt(N)/t) / ureg.s
    N *= ureg.dimensionless  # TODO
    t *= ureg.s

    I_0 = ureg('184/s')  # TODO: Nullmessung pro Material
    y = unp.log((I_0 / I).to('dimensionless').m)

    mu_mat = np.linalg.inv(A.T @ A) @ A.T @ y

    print(f'{material["name"]}:')
    # print(f'\t{mu_mat}')

    μ = mu_mat.mean()

    print(f"{μ=}")

    if 'lit' in material:
        print(tools.fmt_compare_to_ref(material['lit'] * ureg.dimensionless, μ))
