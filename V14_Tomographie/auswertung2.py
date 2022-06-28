import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from numpy.linalg import inv
import pandas as pd
import pint
import rich
from rich.console import Console
import generate_table
import tools
ureg = pint.UnitRegistry()
console = Console()

d_Einzelwürfel = ureg('1 cm')

μ_LIT = {
    # TODO: kopiert aus calc_mu.py
    'Al': ureg('0.2007 / cm'),
    'Pb': ureg('1.1737 / cm'),
    'Fe': ureg('0.5704 / cm'),
    'Messing': ureg('0.6088 / cm'),
    'Delrin': ureg('0.1176 / cm'),
}


def A_row_from_indices(indices):
    """
    Gibt eine Zeile der Matrix A zurück.
    Die Einträge stehen dabei für die je Elementarwürfel (µ_i) passierte Wegstrecke.

    Beispiel:
    '2/4' → [0, √2·d, 0, √2·d, 0, 0, …]
    """
    d_Einzelwürfel = 1  # NOTE: Einheiten + Matrizen = 🗲
    indices = indices.strip()  # TODO Zwischenlösung
    row = np.zeros(9)  # 9 = Anz. Würfel
    if '/' in indices:
        assert '|' not in indices
        # → diagonal
        for i in map(int, indices.split('/')):
            row[i-1] = np.sqrt(2) * d_Einzelwürfel
    elif '|' in indices:
        assert '/' not in indices
        # → parallel
        for i in map(int, indices.split('|')):
            row[i-1] = d_Einzelwürfel
    return row


def d_row_from_indices(indices):
    # NOTE: `row` steht hier für ein einzelnes Element…
    """
    Gibt die gesamte Wegstrecke im Würfel zurück.
    """
    return A_row_from_indices(indices).sum() * ureg.cm  # TODO: Einheiten hübscher…


def A_from_indices(all_indices):
    return np.row_stack(list(map(A_row_from_indices, all_indices)))


def d_from_indices(all_indices):
    return tools.pintify(list(map(d_row_from_indices, all_indices)))


def I_0_from_indices(all_indices):
    return tools.pintify(list(map(I_0_row_from_indices, all_indices)))


def I_0_row_from_indices(indices):
    """
    Gibt die zu den Indizes gehörige Nullmessung zurück.

    Beispiele:
    • '4|5|6' → I_0_Parallel
    • '2/4' → I_0_Nebendiagonale
    """
    if '/' in indices:
        if len(indices.split('/')) == 3:
            return I_0_hauptdiag
        else:
            assert len(indices.split('/')) == 2
            return I_0_nebendiag
    else:
        assert '|' in indices
        return I_0_parallel


@ureg.wraps(ureg.s**(-1), [ureg.dimensionless, ureg.s])
def calc_I(N, T):
    """
    Zählrate aus Anzahl N und Zeit T.
    """
    I = unp.uarray(N / T, np.sqrt(N)/T) / ureg.s
    return I


def analyze_homogen(dat, I_0_getter):
    I_0 = I_0_getter(dat['indices'])
    d = d_from_indices(dat['indices'])
    μ_vec = (
        unp.log((I_0 / dat['I']).to('dimensionless').m) /
        d
    )

    # print(f'{μ_vec=}')
    µ = µ_vec.mean()  # TODO: SEM
    # print(f'{μ=}')
    return µ


def analyze_inhomogen(dat, I_0_getter):
    A = A_from_indices(dat['indices'])
    # print(A.round(1))

    I_0 = I_0_getter(dat['indices'])
    y = unp.log((I_0 / dat['I']).to('dimensionless').m)
    μ_vec = (np.linalg.inv(A.T @ A) @ A.T @ y)
    μ_vec /= ureg.cm  # TODO: hübscher
    return µ_vec


def get_data(filename):
    dat = pd.read_csv(filename, comment='#')

    N = dat['N'].to_numpy() * ureg.dimensionless
    T = dat['T'].to_numpy() * ureg.s

    I = calc_I(N, T)

    # NOTE: Ohne pint-pandas wird implizit zu einheitenlosen Arrays konvertiert,
    # daher verwenden wir hier ein normales dict.
    return {
        'N': N * ureg.dimensionless,
        'T': T * ureg.s,
        'I': I,
        'indices': dat['indices'],
    }


def get_closest_material(µ):
    diff_tuples = [(abs(µ - µ_lit_single), name) for name, µ_lit_single in µ_LIT.items()]
    diff_tuples.sort()
    return diff_tuples[0][1]


console.rule("Nullmessung")
dat_Nullmessung = get_data('dat/Nullmessung.csv')
assert dat_Nullmessung['I'].check('1/[time]')
I_0_parallel, I_0_hauptdiag, I_0_nebendiag = dat_Nullmessung['I']
print(f"I_0 = {tools.nominal_values(dat_Nullmessung['I']):.2f}")

WÜRFEL = [
    {
        'num': 2,
        'material': 'Fe',  # Mampfzwerg, SanjoR
        'homogen': True,
    },
    {
        'num': 3,
        'material': 'Delrin',  # Insider-Tipp: eigentlich Holz
        'homogen': True,
    },
    # {
    #     'num': 4,
    #     # 'material': [1, 2, 3, 4, 5, 6, 7, 8, 9],  # TODO
    #     'homogen': False,
    # },
]

for würfel in WÜRFEL:
    console.rule(f"Würfel {würfel['num']}")
    dat = get_data(f'dat/Würfel{würfel["num"]}.csv')
    if würfel['homogen']:
        µ = analyze_homogen(dat, I_0_from_indices)
    else:
        µ = analyze_inhomogen(dat, I_0_from_indices)
    print(f"μ = {μ:.3f}")
    print(f"rel. Unsicherheit: {µ.s/µ.n:.1%}")

    mat = get_closest_material(µ)
    print(f"Best fit: {mat}")
    if mat == würfel['material']:
        print("→ wie erwartet :)")
    else:
        print(f"→ sollte {würfel['material']} sein :(")
    print(
        "Abweichung µ vs. µ_lit (best fit):\n" +
        tools.fmt_compare_to_ref(µ, µ_LIT[mat])
    )
    print(
        f"Abweichung best-fit vs. µ_lit des tatsächlichen Materials ({würfel['material']}):\n" +
        tools.fmt_compare_to_ref(µ, µ_LIT[würfel['material']])
    )

    generate_table.generate_table_pint(
        f'build/tab/wuerfel{würfel["num"]}.tex',
        ('Projektion', None, dat['indices']),
        ('I', ureg.second**(-1), dat['I']),
        # ('µ', ureg.centimeter**(-1), µ_vec), # TODO
    )


console.rule("Würfel 4")
dat = get_data(f'dat/Würfel4.csv')
µ = analyze_inhomogen(dat, I_0_from_indices)
for y in range(3):
    print('\t'.join([f"{x.n:.2f}" for x in µ[3*y:3*y+3]]))


x, y = np.arange(3), np.arange(3)
µ_plt = unp.nominal_values(µ).reshape(3, 3)
plt.pcolormesh(x, y, µ_plt)
plt.gca().invert_yaxis()
plt.gca().set_aspect('equal')
plt.xticks([])
plt.yticks([])
plt.colorbar()
for i, (y, x) in enumerate(np.ndindex(3, 3), 1):
    plt.text(x, y, r'$\mu_' f'{i}$', ha='center', va='center')
plt.savefig('build/plt/wuerfel4.pdf')
# plt.show()
