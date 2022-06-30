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
    'Holz': ureg('0.0612 / cm'),
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

    return µ_vec


def analyze_inhomogen(dat, I_0_getter):
    """
    Berechne die µ-Werte eines inhomogenen Würfels mittels der Methode der gewichteten kleinsten Quadrate.

    Mehr Infos zum Verfahren: https://en.wikipedia.org/wiki/Weighted_least_squares
    """
    A = A_from_indices(dat['indices'])
    # print(A.round(1))

    I_0 = I_0_getter(dat['indices'])

    # alt ↓
    # y = unp.log((I_0 / dat['I']).to('dimensionless').m)
    # μ_vec = (np.linalg.inv(A.T @ A) @ A.T @ y)

    # neu ↓
    y = unp.log((I_0 / dat['I']).to('dimensionless').m)
    y_var = unp.std_devs(y)**2
    W = np.diag(1 / y_var)
    μ_vec = (np.linalg.inv(A.T @ W @ A) @ A.T @ W @ y)
    µ_vec /= ureg.cm

    # ↓ Das ist erstmal noch kein Vektor, sondern eine Kovarianzmatrix! Wir betrachten die Diagonale.
    µ_vec_var = np.diag(np.linalg.inv(A.T @ W @ A))
    µ_vec_std = np.sqrt(µ_vec_var)

    # ↓ smoke test
    assert np.isclose(unp.std_devs(µ_vec), µ_vec_std, rtol=0.2).all()

    # NOTE: Für den Augenblick sollen die Unsicherheiten, die UNumPy mitgeschleppt hat, genügen.

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


def get_closest_material(µ, µ_map):
    """
    µ sind die Messwerte, denen die Materialien zugeordnet werden sollen.
    µ_map ist ein dict mit Materialnamen als Keys und µ-Werten als Values.
    """
    diff_tuples = [(abs(µ - µ_map_single), name) for name, µ_map_single in µ_map.items()]
    diff_tuples.sort()
    return diff_tuples[0][1]


def indices_string_to_list(indices_string):
    """
    Beispiel: '4|5|6' → [4, 5, 6]
    """
    str_list = indices_string.split('/') if '/' in indices_string else indices_string.split('|')
    int_list = [int(s) for s in str_list]
    return int_list


def visualize_indices(indices_list):
    NUMCELLS = 3

    def i_to_coords(i):
        """Inputs and outputs 0-indexed values."""
        return (i % NUMCELLS, i // NUMCELLS)

    assert i_to_coords(1-1) == (0, 0)
    assert i_to_coords(9-1) == (2, 2)

    TEMPLATE_OUTER = r'''
    \begin{tikzpicture}
    [
        box/.style={rectangle, draw=gray!40, inner sep=0pt },
    ]
    \newcommand*{\numcells}{3}%
    \newcommand*{\cellsize}{(2ex/\numcells)}%

\foreach \x in {1,...,\numcells}{
    \foreach \y in {1,...,\numcells}
        \node[box, minimum size=\cellsize] at ({\x*\cellsize},{-\y*\cellsize}){};
}

{FILL}

\end{tikzpicture}'''

    TEMPLATE_FILL = r'\node[box,minimum size=\cellsize,fill=red] at ({#1*\cellsize},{-#2*\cellsize}){};'

    fill_positions = [i_to_coords(int(i)-1) for i in indices_list]

    fill = '\n'.join(TEMPLATE_FILL.replace('#1', str(x+1)).replace('#2', str(y+1)) for x, y in fill_positions)

    return TEMPLATE_OUTER.replace('{FILL}', fill)


console.rule("Nullmessung")
dat_Nullmessung = get_data('dat/Würfel1.csv')
assert dat_Nullmessung['I'].check('1/[time]')
I_0_parallel, I_0_hauptdiag, I_0_nebendiag = dat_Nullmessung['I']
print(f"I_0 = {tools.nominal_values(dat_Nullmessung['I']):.2f}")

generate_table.generate_table_pint(
    f'build/tab/wuerfel1.tex',
    ('Projektion', None, dat_Nullmessung['indices']),
    ('', None, [visualize_indices(indices_string_to_list(i)) for i in dat_Nullmessung['indices']]),
    ('I_0', ureg.second**(-1), dat_Nullmessung['I']),
    # kein µ hier; dafür brauchen wir ja gerade die Nullmessung
)


WÜRFEL = [
    {
        'num': 2,
        'material': 'Fe',  # Mampfzwerg, SanjoR
    },
    {
        'num': 3,
        # 'material': 'Delrin',  # Insider-Tipp: eigentlich Holz
        'material': 'Holz',
    },
]

µ_mess_dict = {}

for würfel in WÜRFEL:
    console.rule(f"Würfel {würfel['num']}")
    dat = get_data(f'dat/Würfel{würfel["num"]}.csv')
    µ_vec = analyze_homogen(dat, I_0_from_indices)
    µ = µ_vec.mean()  # TODO: SEM
    print(f"μ = {μ:.3f}")
    print(f"rel. Unsicherheit: {µ.s/µ.n:.1%}")

    µ_mess_dict[f"Würfel {würfel['num']}"] = µ

    mat = get_closest_material(µ, µ_LIT)
    print(f"Best fit: {mat}")
    if mat == würfel['material']:
        print("→ wie erwartet :)")
    else:
        print(f"→ sollte {würfel['material']} sein :(")
    print(
        "Abweichung µ vs. µ_lit (best fit):\n" +
        tools.fmt_compare_to_ref(µ, µ_LIT[mat])
    )
    # print(
    #     f"Abweichung best-fit vs. µ_lit des tatsächlichen Materials ({würfel['material']}):\n" +
    #     tools.fmt_compare_to_ref(µ, µ_LIT[würfel['material']])
    # )

    generate_table.generate_table_pint(
        f'build/tab/wuerfel{würfel["num"]}.tex',
        ('Projektion', None, dat['indices']),
        ('', None, [visualize_indices(indices_string_to_list(i)) for i in dat['indices']]),
        ('I', ureg.second**(-1), dat['I']),
        (r'\mu', ureg.centimeter**(-1), µ_vec, 3),
    )


console.rule("Würfel 4")
dat = get_data(f'dat/Würfel4.csv')
µ_vec = analyze_inhomogen(dat, I_0_from_indices)
for y in range(3):
    print('\t'.join([f"{x.n:.2f}" for x in µ_vec[3*y:3*y+3]]))

mat_closest_lit_vec = [get_closest_material(µ, µ_LIT) for µ in µ_vec]
µ_closest_lit_vec = tools.pintify([µ_LIT[mat] for mat in mat_closest_lit_vec])

mat_closest_mess_vec = [get_closest_material(µ, µ_mess_dict) for µ in µ_vec]
µ_closest_mess_vec = tools.pintify([µ_mess_dict[mat] for mat in mat_closest_mess_vec])

generate_table.generate_table_pint(
    f'build/tab/wuerfel4.tex',
    ('Projektion', None, dat['indices']),
    ('', None, [visualize_indices(indices_string_to_list(i)) for i in dat['indices']]),
    ('I', ureg.second**(-1), dat['I']),
    # kein µ hier; für einzelne Projektionen ist µ_vec unterbestimmt
)

ureg.define('percent = 1 / 100')

generate_table.generate_table_pint(
    f'build/tab/wuerfel4_mu.tex',
    ('Index $i$', None, np.arange(len(µ_vec))+1),
    ('', None, [visualize_indices([i]) for i in np.arange(len(µ_vec))+1]),
    (r'\mu', ureg.centimeter**(-1), µ_vec, 3),
    (r'\text{Material}', None, mat_closest_lit_vec),
    (r'\text{Abweichung}', ureg.percent, tools.nominal_values(abs((µ_vec-µ_closest_lit_vec)/µ_closest_lit_vec))),
    (r'\text{Material}', None, mat_closest_mess_vec),
    (r'\text{Abweichung}', ureg.percent, tools.nominal_values(abs((µ_vec-µ_closest_mess_vec)/µ_closest_mess_vec))),
)

# FIXME: Während generate_table keine multicolumns unterstützt, mogeln wir das mal rein…
MULTICOLUMN = r'''
&&& \multicolumn{2}{c}{nächster Literaturwert} & \multicolumn{2}{c}{nächster Messwert} \\
\cmidrule(lr){4-5} \cmidrule(lr){6-7}
'''.strip()
with open(f'build/tab/wuerfel4_mu.tex', 'r') as f:
    lines = f.readlines()
lines = lines[:2] + [MULTICOLUMN] + lines[2:]
with open(f'build/tab/wuerfel4_mu.tex', 'w') as f:
    f.writelines(lines)


# Annahme: Der Würfel besteht aus nur zwei Materialien
# → Bestimmte den jeweiligen durchschnittlichen Absorptionskoefizienten
# TODO: negative µ ignorieren
# material_mask = µ_vec < np.median(µ_vec)
# µ_vec_material_1 = µ_vec[material_mask].mean()
# µ_vec_material_2 = µ_vec[~material_mask].mean()
# print(f"µ_1 = {µ_vec_material_1:.3f}")
# print(f"µ_2 = {µ_vec_material_2:.3f}")


# █ Visualisierung der µ in Würfel 4
x, y = np.arange(3), np.arange(3)
µ_plt = unp.nominal_values(µ_vec).reshape(3, 3)
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
