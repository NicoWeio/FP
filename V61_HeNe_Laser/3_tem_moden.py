import matplotlib.pyplot as plt
import numpy as np
import pint
import rich
from rich.console import Console
import tools
from generate_table import generate_table_pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()


def calc_I_TEM00(
    r,  # Parameter
    I_max, I0, r0, ω,  # Konstanten
):
    """TEM₀₀: Theoretische Intensität in Abhängigkeit vom Abstand r von der Modenmitte."""
    return I_max * np.exp(- (r - r0)**2 / (2 * ω**2)) + I0


def calc_I_TEM01(
    r,  # Parameter
    I_max, I0, r0, ω,  # Konstanten
):
    # Hermite-Polynome (zum Quadrat), nicht Laguerre-Polynome
    """TEM₀₁: Theoretische Intensität in Abhängigkeit vom Abstand r von der Modenmitte."""
    return I_max * (4 * (r - r0)**2 / ω**2) * np.exp(- (r - r0)**2 / (2 * ω**2)) + I0


DATA = [
    {
        'path': 'dat/3_tem_00.csv',
        'name': 'tem_00',
        'calc_I': calc_I_TEM00,
    },
    {
        'path': 'dat/3_tem_01.csv',
        'name': 'tem_01',
        'calc_I': calc_I_TEM01,
    },
]

for setup in DATA:
    console.rule(setup['name'])
    # Daten einlesen
    r, I = np.genfromtxt(setup['path'], delimiter=',', skip_header=1, unpack=True)
    r *= ureg('mm')
    I *= ureg.microwatt  # TODO
    calc_I = setup['calc_I']

    # Tabelle erzeugen
    generate_table_pint(f"build/tab/3_{setup['name']}.tex", ('r', ureg.mm, r, 0), ('I', ureg.microwatt, I, 1))

    # Fit berechnen
    params = tools.pint_curve_fit(calc_I, r, I, (ureg.microwatt, ureg.microwatt, ureg.mm, ureg.mm))
    print(f"params: {params}")
    r_linspace = tools.linspace(*tools.bounds(r), 500)

    # Plotten
    plt.figure()
    with tools.plot_context(plt, 'mm', 'microwatt', 'r', 'I') as plt2:
        plt2.plot(r, I, 'x', zorder=5, label='Messwerte')
        plt2.plot(r_linspace, calc_I(r_linspace, *[tools.nominal_value(p) for p in params]), label='Fit')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/3_{setup['name']}.pdf")
    # plt.show()
