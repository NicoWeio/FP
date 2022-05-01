import tools
import pint
import numpy as np
import matplotlib.pyplot as plt
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()


def calc_I(
    r,  # Parameter
    I0, r0, ω,  # Konstanten
):
    """Theoretische Intensität in Abhängigkeit vom Abstand r von der Modenmitte."""
    return I0 * np.exp(- (r - r0)**2 / (2 * ω**2))


DATA = [
    {
        'path': 'dat/3_tem_00.csv',
    },
]

with tools.plot_context(plt, 'mm', 'microwatt', 'r', 'I') as plt2:
    for setup in DATA:
        # Daten einlesen
        r, I = np.genfromtxt(setup['path'], delimiter=',', skip_header=1, unpack=True)
        r *= ureg('mm')
        I *= ureg.microwatt  # TODO

        # Fit berechnen
        params = tools.pint_curve_fit(calc_I, r, I, (ureg.microwatt, ureg.mm, ureg.mm))
        print(f"params: {params}")
        nominal_params = [tools.nominal_value(p) for p in params]
        print(f"nominal_params: {nominal_params}")
        r_linspace = tools.linspace(*tools.bounds(r))

        # Plotten
        plt2.plot(r, I, 'x', zorder=5, label='Messwerte')
        plt2.plot(r_linspace, calc_I(r_linspace, *[tools.nominal_value(p) for p in params]), label='Fit')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        # plt.savefig('build/plt/….pdf')
        plt.show()
