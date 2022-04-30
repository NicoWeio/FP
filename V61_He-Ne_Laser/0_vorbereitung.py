import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()




# with tools.plot_context(plt, 'A', 'dimensionless', 'I', 'dimensionless') as plt2:
    # plt2.plot(I, B, 'x', zorder=5, label='Messwerte')
    # plt2.plot(I, tools.nominal_values(slope*I+intercept), label='Regressionsgerade')


def calc_g(L, r):
    """Resonatorparameter"""
    return 1 - L/r

MIRRORS = [
    {
        'r1': ureg('1400 mm'),
        'r2': ureg('1400 mm'),
    },
    {
        'r1': ureg('1000 mm'),
        'r2': ureg('1400 mm'),
    },
    {
        'r1': ureg('1400 mm'),
        'r2': np.infty * ureg.mm,
    },
]

L = np.arange(-0.5, 3, 0.01) * ureg.m


for num, mirror in enumerate(MIRRORS):
    r1, r2 = mirror['r1'], mirror['r2']
    g1 = calc_g(L, r1)
    g2 = calc_g(L, r2)

    plt.plot(L, g1*g2, label=f"Anordnung {num} ($r_1 = {r1.to('m').m:.1f}, r_2 = {r2.to('m').m:.1f}$)")

    # TODO: Bereich farblich hervorheben
    plt.axhline(0, color='grey')
    plt.axhline(1, color='grey')

    plt.grid()
    plt.legend()
    plt.tight_layout()
    # plt.savefig('build/plt/….pdf')

plt.show()
