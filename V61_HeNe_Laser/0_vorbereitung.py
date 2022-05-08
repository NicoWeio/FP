import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()


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
    {
        'r1': np.infty * ureg.mm,
        'r2': np.infty * ureg.mm,
    },
]

L = np.arange(0, 3 + 0.01, 0.01) * ureg.m


for num, mirror in enumerate(MIRRORS):
    r1, r2 = mirror['r1'], mirror['r2']
    g1 = calc_g(L, r1)
    g2 = calc_g(L, r2)

    def fmt_r(r):
        if r == np.infty * ureg.mm:
            return '∞'
        else:
            # return f"{r.to('m').m:.1f}"
            return r'\SI{' f"{r.to('m').m:.1f}" r'}{\meter}'

    with tools.plot_context(plt, 'm', 'dimensionless', 'L', 'g_1 \cdot g_2') as plt2:
        plt2.plot(L, g1*g2, label=f"$r_1 = {fmt_r(r1)}$, $r_2 = {fmt_r(r2)}$")


plt.fill_between(L, 0, 1, color='grey', alpha=0.5, label='stabil')
plt.xlim(0, 3)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('build/plt/2_stabilitaetsbedingung_theorie.pdf')
# plt.show()
