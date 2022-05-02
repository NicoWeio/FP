import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

DATA = [
    {
        'path': 'dat/2_stabilitaetsbedingung_plan_konkav.csv',
        'label': 'plan + konkav',
    },
    {
        'path': 'dat/2_stabilitaetsbedingung_konkav_konkav.csv',
        'label': 'konkav + konkav',
    },
]


def sort_xy(x, y):  # TODO: Auslagern in Tools
    tuples = list(zip(x, y))
    tuples.sort(key=lambda t: t[0])
    x = tools.pintify([t[0] for t in tuples])
    y = tools.pintify([t[1] for t in tuples])
    return x, y


with tools.plot_context(plt, 'cm', 'mW', 'L', 'I') as plt2:
    for setup in DATA:
        L, I = np.genfromtxt(setup['path'], delimiter=',', skip_header=1, unpack=True)
        L *= ureg('cm')
        I *= ureg.mW

        L, I = sort_xy(L, I)

        plt2.plot(L, I, '--o', zorder=5, label=setup['label'])

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('build/plt/2_stabilitaetsbedingung.pdf')
# plt.show()
