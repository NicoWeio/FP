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

with tools.plot_context(plt, 'cm', 'dimensionless', 'L', 'I') as plt2:
    for setup in DATA:
        L, I = np.genfromtxt(setup['path'], delimiter=',', skip_header=1, unpack=True)
        L *= ureg('cm')
        I *= ureg.dimensionless  # TODO

        plt2.plot(L, I, '--o', zorder=5, label=setup['label'])

plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig('build/plt/….pdf')
plt.show()
