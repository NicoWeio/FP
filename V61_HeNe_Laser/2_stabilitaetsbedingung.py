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
        'max_theo': 1400 * ureg.mm,
    },
    {
        'path': 'dat/2_stabilitaetsbedingung_konkav_konkav.csv',
        'label': 'konkav + konkav',
        'max_theo': (1400+1400) * ureg.mm,
    },
]


def sort_xy(x, y):  # TODO: Auslagern in Tools
    tuples = list(zip(x, y))
    tuples.sort(key=lambda t: t[0])
    x = tools.pintify([t[0] for t in tuples])
    y = tools.pintify([t[1] for t in tuples])
    return x, y


with tools.plot_context(plt, 'cm', 'mW', 'L', 'I') as plt2:
    for index, setup in enumerate(DATA):
        L, I = np.genfromtxt(setup['path'], delimiter=',', skip_header=1, unpack=True)
        L *= ureg('cm')
        I *= ureg.mW

        L, I = sort_xy(L, I)

        plt2.plot(L, I, '--o', zorder=5, label=setup['label'])
        plt.axvline(x=setup['max_theo'].to('cm'), linestyle='-', color=f'C{index}', label='↪ theoretisches Maximum')

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('build/plt/2_stabilitaetsbedingung.pdf')
# plt.show()
