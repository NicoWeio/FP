import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

I, B = np.genfromtxt('1_magnet.csv', delimiter=',', skip_header=1, unpack=True)
I *= ureg('A')
B *= ureg('mT')

slope, intercept = tools.linregress(I, B)
print(f"{slope=}, {intercept=}")

with tools.plot_context(plt, 'A', 'mT', 'I', 'B') as plt2:
    plt2.plot(I, B, 'x', zorder=5, label='Messwerte')
    plt2.plot(I, tools.nominal_values(slope*I+intercept), label='Regressionsgerade')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('build/plt/1_magnet.pdf')
plt.plot()
