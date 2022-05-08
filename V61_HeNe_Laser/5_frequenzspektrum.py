import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
import rich
from rich.console import Console
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

L = ureg('162.9 cm')
# L = ureg('50 cm') # TODO: TEST!

f, I = np.genfromtxt('dat/5_frequenzspektrum.csv', delimiter=',', skip_header=1, unpack=True)
f *= ureg('MHz')
I *= ureg('dBm')

Δf_theo = (ureg.c / (2 * L)).to('MHz')

Δf_all = np.diff(f)
Δf = tools.ufloat_from_list(Δf_all)
print(tools.fmt_compare_to_ref(Δf, Δf_theo))
