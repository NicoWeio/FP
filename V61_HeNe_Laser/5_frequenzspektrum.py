import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
import rich
from rich.console import Console
from generate_table import generate_table_pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

L = ureg('162.9 cm')
# L = ureg('50 cm') # TODO: TEST!

f, I = np.genfromtxt('dat/5_frequenzspektrum.csv', delimiter=',', skip_header=1, unpack=True)
f *= ureg('MHz')
I *= ureg('dBm')

# Tabelle erzeugen
generate_table_pint(
    f'build/tab/5_frequenzspektrum.tex',
    ('f', ureg.MHz, f, 0),
    ('I', ureg.dBm, I, 1),
    ('I', ureg.microwatt, I, 1),
)


Δf_theo = (ureg.c / (2 * L)).to('MHz')

f_theo = Δf_theo * np.arange(1, (max(f) // Δf_theo).m + 1)
print(f'f_mess = {f:.2f}')
print(f'f_theo = {f_theo:.2f}')

for single_f in f:
    best_n = np.argmin(np.abs(f_theo - single_f))  # off-by-one!
    print(f'Mess: {single_f} ↔ Theo: {f_theo[best_n]} (n={best_n+1}) → Abweichung: {(f_theo[best_n] - single_f):.2f}')


# Ansatz siehe Mampfzwerg. Keine schönen Ergebnisse…
# Δf_all = np.diff(f)
# print(f'Δf_all = {Δf_all}')
# Δf = tools.ufloat_from_list(Δf_all)
# print(tools.fmt_compare_to_ref(Δf, Δf_theo))
