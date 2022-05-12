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

f, I = np.genfromtxt('dat/5_frequenzspektrum.csv', delimiter=',', skip_header=1, unpack=True)
f *= ureg('MHz')
I *= ureg('dBm')


Δf_theo = (ureg.c / (2 * L)).to('MHz')
print(f'Δf_theo = {Δf_theo:.2f}')

f_theo = Δf_theo * np.arange(0, (max(f) // Δf_theo).m + 1)
print(f'f_mess = {f:.2f}')
print(f'f_theo = {f_theo:.2f}')


data = []
for single_f in f:
    best_n = np.argmin(np.abs(f_theo - single_f))
    best_f_theo = f_theo[best_n]
    deviation = f_theo[best_n] - single_f
    data.append((best_n, best_f_theo, deviation))
    print(f'Mess: {single_f:.2f} ↔ Theo: {best_f_theo:.2f} (n={best_n}) → Abweichung: {deviation:.2f}, rel. Abw.: {(deviation / best_f_theo).m:.2%}')

best_n, best_f_theo, deviation = zip(*data)

# Tabelle erzeugen
generate_table_pint(
    f'build/tab/5_frequenzspektrum.tex',
    ('n', ureg.dimensionless, best_n * ureg.dimensionless, 0),
    ('f', ureg.MHz, f, 0),
    ('I', ureg.dBm, I, 1),
    # ('I', ureg.microwatt, I, 1),
    (r'f_\text{theo}', ureg.MHz, best_f_theo, 1),
    (r'(f_\text{theo} - f)', ureg.MHz, deviation, 1),
)


# Ansatz siehe Mampfzwerg. Keine schönen Ergebnisse…
# Δf_all = np.diff(f)
# print(f'Δf_all = {Δf_all}')
# Δf = tools.ufloat_from_list(Δf_all)
# print(tools.fmt_compare_to_ref(Δf, Δf_theo))
