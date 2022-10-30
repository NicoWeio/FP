# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import uncertainties.unumpy as unp
from rich.console import Console

# import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

T = ureg('1 s')

# █ Daten einlesen
α, N = np.genfromtxt("data/3_rockingscan1.txt", unpack=True)  # skip_header=1

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

α *= ureg.degree
N *= ureg.dimensionless

I = N / T


# █ Tabelle generieren
# generate_table.generate_table_pint(
#     'build/tab/1_koinzidenz.tex',
#     (r't_\text{diff}', ureg.degree, Δt),
#     ('I', ureg.second**-1, tools.nominal_values(I)),  # TODO: make ufloats work with tables (again)
# )


I_thresh = tools.nominal_value(20 * min(I))  # TODO: optimize
# find the first index where I > I_thresh and the last index where I > I_thresh
i_min = next(i for i, x in enumerate(I) if x > I_thresh)
i_max = next(i for i, x in reversed(list(enumerate(I))) if x > I_thresh)

α_mean = np.mean(np.abs(α[[i_min, i_max]]))
print(f'α_mean = {α_mean}')

# █ alternative Berechnung
d_0 = ureg('0.28 mm')  # Strahlbreite (@Mampfzwerg)
D = ureg('20 mm') # Probendicke (@Mampfzwerg)
α_alt = np.arcsin(d_0 / D)
print(f'α_alt = {α_alt.to("°")}')


# █ Plot
# α_linspace = tools.linspace(*tools.bounds(α), 1000)

plt.figure()
with tools.plot_context(plt, '°', '1/s', r"\alpha", "I") as plt2:
    plt2.plot(α, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

    plt.axhline(I_thresh, color='C1', alpha=0.5, zorder=0, label="Schwellwert")
    plt.axvline(α[i_min], color='C2', alpha=0.5, zorder=0, label="Geometriewinkel")
    plt.axvline(α[i_max], color='C2', alpha=0.5, zorder=0)

plt.yscale('log')  # TODO
plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig("build/plt/3_rockingscan1.pdf")
plt.show()
