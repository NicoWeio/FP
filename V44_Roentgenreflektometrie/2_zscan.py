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
z, N = np.genfromtxt("data/2_zscan1.txt", unpack=True)  # skip_header=1

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

z *= ureg.mm
N *= ureg.dimensionless

I = N / T


# █ Tabelle generieren
# generate_table.generate_table_pint(
#     'build/tab/1_koinzidenz.tex',
#     (r't_\text{diff}', ureg.degree, Δt),
#     ('I', ureg.second**-1, tools.nominal_values(I)),  # TODO: make ufloats work with tables (again)
# )


# TODO: Auto-detect?
flank_bound_indices = (17, 23)
flank_bounds = tuple(z[list(flank_bound_indices)])


# █ Plot
# z_linspace = tools.linspace(*tools.bounds(z), 1000)

plt.figure()
with tools.plot_context(plt, 'mm', '1/s', "z", "I") as plt2:
    plt2.plot(z, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

    plt.axvspan(*flank_bounds, color='C1', alpha=0.5, zorder=0, label="Strahlbreite")

plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig("build/plt/1_detektorscan.pdf")
plt.show()
