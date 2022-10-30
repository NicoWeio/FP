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

T = ureg('5 s')

# █ Daten einlesen
α, N = np.genfromtxt("data/4_reflektivitätsscan.txt", unpack=True)  # skip_header=1

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

λ = ureg('1,54 Å')  # ? (@Mampfzwerg)
q = 4 * np.pi / λ * np.sin(np.pi / 180 * α)  # TODO: Blind übernommen aus @Mampfzwerg


# █ Plot
# α_linspace = tools.linspace(*tools.bounds(α), 1000)

plt.figure()
# TODO: Doppelachse mit Intensität und Reflektivität?
with tools.plot_context(plt, '1/m', '1/s', "q", "I") as plt2:
    plt2.plot(q, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig("build/plt/4_reflektivitätsscan.pdf")
plt.show()
