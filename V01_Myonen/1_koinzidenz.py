# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import uncertainties.unumpy as unp
from rich.console import Console

import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

T = ureg('20 s')

# █ Daten einlesen
Δt, N = np.genfromtxt("1_koinzidenz.csv", unpack=True, delimiter=",", skip_header=1)

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

Δt *= ureg.nanosecond
N *= ureg.dimensionless

I = N / T


# █ Tabelle generieren
# TODO

# █ lineare Regression
# m, b = tools.linregress(Δt, I)
# print(f"{(m,b)=}")

Δt_peak = Δt[np.argmax(I)]

m_l, b_l = tools.linregress(Δt[Δt < Δt_peak], tools.nominal_values(I[Δt < Δt_peak]))
m_r, b_r = tools.linregress(Δt[Δt > Δt_peak], tools.nominal_values(I[Δt > Δt_peak]))

# █ Plot
# Δt_linspace = tools.linspace(*tools.bounds(Δt))
Δt_linspace_l = tools.linspace(*tools.bounds(Δt[Δt < Δt_peak]))
Δt_linspace_r = tools.linspace(*tools.bounds(Δt[Δt > Δt_peak]))

plt.figure()
with tools.plot_context(plt, 'ns', '1/s', "Δt", "I") as plt2:  # TODO
    plt2.plot(Δt, I, fmt='x', zorder=5, label="Messwerte")
    plt2.plot(Δt_linspace_l, tools.nominal_values(
        m_l * Δt_linspace_l + b_l), label="Ausgleichsgerade links")
    plt2.plot(Δt_linspace_r, tools.nominal_values(
        m_r * Δt_linspace_r + b_r), label="Ausgleichsgerade rechts")

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/1_koinzidenz.pdf")
plt.show()
