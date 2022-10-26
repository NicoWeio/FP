# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import uncertainties.unumpy as unp
from rich.console import Console

import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

T = ureg('100000 s')  # TODO
t_per_channel = ureg('0.0217 µs')

# █ Daten einlesen
N = np.genfromtxt("3_lebensdauer.csv", unpack=True, skip_header=1)  # delimiter=",",
channel = np.arange(len(N))

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

channel *= ureg.dimensionless
N *= ureg.dimensionless

# I = N / T

# Kanal → Zeit
t = channel * t_per_channel

# █ Tabelle generieren
generate_table.generate_table_pint(
    'build/tab/3_lebensdauer.tex',
    ('channel', ureg.dimensionless, channel, 0),
    ('N', ureg.dimensionless, N, 0),
)


# █ Fit


def fit_fn(t, N_0, λ, U_2):
    return N_0 * np.exp(-λ * t) + U_2


# N_0, λ, U_2 = tools.pint_curve_fit(
#     fit_fn, t, tools.nominal_values(N), (ureg.dimensionless, ureg('1/s'), ureg.dimensionless),
#     p0=(tools.nominal_value(max(N)), 0.5E6 / ureg.s, ureg('0 dimensionless')),
# )  # TODO: respect Poisson error
# print(f"{(N_0, λ, U_2)=}")

N_0, λ, U_2 = (tools.nominal_value(max(N)), 0.9E6 / ureg.s, ureg('0 dimensionless'))


# █ Plot
t_linspace = tools.linspace(*tools.bounds(t))


plt.figure()
with tools.plot_context(plt, 'µs', 'dimensionless', "t", "N") as plt2:  # TODO
    plt2.plot(t, N, fmt='x', zorder=5, label="Messwerte")
    # plt2.plot(t_linspace, fit_fn(t_linspace, tools.nominal_value(N_0), tools.nominal_value(λ), tools.nominal_value(U_2)), label="Fit")
    plt2.plot(t_linspace, fit_fn(t_linspace, N_0, λ, U_2), label="Fit")
plt.yscale('symlog')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/3_lebensdauer.pdf")
plt.show()
