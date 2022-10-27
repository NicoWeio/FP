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

T = ureg('159914 s')
print(f"Messzeit: {T.to('hours'):.2f} bzw. {T.to('days'):.2f}")

t_per_channel = ureg('0.0217 µs')  # TODO: Übernommen aus 2_mca.py

# █ Daten einlesen
N = np.genfromtxt("3_lebensdauer.csv", unpack=True, skip_header=1)
channel = np.arange(len(N))

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

channel *= ureg.dimensionless
N *= ureg.dimensionless

# I = N / T

# Kanal → Zeit
t = channel * t_per_channel

# █ Tabelle generieren
print("Tabelle generieren…")
generate_table.generate_table_pint(
    'build/tab/3_lebensdauer.tex',
    (r'\text{Kanal}', ureg.dimensionless, channel, 0),
    ('N', ureg.dimensionless, tools.nominal_values(N), 0),  # TODO: make ufloats work with tables (again)
)


# █ Fit


def fit_fn(t, N_0, τ, U_2):
    return N_0 * np.exp(-t / τ) + U_2


# fit_mask = (t < ureg('4 µs')) & (N > 10)
fit_mask = N > 0

N_0, τ, U_2 = tools.pint_curve_fit(
    fit_fn,
    # t, tools.nominal_values(N),
    t[fit_mask], N[fit_mask],
    (ureg.dimensionless, ureg.microsecond, ureg.dimensionless),
    p0=(tools.nominal_value(max(N)), ureg('2 µs'), ureg('3 dimensionless')),
)  # TODO: respect Poisson error
print(f"{(N_0, τ, U_2)=}")

# N_0, τ, U_2 = (tools.nominal_value(max(N)), ureg('2 µs'), ureg('3 dimensionless'))


# █ Plot
t_linspace = tools.linspace(*tools.bounds(t))


plt.figure()
with tools.plot_context(plt, 'µs', 'dimensionless', "t", "N") as plt2:  # TODO
    plt2.plot(
        t[fit_mask], N[fit_mask], fmt='x', zorder=5,
        elinewidth=0.5, markeredgewidth=0.5,
        label="Messwerte",
    )
    plt2.plot(
        t[~fit_mask], N[~fit_mask], 'xr', zorder=5,
        markeredgewidth=0.5,
        label="Messwerte (nicht berücksichtigt)",
    )
    plt2.plot(
        t_linspace,
        fit_fn(
            t_linspace,
            tools.nominal_value(N_0),
            tools.nominal_value(τ),
            tools.nominal_value(U_2)
        ),
        zorder=10,
        label="Fit",
    )
# plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/3_lebensdauer.pdf")
plt.show()
