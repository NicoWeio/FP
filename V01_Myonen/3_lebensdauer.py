# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import uncertainties.unumpy as unp
from rich.console import Console
from uncertainties import ufloat

import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

T = ureg('159914 s')
print(f"Messzeit: {T.to('hours'):.2f} bzw. {T.to('days'):.2f}")

T_search = ureg('10 µs')

t_per_channel = ureg('0.0217 µs')  # TODO: Übernommen aus 2_mca.py

# █ Daten einlesen
N = np.genfromtxt("3_lebensdauer.csv", unpack=True, skip_header=1)
channel = np.arange(len(N))

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

N_stop = np.sum(N) * ureg.dimensionless
# N_stop = ufloat(N_stop, np.sqrt(N_stop))

# NOTE: Aufgrund eines Fehlers Dritter müssen wir N_start anhand der vorherigen Messungen abschätzen.
# Wir rechnen dazu das Maximum aus der Koinzidenz-Messung auf die längere Messzeit hoch.
I_start = ufloat(19.8, 1) * ureg('1/s')  # aus 1_koinzidenz
N_start = I_start * T

print(f"Start-Events: {N_start:.2f}")
print(f"Stop-Events: {N_stop:.2f}")

channel *= ureg.dimensionless
N *= ureg.dimensionless

# I = N / T

# Kanal → Zeit
t = channel * t_per_channel

# █ Tabelle generieren
print("Tabelle generieren…")
# generate_table.generate_table_pint(
#     'build/tab/3_lebensdauer.tex',
#     (r'\text{Kanal}', ureg.dimensionless, channel, 0),
#     ('N', ureg.dimensionless, tools.nominal_values(N), 0),  # TODO: make ufloats work with tables (again)
# )


# █ Untergrund & Fit
def fit_fn(t, N_0, τ, U_2):
    return N_0 * np.exp(-t / τ) + U_2


# fit_mask = (t < ureg('4 µs')) & (N > 10)
fit_mask = N > 0

underground_mask = fit_mask.copy()
U_1 = I_start * T_search * np.exp((I_start * T_search).to('dimensionless').n) * N_start
U_1_per_channel = U_1 / sum(fit_mask)
print(f"Untergrund 1: {U_1.to('dimensionless'):.2f}")
print(f"→ pro Kanal: {U_1_per_channel.to('dimensionless'):.2f}")
# Den Untergrund gleichmäßig von allen Nicht-Null-Kanälen abziehen (!)
N[fit_mask] -= U_1_per_channel

N_0, τ, U_2 = tools.pint_curve_fit(
    fit_fn,
    # t, tools.nominal_values(N),
    t[fit_mask], N[fit_mask],
    (ureg.dimensionless, ureg.microsecond, ureg.dimensionless),
    p0=(tools.nominal_value(max(N)), ureg('2 µs'), ureg('3 dimensionless')),
)
print(f"{(N_0, τ, U_2)=}")

# N_0, τ, U_2 = (tools.nominal_value(max(N)), ureg('2 µs'), ureg('3 dimensionless'))

τ_lit = ureg('2.1969811 µs')
print(tools.fmt_compare_to_ref(τ, τ_lit))


# █ Plot
t_linspace = tools.linspace(*tools.bounds(t), 1_000)

for yscale in ['linear', 'log']:
    plt.figure()
    with tools.plot_context(plt, 'µs', 'dimensionless', "t", "N") as plt2:
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
    plt.yscale(yscale)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/3_lebensdauer_{yscale}.pdf")
    # plt.show()
