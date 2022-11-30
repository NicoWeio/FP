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

# TODO: Übernommen aus 2_mca.py
t_per_channel = ureg('0.0217 µs')
t_offset = ureg('0.13913043478260878 µs')

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
t = channel * t_per_channel + t_offset

# █ Tabelle generieren
print("Tabelle generieren…")
# generate_table.generate_table_pint(
#     'build/tab/3_lebensdauer.tex',
#     (r'\text{Kanal}', ureg.dimensionless, channel, 0),
#     ('N', ureg.dimensionless, tools.nominal_values(N), 0),  # TODO: make ufloats work with tables (again)
# )


# █ Untergrund & Fit
def fit_fn(t, N_0, τ, U):
    return N_0 * np.exp(-t / τ) + U


first_nonzero = np.argmax(N)
last_nonzero = 448  # TODO: automatisieren
num_channels_support = last_nonzero - first_nonzero + 1
print(f"Angesprochener Kanal-Bereich: {first_nonzero} bis {last_nonzero} ({last_nonzero - first_nonzero} Kanäle)")

# fit_mask = (t < ureg('4 µs')) & (N > 10)
fit_mask = N > 0
print(f"Kanäle mit Ereignissen: {np.sum(fit_mask)} (von {len(fit_mask)})")

underground_mask = fit_mask.copy()
U_1_all_channels = I_start * T_search * np.exp((I_start * T_search).to('dimensionless').n) * N_start
U_1 = U_1_all_channels / np.sum(fit_mask)
print(f"Untergrund 1 (berechnet): {U_1.to('dimensionless'):.2f}")
print(f"→ für alle Kanäle: {U_1_all_channels.to('dimensionless'):.2f}")

# VARIANT: Den Untergrund gleichmäßig von allen Nicht-Null-Kanälen abziehen (!)
N[fit_mask] -= U_1

N_0, τ, U_2 = tools.pint_curve_fit(
    fit_fn,
    # t, tools.nominal_values(N),
    t[fit_mask], N[fit_mask],
    (ureg.dimensionless, ureg.microsecond, ureg.dimensionless),
    p0=(tools.nominal_value(max(N)), ureg('2 µs'), ureg('0 dimensionless')),
)
print(f"{(N_0, τ, U_2)=}")
print(f"Untergrund 2 (Fit): {U_2:.2f}")
print(f"→ für alle Kanäle: {(U_2 * num_channels_support).to('dimensionless'):.2f}")

# N_0, τ, U_2 = (tools.nominal_value(max(N)), ureg('2 µs'), ureg('3 dimensionless'))

τ_lit = ureg('2.1969811 µs')
print(tools.fmt_compare_to_ref(τ, τ_lit, name='τ'))


# █ Plot
t_linspace = tools.linspace(*tools.bounds(t), 1_000)

for yscale in ['linear', 'log']:
    plt.figure()
    with tools.plot_context(plt, 'µs', 'dimensionless', "t", "N") as plt2:
        plt2.plot(
            t[fit_mask], N[fit_mask], fmt='.', zorder=5,
            elinewidth=0.5, markeredgewidth=0.5,
            label="Messwerte",
        )
        plt2.plot(
            t[~fit_mask], N[~fit_mask], '.r', zorder=5,
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
    # if yscale == 'log':
    #     plt.ylim(bottom=np.min(N[N>0])/2)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/3_lebensdauer_{yscale}.pdf")
    # plt.show()
