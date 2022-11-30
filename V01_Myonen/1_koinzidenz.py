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

T = ureg('20 s')

# █ Daten einlesen
Δt, N = np.genfromtxt("1_koinzidenz.csv", unpack=True, delimiter=",", skip_header=1)
sorting_i = Δt.argsort()
Δt = Δt[sorting_i]
N = N[sorting_i]

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

Δt *= ureg.nanosecond
N *= ureg.dimensionless

I = N / T


# █ Tabelle generieren
generate_table.generate_table_pint(
    'build/tab/1_koinzidenz.tex',
    (r't_\text{diff}', ureg.nanosecond, Δt),
    ('I', ureg.second**-1, tools.nominal_values(I)),  # TODO: make ufloats work with tables (again)
)

Δt_peak = Δt[np.argmax(I)]

m_l, b_l = tools.linregress(Δt[Δt <= Δt_peak], tools.nominal_values(I[Δt <= Δt_peak]))
m_r, b_r = tools.linregress(Δt[Δt >= Δt_peak], tools.nominal_values(I[Δt >= Δt_peak]))

print(f"m_l = {m_l:.2f}")
print(f"b_l = {b_l:.2f}")
print(f"m_r = {m_r:.2f}")
print(f"b_r = {b_r:.2f}")

# widths, width_heights, left_ips, right_ips = sp.signal.peak_widths(tools.nominal_values(I), [np.argmax(I)])

print(f"Peak: {max(I)} bei {Δt[np.argmax(I)]}")
# print(f"Halbwertsbreite: {(width_heights[0] * I.units):.2f}")  # TODO: könnte falsch sein…


# █ Plot
Δt_linspace = tools.linspace(*tools.bounds(Δt))
Δt_linspace_l = tools.linspace(*tools.bounds(Δt[Δt <= Δt_peak]))
Δt_linspace_r = tools.linspace(*tools.bounds(Δt[Δt >= Δt_peak]))

plt.figure()
with tools.plot_context(plt, 'ns', '1/s', r"t_\text{diff}", "I") as plt2:
    plt2.plot(Δt, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

    plt2.plot(Δt_linspace_l, tools.nominal_values(
        m_l * Δt_linspace_l + b_l), label="Ausgleichsgerade links")
    plt2.plot(Δt_linspace_r, tools.nominal_values(
        m_r * Δt_linspace_r + b_r), label="Ausgleichsgerade rechts")

    # left_interp = left_ips[0]
    # left_index = int(left_interp)
    # left_frac = left_interp - left_index
    # left = Δt[left_index] + left_frac * (Δt[left_index+1] - Δt[left_index])

    # right_interp = right_ips[0]
    # right_index = int(right_interp)
    # right_frac = right_interp - right_index
    # right = Δt[right_index] + right_frac * (Δt[right_index+1] - Δt[right_index])

    # plt2.plot(
    #     tools.pintify([left, right]), [width_heights[0]]*2,
    #     '-', lw=2, label="Halbwertsbreite",
    # )

    # ---

    # get FWHM bounds from linear fits
    left = (max(I)/2 - b_l) / m_l
    right = (max(I)/2 - b_r) / m_r

    print(f"left: {left}")
    print(f"right: {right}")
    print(f"height: {max(I)/2}")
    print(f"Halbwertsbreite: {right - left}")

    plt2.plot(
        tools.pintify([left, right]),
        tools.pintify([max(I)/2]*2),
        '-',
        lw=2, label="Halbwertsbreite",
        show_xerr=False,
        show_yerr=False,
    )

# plt.ylim(bottom=0)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/1_koinzidenz.pdf")
plt.show()
