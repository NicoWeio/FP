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


widths, width_heights, left_ips, right_ips = sp.signal.peak_widths(tools.nominal_values(I), [np.argmax(I)])

print(f"Peak: {max(I)} bei {Δt[np.argmax(I)]}")
print(f"Halbwertsbreite: {(width_heights[0] * I.units):.2f}")  # TODO…


# █ Plot
Δt_linspace = tools.linspace(*tools.bounds(Δt))

plt.figure()
with tools.plot_context(plt, 'ns', '1/s', r"\mathrm{\delta}t", "I") as plt2:  # TODO
    plt2.plot(Δt, I, fmt='x--', zorder=5, label="Messwerte")

    left_interp = left_ips[0]
    left_index = int(left_interp)
    left_frac = left_interp - left_index
    left = Δt[left_index] + left_frac * (Δt[left_index+1] - Δt[left_index])

    right_interp = right_ips[0]
    right_index = int(right_interp)
    right_frac = right_interp - right_index
    right = Δt[right_index] + right_frac * (Δt[right_index+1] - Δt[right_index])

    plt2.plot(
        tools.pintify([left, right]), [width_heights[0]]*2,
        '-', lw=2, label="Halbwertsbreite",
    )

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/1_koinzidenz.pdf")
plt.show()
