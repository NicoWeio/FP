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
α, N = np.genfromtxt("data/1_detektorscan.txt", unpack=True)  # skip_header=1

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


def fit_fn(α, I_max, I_0, σ, α_0):
    # ↓ Hier entspricht σ nicht der gängigen Definition
    # return I_max * np.exp(-(σ * (α - α_0))**2) + I_0

    # ↓ Mampfzwerg
    # return (I_max / np.sqrt(2 * np.pi * σ**2)) * np.exp(-((α - α_0)**2) / (2 * σ**2)) + I_0

    # return I_max * np.exp(-((α/σ)**2)) + I_0 # kinda works
    return I_max * np.exp(-((α - α_0)**2) / (2 * σ**2)) + I_0  # works!!!

    # mit uncertainties-Dings
    # return (I_0 / np.sqrt(2 * np.pi * σ**2)) * np.exp((-((α - α_0)**2) / (2 * σ**2)).to('dimensionless')) + I_max


I_max, I_0, σ, α_0 = tools.pint_curve_fit(
    fit_fn,
    # α, I,
    α, tools.nominal_values(I),
    (1 / ureg.s, 1 / ureg.s, ureg.deg, ureg.deg),
    p0=(tools.nominal_value(max(I)), tools.nominal_value(min(I)), ureg('0.05°'), ureg('0°')),
    # return_p0=True,  # TODO
)
print(f'I_max = {I_max}')
print(f'I_0 = {I_0}')
print(f'σ = {σ}')
print(f'α_0 = {α_0}')

# TODO: Halbwertsbreite
# Die Standardabweichung σ {\displaystyle \sigma } \sigma beschreibt die Breite der Normalverteilung.
# Die Halbwertsbreite einer Normalverteilung ist ungefähr das 2,4-Fache $2{\sqrt {2\ln 2}}$ der Standardabweichung.

# █ Plot
α_linspace = tools.linspace(*tools.bounds(α), 1000)

plt.figure()
with tools.plot_context(plt, '°', '1/s', r"\alpha", "I") as plt2:
    plt2.plot(α, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

    # plt2.plot(α_linspace, fit_fn(α_linspace, *[I_max, I_0, σ, α_0]), label="Fit")
    plt2.plot(α_linspace, fit_fn(α_linspace, *map(tools.nominal_value, [I_max, I_0, σ, α_0])), label="Fit")

    # plt2.plot(
    #     tools.pintify([left, right]), [width_heights[0]]*2,
    #     '-', lw=2, label="Halbwertsbreite",
    # )

plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig("build/plt/1_detektorscan.pdf")
plt.show()
