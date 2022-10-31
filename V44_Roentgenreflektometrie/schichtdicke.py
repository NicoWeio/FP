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

# ----------------------------

# Versuchsgrößen
l = 1.54e-10  # Wellenlänge
ai = np.pi/180 * np.arange(6e-2, 1.505, 5e-4)  # Winkel -> x−Werte der Theoriekurve
k = 2*np.pi / l  # Wellenvektor
qz = 2*k * np.sin(ai)  # Wellenvektorübertrag -> y-Werte der Theoriekurve

# Parameter des Parratt-Algorithmus

# Brechungsindizes
d1 = 0.7e-6  # Polysterol -> Amplitude vergrößert + negativer Offset
d2 = 6.7e-6  # Silizium -> Amplitude vergkleinert + positiver Offset
#
n1 = 1  # Luft
n2 = 1 - d1  # Polysterol
n3 = 1 - d2  # Silizium

# Rauigkeit
s1 = 7.9e-10  # Polysterol -> Amplitude verkleinert bei hinteren Oszillationen
s2 = 5.7e-10  # Silizium -> Senkung des Kurvenendes und  Amplitudenverkleinerung der Oszillationen

# Schichtdicke
z = 855e-10  # verkleinert Oszillationswellenlänge


def parratt(z):
    # TODO: Warum keine Rekursion? Weil wir nur zwei Schichten haben?

    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(ai)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(ai)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(ai)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * kz1 * kz2 * s1**2)
    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * kz2 * kz3 * s2**2)
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2
    # Strecke vor Beginn der Oszillationen auf 1 setzen
    for i in np.arange(np.size(par)):
        if (i <= 296):  # 296 manuell angepasst
            par[i] = 1
        else:
            pass
    return par

# ----------------------------


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
plt.savefig("build/plt/4_reflektivitätsscan.pdf")
plt.show()
