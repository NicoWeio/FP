# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
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

# █ Schichtdicke bestimmen (Peaks finden)

# Versuchsgrößen
k = 2*np.pi / λ  # Wellenvektor
qz = 2*k * np.sin(α)  # Wellenvektorübertrag -> y-Werte der Theoriekurve

# Parameter des Parratt-Algorithmus

# Brechungsindizes
d1 = 0.7e-6  # Polysterol -> Amplitude vergrößert + negativer Offset
d2 = 6.7e-6  # Silizium -> Amplitude vergkleinert + positiver Offset
#
n1 = 1  # Luft
n2 = 1 - d1  # Polysterol
n3 = 1 - d2  # Silizium

# Rauigkeit
s1 = 7.9e-10 * ureg.m  # Polysterol -> Amplitude verkleinert bei hinteren Oszillationen
s2 = 5.7e-10 * ureg.m  # Silizium -> Senkung des Kurvenendes und  Amplitudenverkleinerung der Oszillationen

# Schichtdicke
z = 855e-10 * ureg.m  # verkleinert Oszillationswellenlänge

# ----------------------------

q = 4 * np.pi / λ * np.sin(np.pi / 180 * α)  # TODO: Blind übernommen aus @Mampfzwerg

peaks, peak_props = sp.signal.find_peaks(tools.nominal_values(I).to('1/s').m, height=(1E2, 1E4), prominence=5)
Δα_mean = np.mean(np.diff(α[peaks].to('rad').m)) * ureg.rad
Δq_mean = np.mean(np.diff(q[peaks].to('1/m').m)) * ureg['1/m']
# d_estim_a = 2*np.pi / Δq_mean # anderes Resultat als bei @Mampfzwerg
d_estim_b = λ / (2 * Δα_mean)  # korrekte Daten bei @Mampfzwerg
print(f"Δα_mean = {Δα_mean}")
print(f"Δq_mean = {Δq_mean}")
# print(f"d_estim_a = {d_estim_a.to('m')}")
print(f"d_estim_b = {d_estim_b.to('m')}")

# █ Theoriekurve: Fresnelreflektivität Si


def parratt2(z, ai):
    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(ai)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(ai)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(ai)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2
    # Strecke vor Beginn der Oszillationen auf 1 setzen

    cutoff_i = 100  # TODO: anpassen (@Mampfzwerg: 296)
    par[:cutoff_i] = 1
    r13[:cutoff_i] = 1
    return par, r13


# z = 855e-10  # Parameter: Schichtdicke | verkleinert Oszillationswellenlänge
z = d_estim_b
par, r13 = parratt2(z, α.to('rad').m)

# █ Plot
# α_linspace = tools.linspace(*tools.bounds(α), 1000)

    plt.figure()
    # TODO: Doppelachse mit Intensität und Reflektivität?
    with tools.plot_context(plt, '1/m', '1/s', "q", "I") as plt2:
        plt2.plot(q, I_refl, fmt='-', zorder=5, label="Messwerte")  # oder 'x--'?
        plt2.plot(q, I_diff, fmt='-', zorder=5, label="Messwerte (diffuse)")
        plt2.plot(q, I, fmt='-', zorder=5, label="Messwerte (ohne diffuse)")

    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/{name}.pdf")
    plt.show()
