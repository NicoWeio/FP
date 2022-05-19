import matplotlib.pyplot as plt
import numpy as np
import pint
from rich.console import Console
import tools
import uncertainties.unumpy as unp
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()


def calc_B_helmholtz(r, N, I):
    B = 8 * ureg.mu_0 * I * N / (np.sqrt(125) * r)
    return B.to('µT')


def calc_g(m):
    g = 4 * np.pi * ureg.m_e / (m * ureg.e)
    return g.to('dimensionless')


def calc_kernspin(g_F, J):
    return J * (2/g_F - 1)


def calc_quad_zeemann(g_F, B, E_HFS, m_f):
    E = g_F**2 * (ureg.e * ureg.h/(4 * np.pi * ureg.m_e))**2 * B**2 * (1 - 2 * m_f)/E_HFS
    return E.to('J')


# █ Daten einlesen
rf_freq, U1_hor, U1_sweep, U2_hor, U2_sweep = np.genfromtxt("data.txt", unpack=True)
rf_freq *= ureg.Hz  # TODO: Einheit korrekt?
U1_hor *= ureg.V
U2_hor *= ureg.V
U1_sweep *= ureg.V
U2_sweep *= ureg.V

# Offset korrigieren
U1_hor -= ureg('13.8 V')
U2_hor -= ureg('13.8 V')

# Unsicherheiten einfügen
# U1_hor = unp.uarray(U1_hor, 0.02)
# U2_hor = unp.uarray(U2_hor, 0.02)
# U1_sweep = unp.uarray(U1_sweep, 0.01)
# U2_sweep = unp.uarray(U2_sweep, 0.01)


# █ Daten zu den Spulen [Versuchsanleitung]
# Sweepspule:
r_sweep = ureg('16.39 cm')  # Radius
N_sweep = 11  # Windungen
R_sweep = 10 * ureg.ohm  # Widerstand # [Quelle unbekannt]

# Horizontalspule:
r_hor = ureg('15.79 cm')
N_hor = 154
R_hor = 10/3 * ureg.ohm  # [Quelle unbekannt]

# Vertikalspule:
r_ver = ureg('11.735 cm')
N_ver = 20

I1_hor = U1_hor/R_hor
I2_hor = U2_hor/R_hor
I1_sweep = U1_sweep/R_sweep
I2_sweep = U2_sweep/R_sweep


B1_hor = calc_B_helmholtz(r_hor, N_hor, I1_hor)
B1_sweep = calc_B_helmholtz(r_sweep, N_sweep, I1_sweep)
B1_ges = B1_hor + B1_sweep

B2_hor = calc_B_helmholtz(r_hor, N_hor, I2_hor)
B2_sweep = calc_B_helmholtz(r_sweep, N_sweep, I2_sweep)
B2_ges = B2_hor + B2_sweep


dip1 = np.argmax(B1_ges)
dip2 = np.argmax(B2_ges)
print(f"Dip1: {max(B1_ges), dip1}, Dip2: {max(B2_ges), dip2}")

# TODO: Plot, weil ich nicht verstehe, was hier abgeht.
plt.figure()
plt.stackplot(rf_freq, B1_hor, B1_sweep, labels=['Horizontalspule', 'Sweepspule'])
plt.axvline(rf_freq[dip1], color='k', linestyle='--', label='Dip-Dings')
# plt.stackplot(rf_freq, B2_hor, B2_sweep, labels=['Horizontalspule', 'Sweepspule'])
# plt.axvline(rf_freq[dip2], color='k', linestyle='--', label='Dip-Dings')
plt.legend()
# plt.show()


console.rule("g-Faktoren [d) in der Versuchsanleitung]")
# █ lineare Regression
params_87 = tools.linregress(rf_freq, B1_ges)
params_85 = tools.linregress(rf_freq, B2_ges)
print(f"87Rb: (m, b)={params_87}")
print(f"85Rb: (m, b)={params_85}")

# █ Plot
rf_freq_linspace = tools.linspace(*tools.bounds(rf_freq))
plt.figure()
with tools.plot_context(plt, 'kHz', 'µT', 'f', 'B') as plt2:
    plt2.plot(rf_freq, B2_ges, 'x', zorder=5, label=r"Messwerte zu $^{85}$Rb")  # TODO: Errorbars sollten erscheinen, sobald U eine Unsicherheit erhält.
    plt2.plot(rf_freq_linspace, tools.nominal_values(params_85[0]*rf_freq_linspace + params_85[1]), label=r"Ausgleichsgerade zu $^{85}$Rb")

    plt2.plot(rf_freq, B1_ges, 'x', zorder=5, label=r"Messwerte zu $^{87}$Rb")  # TODO: Errorbars sollten erscheinen, sobald U eine Unsicherheit erhält.
    plt2.plot(rf_freq_linspace, tools.nominal_values(params_87[0]*rf_freq_linspace + params_87[1]), label=r"Ausgleichsgerade zu $^{87}$Rb")
plt.legend()
plt.savefig("build/plt/g_F.pdf")
# plt.show()

# █ Berechnung
g_F_87 = calc_g(params_87[0])
g_F_85 = calc_g(params_85[0])
print(f"87Rb: g_F={g_F_87}")
print(f"85Rb: g_F={g_F_85}")

console.rule("Kernspins [e) in der Versuchsanleitung]")
I_1 = calc_kernspin(g_F_87, J=1/2)
I_2 = calc_kernspin(g_F_85, J=1/2)
print(f"87Rb: I={I_1}")
print(f"85Rb: I={I_2}")


console.rule("Erdmagnetfeld: vertikal")
# Vertikale Magnetfeldkomponente in µT aus Spannung der vertikalen Spule
U_vert = unp.uarray(2.26, 0.01) * ureg.V  # TODO: eigene Werte!
I_vert = U_vert/R_sweep
B_vert = calc_B_helmholtz(r_ver, N_ver, I_vert)
print(f"B_vert = {B_vert}")
print(tools.fmt_compare_to_ref(B_vert, ureg('40 µT')))  # TODO: Quelle unbekannt


console.rule("Erdmagnetfeld: horizontal")
# Horizontale Magnetfeldkomponente aus y-Achsenabschnitt des lin. Zusammenhangs zwischen B-Feld und RF-Frequenz
B_hor = (params_87[1] + params_85[1])/2
print(f"B_hor = {B_hor}")
print(tools.fmt_compare_to_ref(B_hor, ureg('20 µT')))  # TODO: Quelle unbekannt


console.rule("Isotopenverhältnis")
# Das Verhältnis der beiden Isotope lässt sich anhand des Amplitudenverhältnis der beiden Dips (1. Dip 87Rb und 2. Dip 85Rb) ablesen.
# Es beträgt 87Rb = 0.5 * 85Rb
# print(f"Isotopenverhältnis…", max(B1_ges) / max(B2_ges)) # random guess
print("[…]")


console.rule("Abschätzung des quadratischen Zeemann-Effekts")

# HFS = Hyperfeinstruktur
# Hyperfeinstrukturaufspaltung des Grundzustandes [Versuchsanleitung, Messprogramm h)]:
E_HFS_87 = 4.53e-24 * ureg.J
E_HFS_85 = 2.01e-24 * ureg.J

# QZ = quadratischer Zeeman-Effekt
E_QZ_87 = calc_quad_zeemann(g_F_87, max(B1_ges), E_HFS_87, 1)
E_QZ_85 = calc_quad_zeemann(g_F_85, max(B2_ges), E_HFS_85, 1)

#print("Aufspaltung durch den quadratischen Zeemann-Effekt bei höheren Magnetfeldstärken")
print(f"87Rb: {abs(E_QZ_87)} bei {max(B1_ges):.2f}")
print(f"85Rb: {abs(E_QZ_85)} bei {max(B2_ges):.2f}")
