import matplotlib.pyplot as plt
import numpy as np
import pint
import tools
import uncertainties.unumpy as unp
ureg = pint.UnitRegistry()


def b_helmholtz(r, N, I):
    B = 8 * ureg.mu_0 * I * N / (np.sqrt(125) * r)
    return B


def g(m):
    return 4 * np.pi * ureg.m_e / (m * ureg.e)


def calc_kernspin(g_f, J):
    return J * (2/g_f - 1)


def quad_zeemann(g_f, B, E_HFS, m_f):
    return g_f**2 * (ureg.e * ureg.h/(4 * np.pi * ureg.m_e))**2 * B**2 * (1 - 2 * m_f)/E_HFS


# Daten einlesen
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
R_sweep = 10 * ureg.ohm  # Widerstand

# Horizontalspule:
r_hor = ureg('15.79 cm')
N_hor = 154
R_hor = 10/3 * ureg.ohm

# Vertikalspule:
r_ver = ureg('11.735 cm')
N_ver = 20


I1_hor = U1_hor/R_hor
I2_hor = U2_hor/R_hor
I1_sweep = U1_sweep/R_sweep
I2_sweep = U2_sweep/R_sweep


# ---


B1_hor = b_helmholtz(r_hor, N_hor, I1_hor)
B1_sweep = b_helmholtz(r_sweep, N_sweep, I1_sweep)
B1_ges = B1_hor + B1_sweep

B2_hor = b_helmholtz(r_hor, N_hor, I2_hor)
B2_sweep = b_helmholtz(r_sweep, N_sweep, I2_sweep)
B2_ges = B2_hor + B2_sweep
print(B2_ges)
#print("Gesamte horizontale Magnetfeldstärke in Tesla")
print(f"Dip1: {max(B1_ges), np.argmax(B1_ges)}, Dip2: {max(B2_ges), np.argmax(B2_ges)}")


# ---


params_87 = tools.linregress(rf_freq, B1_ges)
params_85 = tools.linregress(rf_freq, B2_ges)

print(f"m, b von 85Rb: {params_85}")
print(f"m, b von 87Rb: {params_87}")

# Kernspin aus Steigung berechnen

rf_freq_linspace = tools.linspace(*tools.bounds(rf_freq))
with tools.plot_context(plt, 'kHz', 'µT', 'f', 'B') as plt2:
    plt2.plot(rf_freq_linspace, tools.nominal_values(params_85[0]*rf_freq_linspace + params_85[1]), label=r"$^{85}$Rb")
    plt2.plot(rf_freq, B2_ges, 'x')  # TODO: Errorbars sollten erscheinen, sobald U eine Unsicherheit erhält.

    plt2.plot(rf_freq_linspace, tools.nominal_values(params_87[0]*rf_freq_linspace + params_87[1]), label=r"$^{87}$Rb")
    plt2.plot(rf_freq, B1_ges, 'x')  # TODO: Errorbars sollten erscheinen, sobald U eine Unsicherheit erhält.
plt.legend()
# plt.savefig("build/plt/Rb_85_87.pdf")
# plt.show()


# g-Faktoren berechnen

g_f1 = g(params_87[0])
g_f2 = g(params_85[0])
print(f"87Rb: g_f={g_f1}")
print(f"85Rb: g_f={g_f2}")

# ---

I_1 = calc_kernspin(g_f1, 1/2)
I_2 = calc_kernspin(g_f2, 1/2)

print(f"87Rb: I={I_1}")
print(f"85Rb: I={I_2}")


# Erdmagnetfeld in vertikaler und horizontaler Richtung bestimmen

U_vert = unp.uarray(2.26, 0.01) * ureg.V
I_vert = U_vert/R_sweep
B_vert = b_helmholtz(r_ver, N_ver, I_vert)
print("Vertikale Magnetfeldkomponente in µT aus Spannung der vertikalen Spule")
print(B_vert)


# ---


B_hor = (params_87[1] + params_85[1])/2
print("Horizontale Magnetfeldkomponente in µT aus y-Achsenabschnitt des lin. Zusammenhangs zwischen B-Feld und RF-Frequenz")
print(B_hor)


# Vergleichswerte:
# horizontale Magntefledkomponente ~20 µT; vertikale Magnetfeldkomponente ~40 µT

# Das Verhältnis der beiden Isotope lässt sich anhand des Amplitudenverhältnis der beiden Dips (1. Dip 87Rb und 2. Dip 85Rb) ablesen.
# Es beträgt 87Rb = 0.5 * 85Rb

# Abschätzung des quadratischen Zeemann-Effekts:


# HFS = Hyperfeinstruktur
# Hyperfeinstrukturaufspaltung des Grundzustandes [Versuchsanleitung, Messprogramm h)]:
E_HFS_87 = 4.53e-24
E_HFS_85 = 2.01e-24

# QZ = quadratischer Zeeman-Effekt
# TODO: Argmax (9) nicht hardcoden ↓
E_QZ_87 = quad_zeemann(g_f1, B1_ges[9], E_HFS_87, 1)
E_QZ_85 = quad_zeemann(g_f2, B2_ges[9], E_HFS_85, 1)

#print("Aufspaltung durch den quadratischen Zeemann-Effekt bei höheren Magnetfeldstärken")
print(f"87Rb:{abs(E_QZ_87)} J bei {max(B1_ges)}, 85Rb: {abs(E_QZ_85)} J bei {max(B2_ges)}")