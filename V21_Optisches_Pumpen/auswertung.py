import generate_table
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
    """Magnetfeld im Zentrum eines Helmholtz-Spulenpaars aus Radius `r`, Windungszahl `N` und Strom `I`"""
    B = 8 * ureg.mu_0 * I * N / (np.sqrt(125) * r)
    return B.to('µT')


def calc_g(m):
    g = 4 * np.pi * ureg.m_e / (m * ureg.e)
    return g.to('dimensionless')


def calc_kernspin(g_F, J):
    # TODO: BenediktSan hat die 2 im Nenner
    return J * (2/g_F - 1)


def calc_quad_zeemann(g_F, B, E_HFS, m_f):
    E = g_F**2 * (ureg.e * ureg.h/(4 * np.pi * ureg.m_e))**2 * B**2 * (1 - 2 * m_f)/E_HFS
    return E.to('J')


# █ Daten einlesen
rf_freq, I1_hor, U1_sweep, I2_hor, U2_sweep = np.genfromtxt("data.txt", unpack=True)
rf_freq *= ureg.kHz
I1_hor *= ureg.mA
I2_hor *= ureg.mA
U1_sweep *= ureg.V
U2_sweep *= ureg.V

# Unsicherheiten einfügen
# I1_hor = unp.uarray(U1_hor, 0.02)
# I2_hor = unp.uarray(U2_hor, 0.02)
# U1_sweep = unp.uarray(U1_sweep, 0.01)
# U2_sweep = unp.uarray(U2_sweep, 0.01)


# █ Daten zu den Spulen [Versuchsanleitung]
# Sweepspule:
dat_sweep = {
    'r': ureg('16.39 cm'),  # Radius
    'N': 11,  # Windungen
    'R': 10 * ureg.ohm,  # Widerstand # [Quelle unbekannt]
}
# Horizontalspule:
dat_hor = {
    'r': ureg('15.79 cm'),
    'N': 154,
    'R':  10/3 * ureg.ohm,  # [Quelle unbekannt]
}
# Vertikalspule:
dat_vert = {
    'r': ureg('11.735 cm'),
    'N': 20,
}

# Gesamt-Magnetfelder berechnen
# (das Vertikalfeld ist zu 0 kompensiert)
B1_ges = (  # Sweepspule + Horizontalspule
    calc_B_helmholtz(dat_sweep['r'], dat_sweep['N'], I=(U1_sweep / dat_sweep['R'])) +
    calc_B_helmholtz(dat_hor['r'], dat_hor['N'], I1_hor)
)
B2_ges = (  # Sweepspule + Horizontalspule
    calc_B_helmholtz(dat_sweep['r'], dat_sweep['N'], I=(U2_sweep / dat_sweep['R'])) +
    calc_B_helmholtz(dat_hor['r'], dat_hor['N'], I2_hor)
)

# █ Tabelle generieren
# generate_table.generate_table_pint(
#     'build/tab/messwerte.tex',
#     (r'B1_\text{hor}', ureg.microtesla, B1_hor),
#     (r'B2_\text{hor}', ureg.microtesla, B2_hor),
# )

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
    plt2.plot(rf_freq, B2_ges, 'x', zorder=5, label=r"Messwerte zu $^{85}$Rb")
    plt2.plot(rf_freq_linspace, tools.nominal_values(
        params_85[0]*rf_freq_linspace + params_85[1]), label=r"Ausgleichsgerade zu $^{85}$Rb")

    plt2.plot(rf_freq, B1_ges, 'x', zorder=5, label=r"Messwerte zu $^{87}$Rb")
    plt2.plot(rf_freq_linspace, tools.nominal_values(
        params_87[0]*rf_freq_linspace + params_87[1]), label=r"Ausgleichsgerade zu $^{87}$Rb")
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/g_F.pdf")
# plt.show()

# █ Berechnung
g_F_87 = calc_g(params_87[0])
g_F_85 = calc_g(params_85[0])
print(tools.fmt_compare_to_ref(g_F_87, 1/2, name='g_F (87Rb)'))
print(tools.fmt_compare_to_ref(g_F_85, 1/3, name='g_F (85Rb)'))

console.rule("Kernspins [e) in der Versuchsanleitung]")
I_1 = calc_kernspin(g_F_87, J=1/2)
I_2 = calc_kernspin(g_F_85, J=1/2)
print(tools.fmt_compare_to_ref(I_1, 1.5, name='I (87Rb)'))
print(tools.fmt_compare_to_ref(I_2, 2.5, name='I (85Rb)'))


console.rule("Erdmagnetfeld: vertikal")
# Vertikale Magnetfeldkomponente in µT aus Spannung der vertikalen Spule
U_vert = unp.uarray(2.26, 0.01) * ureg.V  # TODO: eigene Werte!
B_vert = calc_B_helmholtz(dat_vert['r'], dat_vert['N'], I=(U_vert / dat_sweep['R'])) # TODO: warum dat_sweep?
print(tools.fmt_compare_to_ref(B_vert, ureg('40 µT')))  # TODO: Quelle unbekannt


console.rule("Erdmagnetfeld: horizontal")
# Horizontale Magnetfeldkomponente aus y-Achsenabschnitt des lin. Zusammenhangs zwischen B-Feld und RF-Frequenz
B_hor = (params_87[1] + params_85[1])/2  # Mittelwert aus den Achsenabschnitten der beiden Regressionsgeraden
print(tools.fmt_compare_to_ref(B_hor, ureg('20 µT')))  # TODO: Quelle unbekannt


console.rule("Isotopenverhältnis")
# Idee: Peak-Tiefe ~ Anteil des Isotops
# Das Verhältnis der beiden Isotope lässt sich anhand des Amplitudenverhältnis der beiden Dips (1. Dip 87Rb und 2. Dip 85Rb) ablesen.
print("[…]")


console.rule("Abschätzung des quadratischen Zeemann-Effekts")

# Hyperfeinstrukturaufspaltung des Grundzustandes [Quelle: Versuchsanleitung, Messprogramm h)]:
E_HFS_87 = 4.53e-24 * ureg.J
E_HFS_85 = 2.01e-24 * ureg.J

# quadratischer Zeeman-Effekt
E_QZ_87 = calc_quad_zeemann(g_F_87, max(B1_ges), E_HFS_87, m_f=1)
E_QZ_85 = calc_quad_zeemann(g_F_85, max(B2_ges), E_HFS_85, m_f=1)

#print("Aufspaltung durch den quadratischen Zeemann-Effekt bei höheren Magnetfeldstärken")
print(f"87Rb: {abs(E_QZ_87)} bei {max(B1_ges):.2f}")
print(f"85Rb: {abs(E_QZ_85)} bei {max(B2_ges):.2f}")
