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


def calc_g(ɑ):
    # g = 4 * np.pi * ureg.m_e / (ɑ * ureg.e) # ✓
    g = ureg.h / (ɑ * ureg.mu_B)
    return g.to('dimensionless')


def calc_kernspin(g_F):
    # S = 1/2
    # L = 0
    # J = 1/2
    g_J = 2.0023  # https://de.wikipedia.org/wiki/Land%C3%A9-Faktor#Elektron
    I = (g_J/g_F - 1) / 2
    return I.to('dimensionless')


def calc_quad_zeemann(g_F, B, E_HFS, m_F):
    E_lin = g_F * ureg.mu_B * B
    E_quad = g_F**2 * ureg.mu_B**2 * B**2 * (1 - 2*m_F) / E_HFS
    # E = E_lin + E_quad
    return (E_lin.to('J'), E_quad.to('J'))


# █ Daten einlesen
rf_freq, I_hor_87, U_sweep_87, I_hor_85, U_sweep_85 = np.genfromtxt("data.txt", unpack=True)
rf_freq *= ureg.kHz
I_hor_87 *= ureg.mA
I_hor_85 *= ureg.mA
U_sweep_87 *= ureg.V
U_sweep_85 *= ureg.V


# █ Daten zu den Spulen [Versuchsanleitung]
# Sweepspule:
dat_sweep = {
    'r': ureg('16.39 cm'),  # Radius
    'N': 11,  # Windungen
    'R': 10 * ureg.ohm,  # Widerstand („1 Umdrehung = 0.1 A“)
}
# Horizontalspule:
dat_hor = {
    'r': ureg('15.79 cm'),
    'N': 154,
    # Nicht verwendet, da das Potenziometer defekt war und auf eine externe Stromversorgung mit Amperemeter zurückgegriffen wurde. ↓
    # 'R': 10/3 * ureg.ohm,  # („1 Umdrehung = 0.3 A“)
}
# Vertikalspule:
dat_vert = {
    'r': ureg('11.735 cm'),
    'N': 20,
    'R': 10 * ureg.ohm,  # („1 Umdrehung = 0.1 A“)
}

console.rule("Erdmagnetfeld: vertikal")
# Vertikale Magnetfeldkomponente aus Spannung der vertikalen Spule
U_vert = unp.uarray(2.3, 0.01) * ureg.V  # TODO: eigene Werte!
B_vert = calc_B_helmholtz(dat_vert['r'], dat_vert['N'], I=(U_vert / dat_vert['R']))
print(tools.fmt_compare_to_ref(B_vert, ureg('44 µT')))  # Quelle: https://de.wikipedia.org/wiki/Erdmagnetfeld


# Gesamt-Magnetfelder berechnen
# (das Vertikalfeld ist hierbei zu 0 kompensiert)
I_sweep_87 = (U_sweep_87 / dat_sweep['R'])
B_ges_87 = (  # Sweepspule + Horizontalspule
    calc_B_helmholtz(dat_sweep['r'], dat_sweep['N'], I_sweep_87) +
    calc_B_helmholtz(dat_hor['r'], dat_hor['N'], I_hor_87)
)

I_sweep_85 = (U_sweep_85 / dat_sweep['R'])
B_ges_85 = (  # Sweepspule + Horizontalspule
    calc_B_helmholtz(dat_sweep['r'], dat_sweep['N'], I_sweep_85) +
    calc_B_helmholtz(dat_hor['r'], dat_hor['N'], I_hor_85)
)

# █ Tabelle generieren
# generate_table.generate_table_pint(
#     'build/tab/messwerte.tex',
#     (r'f', ureg.kHz, rf_freq, 0),
#     (r'I_\text{hor}', ureg.mA, I_hor_87),
#     (r'I_\text{sweep}', ureg.mA, I_sweep_87),
#     (r'B_\text{ges}', ureg.microtesla, B_ges_87),
#     (r'I_\text{hor}', ureg.mA, I_hor_85),
#     (r'I_\text{sweep}', ureg.mA, I_sweep_85),
#     (r'B_\text{ges}', ureg.microtesla, B_ges_85),
# )

console.rule("g-Faktoren [d) in der Versuchsanleitung]")
# █ lineare Regression
params_87 = tools.linregress(rf_freq, B_ges_87)
params_85 = tools.linregress(rf_freq, B_ges_85)
print(f"87Rb: (m, b)={params_87}")
print(f"85Rb: (m, b)={params_85}")

# █ Plot
rf_freq_linspace = tools.linspace(*tools.bounds(rf_freq))
plt.figure()
with tools.plot_context(plt, 'kHz', 'µT', 'f', 'B') as plt2:
    plt2.plot(rf_freq, B_ges_85, 'x', zorder=5, label=r"Messwerte zu $^{85}$Rb")
    plt2.plot(rf_freq_linspace, tools.nominal_values(
        params_85[0]*rf_freq_linspace + params_85[1]), label=r"Ausgleichsgerade zu $^{85}$Rb")

    plt2.plot(rf_freq, B_ges_87, 'x', zorder=5, label=r"Messwerte zu $^{87}$Rb")
    plt2.plot(rf_freq_linspace, tools.nominal_values(
        params_87[0]*rf_freq_linspace + params_87[1]), label=r"Ausgleichsgerade zu $^{87}$Rb")
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/g_F.pdf")
# plt.show()

# █ Berechnung
g_F_87 = calc_g(params_87[0])
g_F_85 = calc_g(params_85[0])
print(tools.fmt_compare_to_ref(g_F_87, 1/2, name='g_F (87Rb)', precision=3))
print(tools.fmt_compare_to_ref(g_F_85, 1/3, name='g_F (85Rb)', precision=3))

console.rule("Kernspins [e) in der Versuchsanleitung]")
I_1 = calc_kernspin(g_F_87)
I_2 = calc_kernspin(g_F_85)
print(tools.fmt_compare_to_ref(I_1, 1.5, name='I (87Rb)'))
print(tools.fmt_compare_to_ref(I_2, 2.5, name='I (85Rb)'))


console.rule("Erdmagnetfeld: horizontal")
# Horizontale Magnetfeldkomponente aus y-Achsenabschnitt des lin. Zusammenhangs zwischen B-Feld und RF-Frequenz
B_hor = (params_87[1] + params_85[1])/2  # Mittelwert aus den Achsenabschnitten der beiden Regressionsgeraden
print(tools.fmt_compare_to_ref(B_hor, ureg('20 µT')))  # Quelle: https://de.wikipedia.org/wiki/Erdmagnetfeld


console.rule("Isotopenverhältnis")
# Idee: Peak-Tiefe ~ Anteil des Isotops
# Das Verhältnis der beiden Isotope lässt sich anhand des Amplitudenverhältnis der beiden Dips (1. Dip 87Rb und 2. Dip 85Rb) ablesen.
print("[…]")


console.rule("Abschätzung des quadratischen Zeemann-Effekts")

# Hyperfeinstrukturaufspaltung des Grundzustandes [Quelle: Versuchsanleitung, Messprogramm h)]:
E_HFS_87 = 4.53e-24 * ureg.J
E_HFS_85 = 2.01e-24 * ureg.J

# quadratischer Zeeman-Effekt
E_QZ_87_parts = calc_quad_zeemann(g_F_87, max(B_ges_87), E_HFS_87, m_F=2)
E_QZ_85_parts = calc_quad_zeemann(g_F_85, max(B_ges_85), E_HFS_85, m_F=3)

#print("Aufspaltung durch den quadratischen Zeemann-Effekt bei höheren Magnetfeldstärken")
print(
    "87Rb:",
    f"linearer Term {E_QZ_87_parts[0]:.2e}",
    f"+ quadratischer Term {E_QZ_87_parts[1]:.2e}",
    f"= {sum(E_QZ_87_parts):.2e}",
    f"| B = {max(B_ges_87):.2f}",
)
print(
    "85Rb:",
    f"linearer Term {E_QZ_85_parts[0]:.2e}",
    f"+ quadratischer Term {E_QZ_85_parts[1]:.2e}",
    f"= {sum(E_QZ_85_parts):.2e}",
    f"| B = {max(B_ges_85):.2f}",
)
