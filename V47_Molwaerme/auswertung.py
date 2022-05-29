import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import tools
from uncertainties import ufloat
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

# █ Formeln


def calc_T(R):
    """Berechnet die Temperatur `T` aus dem Widerstand `R`"""
    # Quelle: Versuchsanleitung
    R_ = R.to('ohm').m
    T_ = 0.00134 * R_**2 + 2.296 * R_ - 243.02
    T = T_ * ureg.degC
    return T


# █ Konstanten
# Masse der Probe
m = 0.342 * ureg('kg')
# Kompressionsmodul Kupfer
κ = 139e9 * ureg('N/m^2')
# Molvolumen Kupfer
V0 = 7.11e-6 * ureg('m^3/mol')
# Molare Masse Kupfer
M = 63.55*1e-3 * ureg('kg/mol')
# Stoffmenge der Probe
n = m / M
# Loschmidtsche Zahl (CODATA)
Nl = sp.constants.value('Loschmidt constant (273.15 K, 100 kPa)') * ureg('1/m^3')
# longitudinale Phasengeschwindigkeit in Kupfer
v_long = ureg('4.7 km/s')
# transversale Phasengeschwindigkeit in Kupfer
v_trans = ureg('2.26 km/s')
# Volumen der Probe
V_probe = V0 * n * ureg('m^3')
# Avogadro-Konstante
# Na = ureg.avogadro_constant


# █ Messwerte einlesen
dt, U, I, iR_probe, iR_zylinder, fR_probe, fR_zylinder = np.genfromtxt('dat/messung_KarlSchiller.txt', unpack=True)
dt *= ureg.s
U *= ureg.V
I *= ureg.mA
# i = init, f = final
iR_probe *= ureg.ohm
iR_zylinder *= ureg.ohm
fR_probe *= ureg.ohm
fR_zylinder *= ureg.ohm

# Temperaturen aus den Widerständen berechnen
iT_probe = calc_T(iR_probe)
iT_zylinder = calc_T(iR_zylinder)
fT_probe = calc_T(fR_probe)
fT_zylinder = calc_T(fR_zylinder)  # NOTE: War vorher auch fTp – scheinbar ein Fehler bei KarlSchiller

# █ a) Man messe die Molwärme C_p von Kupfer in Abhängigkeit von der Temperatur im Bereich von ca 80 bis 300 K.

# elektrische Arbeit pro Zeitintervall
W_el = U * I * dt
assert W_el.check('[energy]')

# C_p aus den Messwerten
C_p = W_el / (abs(fT_probe - iT_probe) * n)
assert C_p.check('J/(mol·K)')


# █ b) Man errechne daraus C_v […]

# Einlesen der Werte von alpha, aus der Tabelle der Anleitung
T_ɑ, ɑ = np.genfromtxt('dat/Tab_2.csv', delimiter=',', skip_header=1, unpack=True)
T_ɑ *= ureg.K
ɑ *= 1e-6 / ureg.delta_degC

# Bestimmung einer allgemeinen Funktion von alpha
params = tools.pint_polyfit(T_ɑ, ɑ, deg=4)
print("Fit-Parameter für T_ɑ(ɑ):", params)


def poly4(T, a, b, c, d, e):
    """Polynom 4. Ordnung"""
    return a * T**4 + b * T**3 + c * T**2 + d * T + e


# Plotten von alpha
T_linspace = tools.linspace(*tools.bounds(T_ɑ), 500)
plt.figure()
with tools.plot_context(plt, '°C', '1/°C', 'T', 'ɑ') as plt2:
    plt2.plot(T_ɑ, ɑ, 'x', zorder=5, label='Messwerte')
    plt2.plot(T_linspace, tools.nominal_values(poly4(T_linspace, *params)), label='Fit')
plt.grid()
plt.tight_layout()
# plt.savefig('build/alpha.pdf')
# plt.show()

# Berechne Cv mittels Korrekturformel
T_avg = (iT_probe.to('K') + fT_probe.to('K'))/2
Cv = C_p - 9 * poly4(T_avg, *params)**2 * κ * V0 * T_avg  # Quelle: Versuchsanleitung
assert Cv.check('J/(mol·K)')

# █ c) Man versuche, die gemessenen (Cv,T)-Wertepaare durch Wahl einer geeigneten Debye-Temperatur θ_D in der universellen Debye-Kurve anzupassen.
# Man berücksichtige hierfür nur Messwerte bis T_max = 170K.
# Welchen Wert für θ_D erhält man?

# Plotten von Cv
plt.figure()
with tools.plot_context(plt, '°C', 'J/(mol·°C)', 'T', 'C_V') as plt2:
    # plt.errorbar(x=noms(Tmittel), xerr=stds(Tmittel), y=noms(Cv), yerr=stds(Cv), color='b', fmt='x', label='Stützstellen')
    plt2.plot(T_avg, tools.nominal_values(Cv), 'x', color='b', label='Stützstellen')
    # „Man berücksichtige hierfür nur Messwerte bis T_max“ […]“
    T_max = ureg('170 K')
    plt.axvline(x=T_max.to('°C'), linestyle='--', color='grey')
plt.grid()
plt.tight_layout()
# plt.savefig('build/plt/cv.pdf')
plt.show()
