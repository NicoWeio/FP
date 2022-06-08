import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import tools
from uncertainties import ufloat
import generate_table
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
dt, U, I, R_probe, R_zylinder = np.genfromtxt('dat/messung.txt', unpack=True)
dt *= ureg.s
U *= ureg.V
I *= ureg.mA
R_probe *= ureg.ohm
R_zylinder *= ureg.ohm

# dt sei die Zeitdifferenz zwischen aktueller und nächster Zeile;
# somit wird der letzte Eintrag für dt verworfen, wenn Differenzen betrachtet werden.
# dt = dt[:-1]

# Temperaturen aus den Widerständen berechnen
T_probe = calc_T(R_probe)
T_zylinder = calc_T(R_zylinder)  # NOTE: War vorher auch fTp – scheinbar ein Fehler bei KarlSchiller

# █ a) Man messe die Molwärme C_p von Kupfer in Abhängigkeit von der Temperatur im Bereich von ca 80 bis 300 K.

# elektrische Arbeit pro Zeitintervall
W_el = U * I * dt
assert W_el.check('[energy]')

# C_p aus den Messwerten
C_p = W_el[:-1] / (np.diff(T_probe) * n)
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
plt.legend()
plt.tight_layout()
# plt.savefig('build/alpha.pdf')
# plt.show()

# Berechne C_V mittels Korrekturformel
# T_avg = (iT_probe.to('K') + T_probe.to('K'))/2
T_avg = T_probe.to('K')[:-1]  # TODO…
C_V = C_p - 9 * poly4(T_avg, *params)**2 * κ * V0 * T_avg  # Quelle: Versuchsanleitung
assert C_V.check('J/(mol·K)')

# █ c) Man versuche, die gemessenen (C_V,T)-Wertepaare durch Wahl einer geeigneten Debye-Temperatur θ_D in der universellen Debye-Kurve anzupassen.
# Man berücksichtige hierfür nur Messwerte bis T_max = 170K.
# Welchen Wert für θ_D erhält man?

# Plotten von C_V
plt.figure()
with tools.plot_context(plt, 'K', 'J/(mol·K)', 'T', 'C_V') as plt2:
    # plt.errorbar(x=noms(Tmittel), xerr=stds(Tmittel), y=noms(C_V), yerr=stds(C_V), fmt='x', label='Stützstellen')
    plt2.plot(T_avg, tools.nominal_values(C_V), 'o--', label='Messwerte')
    plt2.plot(T_avg, C_p, 'x', label='Messwerte $C_p$')
    # „Man berücksichtige hierfür nur Messwerte bis T_max“ […]“
    T_max = ureg('170 K')
    plt.axvline(x=T_max, linestyle='--', color='grey')
plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig('build/plt/cv.pdf')
plt.show()


# → Tabelle
generate_table.generate_table_pint(
    'build/tab/mess.tex',
    (r'\mathrm{\Delta} t', ureg.s, dt),
    ('U', ureg.V, U),
    ('I', ureg.mA, I),
    (r'R_\text{Probe}', ureg.ohm, R_probe),
    (r'R_\text{Zylinder}', ureg.ohm, R_zylinder),
)

