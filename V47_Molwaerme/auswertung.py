import generate_table
import matplotlib.pyplot as plt
import tools
import numpy as np
import pint
from rich.console import Console
import scipy as sp
import tools
console = Console()
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

# █ Formeln


def calc_T(R):
    """Berechnet die Temperatur `T` aus dem Widerstand `R`"""
    # Quelle: Versuchsanleitung
    R_ = R.to('ohm').m
    T_ = 0.00134 * R_**2 + 2.296 * R_ - 243.02
    T = T_ * ureg.degC
    return T.to('K')


# █ Konstanten
# Masse der Probe
m = ureg('342 g')  # Quelle: Versuchsanleitung
# Kompressionsmodul Kupfer
# > Achtung: In der Versuchsanleitung steht für den Kompressionsmodul ein κ,
# > welches eigentlich für die Kompressibilität (=1/Kompressionsmodul) steht!
K = ureg('140 GPa')  # Quelle: https://periodictable.com/Elements/029/data.html → „Bulk Modulus“
# Molvolumen Kupfer
V0 = 7.0922e-6 * ureg('m^3/mol')  # Quelle: https://periodictable.com/Elements/029/data.html → „Molar Volume“
# Molare Masse Kupfer
M = 63.546e-3 * ureg('kg/mol')  # Quelle: https://periodictable.com/Elements/029/data.html → „Atomic Weight“
# Stoffmenge der Probe
n = m / M
print(f"n = {n.to('mol'):.2f}")
# longitudinale Phasengeschwindigkeit in Kupfer
v_long = ureg('4.7 km/s')  # Quelle: Versuchsanleitung
# transversale Phasengeschwindigkeit in Kupfer
v_trans = ureg('2.26 km/s')  # Quelle: Versuchsanleitung


# █ Messwerte einlesen
dt, U, I, R_probe, R_zylinder = np.genfromtxt('dat/messung.txt', unpack=True)
dt *= ureg.s
U *= ureg.V
I *= ureg.mA
R_probe *= ureg.kiloohm
R_zylinder *= ureg.kiloohm

# dt sei die Zeitdifferenz zwischen aktueller und nächster Zeile;
# somit wird der letzte Eintrag für dt verworfen, wenn Differenzen betrachtet werden.
# dt = dt[:-1]

# Temperaturen aus den Widerständen berechnen
T_probe = calc_T(R_probe)
T_zylinder = calc_T(R_zylinder)  # NOTE: War vorher auch fTp – scheinbar ein Fehler bei KarlSchiller

# █ a) Man messe die Molwärme C_p von Kupfer in Abhängigkeit von der Temperatur im Bereich von ca 80 bis 300 K.
console.rule('a)')

# elektrische Arbeit pro Zeitintervall
ΔW_el = U * I * dt
assert ΔW_el.check('[energy]')

# C_p aus den Messwerten
ΔT_probe = np.diff(T_probe)
C_p = ΔW_el[:-1] / (ΔT_probe * n)
assert C_p.check('J/(mol·K)')


# █ b) Man errechne daraus C_v […]
console.rule('b)')

# Einlesen der Werte von alpha, aus der Tabelle der Anleitung
tab_T_ɑ, tab_ɑ = np.genfromtxt('dat/Tab_2.csv', delimiter=',', skip_header=1, unpack=True)
# tab_T_ɑ *= ureg.K
# tab_ɑ *= 1e-6 / ureg.delta_degC


def calc_ɑ(T):
    """
    Berechne den linearen Ausdehnungskoeffizienten `ɑ` aus der Temperatur `T`
    durch lineare Interpolation zwischen den in Tabelle 2 der Versuchsanleitung gegebenen Werten.
    """
    return np.interp(T.to('K').m, tab_T_ɑ, tab_ɑ) * 1e-6 / ureg.delta_degC


# Berechne C_V mittels Korrekturformel
# T_avg = (iT_probe.to('K') + T_probe.to('K'))/2
T_avg = T_probe.to('K')[:-1]  # TODO…
C_V = C_p - 9 * calc_ɑ(T_avg)**2 * K * V0 * T_avg  # Quelle: Versuchsanleitung
assert C_V.check('J/(mol·K)')

# „Man berücksichtige hierfür nur Messwerte bis T_max […]“
T_max = ureg('170 K')

# █ c) Man versuche, die gemessenen (C_V,T)-Wertepaare durch Wahl einer geeigneten Debye-Temperatur θ_D in der universellen Debye-Kurve anzupassen.
console.rule('c)')
# Man berücksichtige hierfür nur Messwerte bis T_max = 170K.
# Welchen Wert für θ_D erhält man?

# Plotten von C_V
plt.figure()
with tools.plot_context(plt, 'K', 'J/(mol·K)', 'T', 'C_i') as plt2:  # C_{(\cdot)}
    # plt.errorbar(x=noms(Tmittel), xerr=stds(Tmittel), y=noms(C_V), yerr=stds(C_V), fmt='x', label='Stützstellen')
    plt2.plot(T_avg, C_p, 'x', zorder=5, label='gemessene Werte $C_p$')
    plt2.plot(T_avg, C_V, 'x', zorder=5, label='berechnete Werte $C_V$')

    plt.axvline(x=T_max, linestyle='--', color='grey', label=f'$T = {T_max.to("K"):Lx}$')
plt.grid()
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('build/plt/c_v.pdf')
# plt.show()


# → Tabelle
generate_table.generate_table_pint(
    'build/tab/mess.tex',
    (r'\symup{\Delta} t', ureg.s, dt, 0),
    ('U', ureg.V, U),
    ('I', ureg.mA, I),
    (r'R_\text{Probe}', ureg.ohm, R_probe),
    (r'T_\text{Probe}', ureg.kelvin, T_probe),
    (r'R_\text{Zylinder}', ureg.ohm, R_zylinder),
    (r'T_\text{Zylinder}', ureg.kelvin, T_zylinder),
)

# print(
#     T_probe.shape,
#     ΔT_probe.shape,
#     ΔW_el.shape,
#     C_p.shape,
#     C_V.shape,
# )

generate_table.generate_table_pint(
    'build/tab/warmekapazitaeten.tex',
    (r'T_\text{Probe}', ureg.kelvin, T_probe[:-1]),
    (r'\symup{\Delta} T_\text{Probe}', ureg.kelvin, ΔT_probe),
    (r'\symup{\Delta} E_\text{Probe}', ureg.joule, ΔW_el[:-1]),
    (r'C_p', (ureg.joule / ureg.mol / ureg.kelvin), C_p),
    (r'C_V', (ureg.joule / ureg.mol / ureg.kelvin), C_V),
)


# Bestimmung der Debye-Temperatur
T_max_mask = T_probe < T_max
C_V_avg = C_V[T_max_mask[:-1]].mean()
print(f"{C_V_avg.to('J/(K·mol)')=}")
T_avg = T_probe[T_max_mask].mean()
print(f"{T_avg=}")

# Versuchsanleitung, Tabelle 1
# Zahlenwerte der Debye-Funktion für R = 8,31439 J/(mol · °C)
# θ_D/T [–], C_V[J/(mol · °C)]
tab_θ_D_div_T, tab_C_V = np.genfromtxt('dat/Tab_1.csv', delimiter=',', skip_header=1, unpack=True)

# np.interp erwartet aufsteigende x-Werte
tab_θ_D_div_T = np.flip(tab_θ_D_div_T)
tab_C_V = np.flip(tab_C_V)

θ_D_div_T = np.interp(C_V_avg.to('J/(mol·K)').m, tab_C_V, tab_θ_D_div_T)
θ_D_div_T *= ureg.dimensionless
print(f"{θ_D_div_T=}")
θ_D = θ_D_div_T * T_avg.to('K')
print(f"{θ_D=}")

θ_D_lit = ureg('343 K')  # Quelle: Gross, Marx – Festkörperpysik (Abb. 6.9)
print(tools.fmt_compare_to_ref(θ_D, θ_D_lit, name='θ_D'))

# --

ω_D = np.cbrt(18 * np.pi**2 * ureg.N_A / V0 / (2/v_trans**3 + 1/v_long**3))
print(f"{ω_D.to('Hz'):.2e}")
