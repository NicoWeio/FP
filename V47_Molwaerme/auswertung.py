import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
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
# Loschmidtsche Zahl von CODATA
Nl = ufloat(2.6516467, 0.0000015)*1e25 * ureg('1/m^3')
# longitudinale Phasengeschwindigkeit in Kupfer
v_long = 4.7*1e3 * ureg('m/s')
# transversale Phasengeschwindigkeit in Kupfer
v_trans = 2.26*1e3 * ureg('m/s')
# Volumen der Probe
Vp = V0 * n * ureg('m^3')
# Avogadro-Konstante
# Na = ureg.avogadro_constant


# █ Messwerte einlesen
dt, U, I, iRp, iRz, fRp, fRz = np.genfromtxt('dat/messung_KarlSchiller.txt', unpack=True)
dt *= ureg.s
U *= ureg.V
I *= ureg.mA
iRp *= ureg.ohm
iRz *= ureg.ohm
fRp *= ureg.ohm
fRz *= ureg.ohm

# Temperaturen aus den Widerständen berechnen
iTp = calc_T(iRp)
iTz = calc_T(iRz)
fTp = calc_T(fRp)
fTp = calc_T(fRz)

# █ a) Man messe die Molwärme C_p von Kupfer in Abhängigkeit von der Temperatur im Bereich von ca 80 bis 300 K.

# elektrische Arbeit pro Zeitintervall
Wel = U * I * dt
assert Wel.check('[energy]')

# Cp aus den Messwerten
Cp = Wel / (abs(fTp - iTp) * n)
assert Cp.check('J/(mol·K)')


# █ b) Man errechne daraus C_v […]

# Einlesen der Werte von alpha, aus der Tabelle der Anleitung
Talpha, alpha = np.genfromtxt('dat/alpha_KarlSchiller.txt', unpack=True)
Talpha -= 273.15                # Umrechnen K in grd
alpha *= 1e-6                    # in 1/grd
Talpha *= ureg.degC
alpha /= ureg.delta_degC
Talpha.ito('K')  # https://pint.readthedocs.io/en/latest/nonmult.html

# Bestimmung einer allgemeinen Funktion von alpha
print('Regression für alpha')
params = tools.pint_polyfit(Talpha, alpha, 4)  # TODO: irgendwie mit Einheiten!
print(params)
a, b, c, d, e = params


def poly4(T, a, b, c, d, e):
    """Polynom 4. Ordnung"""
    return a * T**4 + b * T**3 + c * T**2 + d * T + e


# Plotten von alpha
plt.figure()
with tools.plot_context(plt, '°C', '1/°C', 'T', 'ɑ') as plt2:
    plt2.plot(Talpha, alpha, 'bx', label='Stützstellen')
    Tplot = np.linspace(Talpha[0], Talpha[-1], 500)
    plt2.plot(Tplot, tools.nominal_values(poly4(Tplot, *params)), 'k-', label='Regression')
plt.grid()
plt.tight_layout()
# plt.savefig('build/alpha.pdf')
plt.show()
plt.clf()

# Berechne Cv mittels Korrekturformel
Tmittel = (iTp.to('K') + fTp.to('K'))/2
Cv = Cp - 9 * poly4(Tmittel, *params)**2 * κ * V0 * Tmittel # Quelle: Versuchsanleitung
assert Cv.check('J/(mol·K)')

# Plotten von Cv
plt.figure()
with tools.plot_context(plt, '°C', 'J/(mol·°C)', 'T', 'V') as plt2:
    # plt.errorbar(x=noms(Tmittel), xerr=stds(Tmittel), y=noms(Cv), yerr=stds(Cv), color='b', fmt='x', label='Stützstellen')
    plt2.plot(Tmittel, tools.nominal_values(Cv), 'x', color='b', label='Stützstellen')
    # Tmax = ureg('170 K')
    # plt.axvline(x=Tmax.to('°C'), color='k')
plt.grid()
plt.tight_layout()
# plt.savefig('build/plt/cv.pdf')
plt.show()
plt.clf()
