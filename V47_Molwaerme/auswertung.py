import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
from uncertainties import ufloat
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

# █ Formeln
def calc_T(R):
    """Berechnet die Temperatur `T` aus dem Widerstand `R`"""
    # Quelle: Versuchsanleitung
    R_ = R.to('Ω').m
    T_ = 0.00134 * R_**2 + 2.296 * R_ - 243.02
    T = T_ * ureg.degC
    return T


# █ Konstanten
# Masse der Probe
m = 0.342 * ureg('kg')
# Kompressionsmodul Kupfer
kappa = 139e9 * ureg('N/m^2')
# Molvolumen Kupfer
V0 = 7.11e-6 * ureg('m^3/mol')
# Molare Masse Kupfer
M = 63.55*1e-3 * ureg('kg/mol')
# Stoffmenge der Probe
n = m / M * ureg('mol')
# Loschmidtsche Zahl von CODATA
Nl = ufloat(2.6516467, 0.0000015)*1e25 * ureg('1/m^3')
# longitudinale Phasengeschwindigkeit in Kupfer
vlong = 4.7*1e3 * ureg('m/s')
# transversale Phasengeschwindigkeit in Kupfer
vtrans = 2.26*1e3 * ureg('m/s')
# Volumen der Probe
Vp = V0 * n * ureg('m^3')
# Avogadro-Konstante
Na = ureg.avogadro_constant
