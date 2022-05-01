import numpy as np
import pint
ureg = pint.UnitRegistry()

λ = ureg('632.816 nm')  # (mittlere) Wellenlänge des Übergangs (3s₂ - 2p₄)
f_0 = ureg.c / λ  # (mittlere) Frequenz
m = ureg('20.1797 u')  # molare Masse von Neon
T = ureg.Quantity(20, ureg.degC).to(ureg.kelvin)  # Umgebungstemperatur

# Halbwertsbreite, siehe z.B. https://de.wikipedia.org/wiki/Dopplerverbreiterung
δf = f_0 * np.sqrt(8 * ureg.k_B * T * np.log(2) / m) / ureg.c

print(f"f_0 = {f_0.to('THz'):.3f}")
print(f"δf = {δf.to('GHz'):.3f}")
