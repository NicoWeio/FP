import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

k, a_k = np.genfromtxt('dat/6_wellenlaenge_100.csv', delimiter=',', skip_header=1, unpack=True)
a_k *= ureg.cm  # Distanz vom Hauptmaximum

# https://www.leifiphysik.de/optik/beugung-und-interferenz/grundwissen/vielfachspalt-und-gitter

d = 1 / ureg('100 / mm')  # Spaltabstand
e = ureg('90.5 cm')  # Abstand Gitter—Schirm
λ_lit = ureg('632.816 nm')  # theoretische Wellenlänge


# TODO: exakte Formel

# Näherung:

# 0-Werte rauswerfen
a_k = a_k[k != 0]
k = k[k != 0]

all_λ = d * a_k / (k * e)
λ = abs(all_λ).mean().to('nm')
print(f"λ: {λ}")
print(tools.fmt_compare_to_ref(λ, λ_lit))
