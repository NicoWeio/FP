import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

# █ Berechnung der Dispersionsgebiete

DATA = [
    (ureg('643.8 nm'), 1.4567, 'rot'),
    (ureg('480.0 nm'), 1.4635, 'blau'),
]

# Abmessungen der Lummer-Gehrcke-Platte:
d = ureg('4 mm')  # Dicke
L = ureg('120 mm')

for λ, n, color in DATA:
    # TODO: anderswo steht n**2 + 1
    delta_λ_D = λ**2 / (2*d * np.sqrt(n**2 - 1))
    delta_λ_D.ito('pm')
    print(f'Δλ_D ({color}) = {delta_λ_D:.3f}')


# █ Bestimmung der Wellenlängenaufspaltung

ordnung, Δs, δs = np.genfromtxt('blau_sigma.csv', delimiter=',', skip_header=1, unpack=True)
δλ = δs * delta_λ_D / (2 * Δs)

print(f"{δλ=}")

# █ Bestimmung der Landé-Faktoren

# TODO: für blau_sigma, blau_pi, rot ↓

I = ureg('5 A')  # Stromstärke des Elektromagneten
B = ureg('302.1 mT')  # TODO: Nutze Daten aus 1_magnet.py

# Landé-Faktor:
# g_ij = m_j * g_j - m_i * g_i

μ_B = ureg.e * ureg.hbar / (2 * ureg.m_e)

g_ij = δλ.mean() * ureg.h * ureg.c / (λ**2 * μ_B * B)
g_ij.ito('dimensionless')
print(f"{g_ij=}")
