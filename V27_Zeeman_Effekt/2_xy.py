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


FOO = [
    ('rot', ureg('5 A')),
    ('blau_pi', ureg('5 A')),
    ('blau_sigma', ureg('3.24 A')),
]

for name, I in FOO:
    print(f'█ {name}')
    # █ Bestimmung der Wellenlängenaufspaltung
    ordnung, Δs, δs = np.genfromtxt(f'data/{name}.csv', delimiter=',', skip_header=1, unpack=True)
    δλ = δs * delta_λ_D / (2 * Δs)
    print(f"{δλ.mean()=}")

    # █ Bestimmung der Landé-Faktoren
    # I = Stromstärke des Elektromagneten
    B = ureg('302.1 mT')  # TODO: Nutze Daten aus 1_magnet.py

    # g_ij = m_j * g_j - m_i * g_i
    μ_B = ureg.e * ureg.hbar / (2 * ureg.m_e)
    g_ij = δλ.mean() * ureg.h * ureg.c / (λ**2 * μ_B * B)  # Landé-Faktor
    g_ij.ito('dimensionless')
    print(f"{g_ij=}")
