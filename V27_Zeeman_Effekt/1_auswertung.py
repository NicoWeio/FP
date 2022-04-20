import bildanalyse
import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np
import pint
import tools
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

# █ Eichung des Elektromagneten

I, B = np.genfromtxt('data/1_magnet.csv', delimiter=',', skip_header=1, unpack=True)
I *= ureg('A')
B *= ureg('mT')

fit_params = tools.pint_polyfit(I, B, 3)
print(f"{fit_params=}")


def calc_B(I):
    a, b, c, d = fit_params
    return a*I**3 + b*I**2 + c*I + d


# with tools.plot_context(plt, 'A', 'mT', 'I', 'B') as plt2:
#     plt2.plot(I, B, 'x', zorder=5, label='Messwerte')
#     plt2.plot(I, tools.nominal_values(calc_B(I)), label='Regressionsgerade')
# plt.grid()
# plt.legend()
# plt.tight_layout()
# plt.savefig('build/plt/1_magnet.pdf')
# plt.plot()


# █ Berechnung der Dispersionsgebiete

DATA = [
    # (λ, n, color)
    (ureg('643.8 nm'), 1.4567, 'rot'),
    (ureg('480.0 nm'), 1.4635, 'blau'),
]

# Abmessungen der Lummer-Gehrcke-Platte:
d = ureg('4 mm')  # Dicke
L = ureg('120 mm')

for λ, n, color in DATA:
    Δλ_D = λ**2 / (2*d * np.sqrt(n**2 - 1))
    Δλ_D.ito('pm')
    print(f'Δλ_D ({color}) = {Δλ_D:.3f}')


FOO = [
    ('rot', ureg('8 A')),
    # ('blau_pi', ureg('5 A')),
    # ('blau_sigma', ureg('3.24 A')),
]

for name, I in FOO:
    print(f'█ {name}')
    # █ Bestimmung der Wellenlängenaufspaltung
    # ordnung, Δs, δs = np.genfromtxt(f'data/{name}.csv', delimiter=',', skip_header=1, unpack=True)

    img1 = imread('img/rot_0A.jpg')
    img2 = imread('img/rot_8A.jpg')

    Δs = bildanalyse.get_Δs(img1)
    δs = bildanalyse.get_δs(img2)

    δλ = δs * Δλ_D / (2 * Δs)
    print(f"{δλ.mean()=}")

    # █ Bestimmung der Landé-Faktoren
    B = calc_B(I)

    # g_ij = m_j * g_j - m_i * g_i
    μ_B = ureg.e * ureg.hbar / (2 * ureg.m_e)
    g_ij = δλ.mean() * ureg.h * ureg.c / (λ**2 * μ_B * B)  # Landé-Faktor
    g_ij.ito('dimensionless')
    print(f"{g_ij=}")
