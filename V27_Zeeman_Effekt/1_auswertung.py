import matplotlib.pyplot as plt
import numpy as np
import pint
import tools
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()

# █ Eichung des Elektromagneten

I, B = np.genfromtxt('data/1_magnet.csv', delimiter=',', skip_header=1, unpack=True)
I *= ureg('A')
B *= ureg('mT')

slope, intercept = tools.linregress(I, B)
print(f"{slope=}, {intercept=}")


# def fit_fn(I, a, b, c):
#     return a*I**2 + b*I + c
# a, b, c = tools.pint_curve_fit(fit_fn, I, B, (ureg('mT/A²'), ureg('mT/A'), ureg('mT')))


def calc_B(I):
    return slope * I + intercept
    # return fit_fn(I, a, b, c)


with tools.plot_context(plt, 'A', 'mT', 'I', 'B') as plt2:
    plt2.plot(I, B, 'x', zorder=5, label='Messwerte')
    plt2.plot(I, tools.nominal_values(calc_B(I)), label='Regressionsgerade')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('build/plt/1_magnet.pdf')
plt.plot()


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
    B = calc_B(I)

    # g_ij = m_j * g_j - m_i * g_i
    μ_B = ureg.e * ureg.hbar / (2 * ureg.m_e)
    g_ij = δλ.mean() * ureg.h * ureg.c / (λ**2 * μ_B * B)  # Landé-Faktor
    g_ij.ito('dimensionless')
    print(f"{g_ij=}")
