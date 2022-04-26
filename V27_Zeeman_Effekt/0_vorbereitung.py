import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np
import pint
import tools
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()


def calc_Δλ_D(d, n, λ):
    """
    Damit unterschiedliche Ordnungen sich nicht überlagern,
    dürfen zwei aufgespaltene Wellenlängen maximal eine Wellenlängendifferenz Δλ_D besitzen.
    (Näherung für große Austrittswinkel)

    d: Dicke der LG-Platte
    n: Brechungsindex
    λ: eingestrahlte Wellenlänge
    """
    Δλ_D = λ**2 / (2*d * np.sqrt(n**2 - 1))
    return Δλ_D.to('pm')


def calc_A(L, n, λ):
    """
    Auflösungsvermögen der Lummer-Gehrcke-Platte

    L: Länge der LG-Platte
    n: Brechungsindex
    λ: eingestrahlte Wellenlänge
    """
    A = L/λ * (n**2 - 1)
    return A.to('dimensionless')

# █ 6. Dispersionsgebiet und spektrale Auflösung der Messapparatur

DATA = [
    {
        'λ': ureg('643.8 nm'),
        'n': 1.4567,
        'minmax': (-1, 1),
    },
    {
        'λ': ureg('480.0 nm'),
        'n': 1.4635,
        'minmax': (-2, 2),
    },
]

d = ureg('4 mm')  # Dicke der LG-Platte
L = ureg('120 mm')  # Länge der LG-Platte

for data in DATA:
    Δλ_D = calc_Δλ_D(d, λ=data['λ'], n=data['n'])
    A = calc_A(L, λ=data['λ'], n=data['n'])
    print(f"→ {data['λ']:.1f}:")
    print(f"Δλ_D: {Δλ_D:.2f}")
    print(f"A: {A:.2f}")

# █ 7.


def calc_E0(λ):
    return ureg.h * ureg.c / λ


def calc_p(E0, delta_λ, min_, max_):
    return (E0 * delta_λ * (max_+min_) - ureg.h * ureg.c * (min_-max_)) / (delta_λ * ureg.mu_B * max_ * min_)


def calc_q(E0, min_, max_):
    return E0**2/(ureg.mu_B**2 * max_ * min_)


def calc_B_plus(p, q):
    return -p/2 + np.sqrt((p/2)**2-q)


def calc_B_minus(p, q):
    return -p/2-np.sqrt((p/2)**2-q)

for data in DATA:
    min_, max_ = data['minmax']
    print(f"→ {data['λ']:.1f}:")
    E0 = calc_E0(data['λ'])
    print(f'E0: {E0.to("eV"):.2f}')
    p = calc_p(E0, calc_Δλ_D(d, data['n'], data['λ']), min_, max_)
    print(f'p: {p:.2f}')
    q = calc_q(E0, min_, max_)
    print(f'q: {q:.2f}')
    B_plus = calc_B_plus(p, q)
    print(f'B_plus: {B_plus.to("T"):.3f}')
    B_minus = calc_B_minus(p, q)
    print(f'B_minus: {B_minus.to("T"):.3f}')
