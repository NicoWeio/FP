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
    },
    {
        'λ': ureg('480.0 nm'),
        'n': 1.4635,
    },
]

d = ureg('4 mm')  # Dicke der LG-Platte
L = ureg('120 mm')  # Länge der LG-Platte

for data in DATA:
    Δλ_D = calc_Δλ_D(d, **data)
    A = calc_A(L, **data)
    print(f"→ {data['λ']:.1f}:")
    print(f"Δλ_D: {Δλ_D:.2f}")
    print(f"A: {A:.2f}")

# █ 7.


def calc_E0(λ):
    return ureg.h * ureg.c / λ


def p(E0, delta_λ, min_, max_):
    return (E0*delta_λ*(max_+min_)-h*c*(min_-max_)) / (delta_λ*mu_B*max_*min_)


def q(E0, min_, max_):
    return E0**2/(mu_B**2*max_*min_)


def B_plus(p, q):
    return -p/2 + np.sqrt((p/2)**2-q)


def B_minus(p, q):
    return -p/2-np.sqrt((p/2)**2-q)


max1 = 1
min1 = -1
max2 = 2
min2 = -2


for λ in [ureg('643.8 nm'), ureg('480.0 nm')]:
    E0 = calc_E0(λ)
    p1 = p(E0_1, delta_λ(λ1, A1), min1, max1)
    q1 = q(E0_1, min1, max1)
    B1p = B_p(p1, q1)
    B1m = B_m(p1, q1)
