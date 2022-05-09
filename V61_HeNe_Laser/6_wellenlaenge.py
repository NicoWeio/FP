import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
import rich
from rich.console import Console
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

DATA = [
    {
        'g': ureg('80 / mm'),
        'e': ureg('63.0 cm'),
    },
    {
        'g': ureg('100 / mm'),  # Spalte pro Längeneinheit
        'e': ureg('90.5 cm'),  # Abstand Gitter—Schirm
    },
    {
        'g': ureg('600 / mm'),
        'e': ureg('30.7 cm'),
    },
    {
        'g': ureg('1200 / mm'),
        'e': ureg('33.0 cm'),
    },
]
λ_lit = ureg('632.816 nm')  # theoretische Wellenlänge

global_λ = []
global_all_λ = []

for data in DATA:
    console.rule(f"{data['g']}")
    k, a_k = np.genfromtxt(f"dat/6_wellenlaenge_{data['g'].m:.0f}.csv", delimiter=',', skip_header=1, unpack=True)
    a_k *= ureg.cm  # Distanz vom Hauptmaximum
    print(f"({len(a_k)} Werte)")

    # 0-Werte rauswerfen
    a_k = a_k[k != 0]
    k = k[k != 0]

    # # TEST: Symmetrische Messung simulieren (setzt vorzeichenbehaftetes k voraus)
    # largest_common_k = min(max(k), abs(min(k)))
    # # Hinweis ausgeben
    # discarded_k = k[abs(k) > largest_common_k]
    # if len(discarded_k):
    #     print(f"{len(discarded_k)} k-Werte wurden entfernt: {discarded_k}")
    # # Darüber hinausgehende Werte abschneiden
    # a_k = a_k[abs(k) <= largest_common_k]
    # k = k[abs(k) <= largest_common_k]

    # |k|
    k = np.abs(k)

    # https://www.leifiphysik.de/optik/beugung-und-interferenz/grundwissen/vielfachspalt-und-gitter
    d = 1 / data['g']  # Spaltabstand
    e = ureg('90.5 cm')  # Abstand Gitter—Schirm

    # Näherung:
    # all_λ = d * a_k / (k * e)
    # exakte Formel:
    all_λ = d * a_k / (k * np.sqrt(e**2 + a_k**2))
    print(f"all_λ = {all_λ.to('nm'):.1f}")

    λ = tools.ufloat_from_list(abs(all_λ).to('nm'))
    print(tools.fmt_compare_to_ref(λ, λ_lit))
    global_λ.append(λ)
    global_all_λ.extend(all_λ)

console.rule("Mittelwert A")
λ_mean = tools.pintify(global_λ).to('nm').m.mean() * ureg('nm')  # Workaround: die einzelnen λ haben bereits eine Unsicherheit
print(tools.fmt_compare_to_ref(λ_mean, λ_lit))

console.rule("Mittelwert B")
λ_mean = tools.pintify(global_all_λ).to('nm').m.mean() * ureg('nm')
print(tools.fmt_compare_to_ref(λ_mean, λ_lit))
