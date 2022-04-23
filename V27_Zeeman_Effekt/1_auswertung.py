import bildanalyse
import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np
import pint
from rich.console import Console
import rich
# from rich import print
import tools
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

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

# Abmessungen der Lummer-Gehrcke-Platte:
d = ureg('4 mm')  # Dicke
L = ureg('120 mm')

FOO = [
    {
        'color': 'rot',
        'λ': ureg('643.8 nm'),
        'n': 1.4567,
        'rotate_deg': -2.3,
        'g_lit': 1,
        'images': [
            {
                'I': ureg('0 A'),
                'path': 'Bilder/rot neu/IMG_0001.JPG',
                'min_height': 0.6,  # entspricht der alten Analyse
            },
            {
                'I': ureg('8 A'),
                'path': 'Bilder/rot neu/IMG_0002.JPG',
                # 'min_distance': 70,
                # 'min_height': 1600/2248,
                #
                'min_distance': 1,
                'min_height': 0.4,
                'prominence': 0.2,
            },
        ],
    },
    {
        'color': 'blau',
        'λ': ureg('480.0 nm'),
        'n': 1.4635,
        'rotate_deg': -2.3,
        'g_lit': 2,  # TODO
        'images': [
            {
                'I': ureg('0 A'),
                'polarisation': 0,  # in °
                'path': 'Bilder/blau/3_4.6A_sigma/IMG_0027.JPG',
                'min_distance': 40,
            },
            {
                'I': ureg('3.4 A'),  # TODO
                'polarisation': 0,  # in °
                'path': 'Bilder/blau/3_4.6A_sigma/IMG_0028.JPG',
                'min_distance': 1,
                'min_height': 0.6,
                'prominence': 0.0,
            },
        ],
    },
    {
        'color': 'blau',
        'λ': ureg('480.0 nm'),
        'n': 1.4635,
        'rotate_deg': -2.3,
        'g_lit': 2,  # TODO
        'images': [
            {
                'I': ureg('0 A'),
                'polarisation': 0,  # in °
                'path': 'Bilder/blau/4/IMG_0029.JPG',
                'min_distance': 40,
            },
            {
                'I': ureg('3.4 A'),  # TODO
                'polarisation': 0,  # in °
                'path': 'Bilder/blau/4/IMG_0031.JPG',
                'min_distance': 1,
                'min_height': 0.6,
                'prominence': 0.0,
            },
        ],
    },

]

MESSREIHEN = [
    FOO[0],
    FOO[1],
    # FOO[2],
]

for messreihe in MESSREIHEN:
    console.rule(f"Messreihe {messreihe['color']}")

    img1 = bildanalyse.preprocess_image(messreihe['images'][0]['path'], rotate_deg=messreihe['rotate_deg'])
    img2 = bildanalyse.preprocess_image(messreihe['images'][1]['path'], rotate_deg=messreihe['rotate_deg'])

    Δs = bildanalyse.get_Δs(
        img1,
        min_distance=messreihe['images'][0].get('min_distance', 100),
        min_height=messreihe['images'][0].get('min_height', 0.4),
        prominence=messreihe['images'][0].get('prominence', 0.2),
        show=False,
    )
    δs = bildanalyse.get_δs(
        img2,
        min_distance=messreihe['images'][1].get('min_distance', 100),
        min_height=messreihe['images'][1].get('min_height', 0.4),
        prominence=messreihe['images'][1].get('prominence', 0.2),
        show=False,
    )
    print(f"{Δs=}")
    print(f"{δs=}")

    Δλ_D = messreihe['λ']**2 / (2*d * np.sqrt(messreihe['n']**2 - 1))
    Δλ_D.ito('pm')
    print(f"Δλ_D ({messreihe['color']}) = {Δλ_D:.3f}")

    δλ = δs * Δλ_D / (2 * Δs)
    print(f"{δλ.mean()=}")

    # █ Bestimmung der Landé-Faktoren
    B = calc_B(messreihe['images'][1]['I'])

    # g_ij = m_j * g_j - m_i * g_i
    μ_B = ureg.e * ureg.hbar / (2 * ureg.m_e)
    g_ij = δλ.mean() * ureg.h * ureg.c / (messreihe['λ']**2 * μ_B * B)  # Landé-Faktor
    g_ij.ito('dimensionless')
    print(f"{g_ij=}")

    print(tools.fmt_compare_to_ref(g_ij, messreihe['g_lit']))
