import bildanalyse
import itertools
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

# TODO: Auslagern in tools / env-Variable
# steuert, ob Plots und Tabellen erstellt werden sollen
# BUILD = True

# █ Eichung des Elektromagneten

console.rule("Eichung des Elektromagneten")
I, B = np.genfromtxt('data/1_magnet.csv', delimiter=',', skip_header=1, unpack=True)
I *= ureg('A')
B *= ureg('mT')

fit_params = tools.pint_polyfit(I, B, 3)
print(f"{fit_params=}")


def calc_B(I):
    a, b, c, d = fit_params
    return a*I**3 + b*I**2 + c*I + d


with tools.plot_context(plt, 'A', 'mT', 'I', 'B') as plt2:
    plt2.plot(I, B, 'x', zorder=5, label='Messwerte')
    plt2.plot(I, tools.nominal_values(calc_B(I)), label='Regressionsgerade')
plt.grid()
plt.legend()
plt.tight_layout()
# if tools.BUILD:
plt.savefig('build/plt/1_magnet.pdf')


# █ Bestimmung der Landé-Faktoren

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
        'name': 'rot',
        'images': [
            {
                'I': ureg('0 A'),
                'path': 'Bilder/rot neu/IMG_0001.JPG',
                'find_peaks': {
                    'min_height': 0.6,  # entspricht der alten Analyse
                },
            },
            {
                'I': ureg('8 A'),
                'path': 'Bilder/rot neu/IMG_0002.JPG',
                # 'min_distance': 70,
                # 'min_height': 1600/2248,
                #
                'find_peaks': {
                    'min_distance': 1,
                    'min_height': 0.4,
                    'prominence': 0.2,
                },
            },
        ],
    },
    {
        'color': 'blau',
        'λ': ureg('480.0 nm'),
        'n': 1.4635,
        'rotate_deg': -2.3,
        'g_lit': 1.75,  # Mampfzwerg
        'name': 'blau_sigma',  # TODO
        'images': [
            {
                'I': ureg('0 A'),
                'polarisation': 0,  # in °
                'path': 'Bilder/blau/3_4.6A_sigma/IMG_0027.JPG',
                'find_peaks': {
                    'min_distance': 40,
                }
            },
            {
                'I': ureg('4.6 A'),
                'polarisation': 0,  # in °
                'path': 'Bilder/blau/3_4.6A_sigma/IMG_0028.JPG',
                'find_peaks': {
                    'min_distance': 1,
                    'min_height': 0.6,
                    'prominence': 0.0,
                },
            },
        ],
    },
    {
        'color': 'blau',
        'λ': ureg('480.0 nm'),
        'n': 1.4635,
        'rotate_deg': -2.3,
        'g_lit': 0.5,  # Mampfzwerg
        'name': 'blau_pi',
        'images': [
            {
                'I': ureg('0 A'),
                'polarisation': 0,  # TODO
                'path': 'Bilder/blau/5_8A_pi/IMG_0033.JPG',
                'find_peaks': {
                    'min_distance': 40,
                }
            },
            {
                'I': ureg('8 A'),  # TODO
                'polarisation': 0,  # TODO
                'path': 'Bilder/blau/5_8A_pi/IMG_0034.JPG',
                'find_peaks': {
                    'min_distance': 1,
                    'min_height': 0.6,
                    'prominence': 0.0,
                },
            },
        ],
    },
]


MESSREIHEN = [
    # FOO[0],
    FOO[1],
    FOO[2],
]


def get_δs_clean(img, opts):
    """Analyse mit B-Feld/Aufspaltung"""
    peaks2 = bildanalyse.get_peaks(
        img2,
        **opts,
        name=f"{messreihe['name']}_2",
        show=False,
    )

    # group peaks into pairs of 2
    pairs = list(itertools.pairwise(peaks2))[::2]
    print(f"{pairs=}")

    diffs = [b - a for a, b in pairs]
    print(f"{diffs=}")

    δs = np.array(diffs).mean()
    return δs


for messreihe in MESSREIHEN:
    console.rule(f"Messreihe {messreihe['name']}")

    # Bilder vorbereiten
    img1 = bildanalyse.preprocess_image(messreihe['images'][0]['path'], rotate_deg=messreihe['rotate_deg'])
    img2 = bildanalyse.preprocess_image(messreihe['images'][1]['path'], rotate_deg=messreihe['rotate_deg'])

    # Peaks finden
    peaks1 = bildanalyse.get_peaks(
        img1,
        **messreihe['images'][0]['find_peaks'],
        name=f"{messreihe['name']}_1",
        show=False,
    )

    # Analyse ohne B-Feld/Aufspaltung
    Δs = np.diff(peaks1).mean()
    # Δs /= 2  # !?

    # Analyse mit B-Feld/Aufspaltung
    δs = get_δs_clean(img2, messreihe['images'][1]['find_peaks'])

    print(f"Δs={Δs:.1f}")
    print(f"δs={δs:.1f}")

    Δλ_D = messreihe['λ']**2 / (2*d * np.sqrt(messreihe['n']**2 - 1))
    Δλ_D.ito('pm')
    print(f"Δλ_D={Δλ_D:.3f}")

    δλ = δs * Δλ_D / (2 * Δs)
    print(f"{δλ.mean()=}")

    # █ Bestimmung der Landé-Faktoren
    B = calc_B(messreihe['images'][1]['I'])
    print(f"B={B:.2f}")

    # g_ij = m_j * g_j - m_i * g_i
    μ_B = ureg.e * ureg.hbar / (2 * ureg.m_e)
    g_ij = δλ.mean() * ureg.h * ureg.c / (messreihe['λ']**2 * μ_B * B)  # Landé-Faktor
    g_ij.ito('dimensionless')
    print(f"{g_ij=}")

    print(tools.fmt_compare_to_ref(g_ij, messreihe['g_lit']))
