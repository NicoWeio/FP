import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pint
import rich
import uncertainties.unumpy as unp
from numpy.linalg import inv
from rich.console import Console
from scipy.signal import find_peaks
from uncertainties import ufloat

# import generate_table
import tools

ureg = pint.UnitRegistry()
console = Console()


@ureg.wraps(ureg.s**(-1), [ureg.dimensionless, ureg.s])
def calc_I(N, T):
    """
    Zählrate aus Anzahl N und Zeit T.
    """
    I = unp.uarray(N / T, np.sqrt(N)/T) / ureg.s
    return I


def read_spe(filename):
    # TODO: move to tools.py
    N = np.genfromtxt(filename, skip_header=12, skip_footer=14)
    assert len(N) == 8192

    N *= ureg.dimensionless
    return N


# █ Energiekalibrierung
N_calib = read_spe("data/2022-11-28/1_Eu.Spe")
x_calib = np.arange(0, len(N_calib)) # * ureg.dimensionless

# finde Peaks
peaks, _ = find_peaks(N_calib, prominence=50)

# ▒ Plotte Spektrum
plt.plot(N_calib)
plt.plot(peaks, N_calib[peaks], "x")
plt.xlabel("Kanal")
plt.ylabel("Counts")
# plt.yscale('log')
# plt.show()

def fit_fn(x, a, x_0, sigma, c):
    return a * np.exp(-((x - x_0) ** 2) / (2 * sigma ** 2)) + c

for peak in peaks:
    print(f"Peak {peak} with {N_calib[peak]}")
    area_halfwidth = 10
    area_x = x_calib[(peak - area_halfwidth) : (peak + area_halfwidth)]
    area_N = N_calib[(peak - area_halfwidth) : (peak + area_halfwidth)]

    # ▒ Fitte Gauß-Kurve
    a, x_0, sigma, c = tools.curve_fit(fit_fn,
    area_x, area_N.m,
    p0=[max(area_N), peak, 1, min(area_N)],
    )
    print(f"a={a}, x_0={x_0}, sigma={sigma}, c={c}")

    # ▒ Plotte Peak
    plt.plot(area_x, area_N)
    plt.plot(area_x, fit_fn(
        area_x, *[param.n for param in [a, x_0, sigma, c]]
    ))
    plt.show()
