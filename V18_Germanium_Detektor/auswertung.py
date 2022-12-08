# %%
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
x_calib = np.arange(0, len(N_calib))  # * ureg.dimensionless

# %% finde Peaks
peaks, _ = find_peaks(
    N_calib,
    prominence=28,
    # prominence=70,
    # height=100,
)

peaks_to_ignore = {
    58, 62, 150, 200, 226, 415,
    536, 576,
    # 596, # erster Peak?
    764,
}
peaks = [peak for peak in peaks if peak not in peaks_to_ignore]

peaks.append(4200)  # TODO!

assert len(peaks) == 11, f"Expected 11 peaks, found {len(peaks)}: {peaks}"  # for testing (@Mampfzwerg)
peaks

# %%


def fit_fn(x, a, x_0, sigma, c):
    return a * np.exp(-((x - x_0) ** 2) / (2 * sigma ** 2)) + c


def fit_peak(peak, plot=True):
    print(f"Peak {peak} with {N_calib[peak]}")
    area_halfwidth = 10
    area_x = x_calib[(peak - area_halfwidth): (peak + area_halfwidth)]
    area_N = N_calib[(peak - area_halfwidth): (peak + area_halfwidth)]

    # ▒ Fitte Gauß-Kurve
    a, x_0, sigma, c = tools.curve_fit(fit_fn,
                                       area_x, area_N.m,
                                       p0=[max(area_N), peak, 1, min(area_N)],
                                       )
    print(f"a={a}, x_0={x_0}, sigma={sigma}, c={c}")

    if plot:
        # ▒ Plotte Peak
        plt.plot(area_x, area_N)
        plt.plot(area_x, fit_fn(
            area_x, *[param.n for param in [a, x_0, sigma, c]]
        ))
        plt.show()

    # return a, x_0, sigma, c
    return x_0


peak_channels = [fit_peak(peak, plot=False) for peak in peaks]


# %%
# Lade Emissionsspektrum aus der Literatur
# http://www.lnhb.fr/Laraweb/index.php
df = pd.read_csv("Eu-152.lara.txt", sep=" ; ", skiprows=13, skipfooter=1, engine='python', index_col=False)
df

# %% Lineare Regression (Zuordnung Energie ↔ Kanal)
# select the most intense energies from the df, so that len(literature_energies) == len(peak_positions)
# assumes the peaks are already sorted by intensity

# TODO redo:
# - sort(ed) by intensity
# - ignore peaks with E < 100 keV
# - sort by energy
# - cut to len(peak_channels)

# ↓ selected values
# mask = df['Energy (keV)'] > 100  # ignore that first peak…
lit_energies = df['Energy (keV)'][mask].values[:len(peak_channels)] * ureg.keV
lit_intensities = df['Intensity (%)'][mask].values[:len(peak_channels)]
# ↓ all values
lit_energies_all = df['Energy (keV)'].values * ureg.keV
lit_intensities_all = df['Intensity (%)'].values

# sort by energy

peak_channels_n = tools.nominal_values(peak_channels * ureg.dimensionless)
# m, b = tools.linregress(peak_channels_n, lit_energies)
# m, b = ureg('0.403 keV'), ufloat(-2.68, 0.05) * ureg.keV # @Mampfzwerg
m, b = ureg('0.35 keV'), ufloat(-200, 0.05) * ureg.keV  # custom

# energy = m*channel + b
# channel = (energy - b)/m
lit_channels = tools.nominal_values((lit_energies - b) / m)
lit_channels_all = tools.nominal_values((lit_energies_all - b) / m)

# Plot
plt.plot(peak_channels_n, lit_energies, 'x')
plt.plot(peak_channels_n, tools.nominal_values(m * peak_channels_n + b))
plt.xlabel("Kanal")
plt.ylabel("Energie [keV]")
# plt.show()
# %%

# ▒ Plotte Spektrum
plt.figure()
plt.plot(N_calib)
plt.plot(peaks, N_calib[peaks], "x")
for lit_channel, lit_intensity in zip(lit_channels, lit_intensities):
    # for lit_channel, lit_intensity in zip(lit_channels_all, lit_intensities_all):
    plt.axvline(lit_channel, color='red', alpha=min(1, lit_intensity/100*3))
plt.xlabel("Kanal")
plt.ylabel("Counts")
# plt.xlim(right=4000)
# plt.yscale('log')
plt.show()

# %%
