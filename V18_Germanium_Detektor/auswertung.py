# %%
import datetime
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
peaks = list(sorted(peaks))

assert len(peaks) == 11, f"Expected 11 peaks, found {len(peaks)}: {peaks}"  # for testing (@Mampfzwerg)
peaks

# %%


def fit_fn(x, a, x_0, sigma, c):
    return a * np.exp(-((x - x_0) ** 2) / (2 * sigma ** 2)) + c


def fit_peak(peak, plot=True):
    print(f"Peak bei ({peak}, {N_calib[peak].m})")
    AREA_HALFWIDTH = 40
    area_x = x_calib[(peak - AREA_HALFWIDTH): (peak + AREA_HALFWIDTH)]
    area_N = N_calib[(peak - AREA_HALFWIDTH): (peak + AREA_HALFWIDTH)]

    # ▒ Fitte Gauß-Kurve
    a, x_0, sigma, c = tools.curve_fit(fit_fn,
                                       area_x, area_N.m,
                                       p0=[max(area_N), peak, 1, min(area_N)],
                                       )
    print(f"→ a={a}, x_0={x_0}, sigma={sigma}, c={c}")

    if plot:
        # ▒ Plotte Peak
        plt.plot(area_x, area_N)
        plt.plot(area_x, fit_fn(
            area_x, *[param.n for param in [a, x_0, sigma, c]]
        ))
        plt.show()

    # return a, x_0, sigma, c
    # return x_0
    return {
        'a': a,
        'x_0': x_0,
        'sigma': sigma,
        'c': c,
    }


peak_fits = [fit_peak(peak, plot=True) for peak in peaks]
peak_channels = [fit['x_0'] for fit in peak_fits]


# %%
# Lade Emissionsspektrum aus der Literatur
# http://www.lnhb.fr/Laraweb/index.php
df = pd.read_csv("data/Eu-152.lara.txt", sep=" ; ", skiprows=13, skipfooter=1, engine='python', index_col=False)
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
mask = df['Energy (keV)'] > 100  # ignore that first peak…
lit_energies = df['Energy (keV)'][mask].values[:len(peak_channels)]
lit_intensities = df['Intensity (%)'][mask].values[:len(peak_channels)]
# ↓ all values
lit_energies_all = df['Energy (keV)'].values
lit_intensities_all = df['Intensity (%)'].values

# sort by energy
lit_energies, lit_intensities = zip(*sorted(zip(lit_energies, lit_intensities)))

lit_energies *= ureg.keV
lit_energies_all *= ureg.keV

peak_channels_n = tools.nominal_values(peak_channels * ureg.dimensionless)
m, b = tools.linregress(peak_channels_n, lit_energies)
# m, b = ureg('0.403 keV'), ufloat(-2.68, 0.05) * ureg.keV # @Mampfzwerg
# m, b = ureg('0.35 keV'), ufloat(-200, 0.05) * ureg.keV  # custom

# energy = m*channel + b
# channel = (energy - b)/m
lit_channels = tools.nominal_values((lit_energies - b) / m)
lit_channels_all = tools.nominal_values((lit_energies_all - b) / m)

# Plot
plt.figure()
with tools.plot_context(plt, 'dimensionless', 'keV', "Kanal", "Energie") as plt2:
    plt2.plot(peak_channels_n, m * peak_channels_n + b, label="Regressionsgerade")
    plt2.plot(peak_channels_n, lit_energies, 'x', label="Literaturwerte")
# plt.show()
plt.legend()
# %%

# ▒ Plotte Spektrum
plt.figure()
with tools.plot_context(plt, 'dimensionless', 'dimensionless', "Kanal", "Counts") as plt2:
    plt.plot(N_calib, label="Messwerte")
    plt.plot(peaks, N_calib[peaks], 'x', label="Peaks")
    # for lit_channel, lit_intensity in zip(lit_channels, lit_intensities):
    for lit_channel, lit_intensity in zip(lit_channels_all, lit_intensities_all):
        # TODO: Label
        plt.axvline(lit_channel, color='C2', alpha=min(1, lit_intensity/100*3))
# plt.xlim(right=4000)
plt.yscale('log')
plt.legend()
plt.show()

# %% ███ Effizienz ███
# TODO: Aus .Spe-Datei lesen
t_mess = ureg('3388 s')  # Messzeit

# TODO: Quelle?
r = ureg('2.25 cm')  # Radius
l = ureg('9.5 cm')  # Abstand Probe–Detektor
Ω = 2*np.pi*(1 - l/np.sqrt(l**2 + r**2))
print(f"Ω={Ω.to('pi'):.4f}")

probe_creationdate = datetime.datetime(2000, 10, 1, 0, 0, 0)
durchfuehrung_date = datetime.datetime(2022, 11, 28, 0, 0, 0)
age_probe = (durchfuehrung_date - probe_creationdate).total_seconds() * ureg.s
print(f"Alter Probe: {age_probe.to('s'):.2e} = {age_probe.to('a'):.2f}")

t_hw = ureg('4943 d')  # TODO: Unsicherheit
A_0 = ufloat(4130, 60) * ureg.Bq
A = A_0 * np.exp(-np.log(2) / t_hw * age_probe)
W = 1  # Emissionswahrscheinlichkeit; TODO
# W = df['Intensity (%)'] / 100 # Emissionswahrscheinlichkeit; TODO: untested
print(f"A={A.to('Bq'):.2f}")


# %% Flächeninhalt der Gaußkurven
a_np = np.array([fit['a'].n for fit in peak_fits])  # Amplituden
sigma_np = np.array([fit['sigma'].n for fit in peak_fits])  # Breiten (Standardabweichungen)
Z = a_np * np.sqrt(2*np.pi * sigma_np**2)  # Flächeninhalte
Z

# %% Effizienz
Q = 4*np.pi*Z / (Ω * A * W * t_mess)  # Effizienz
# assert Q.check('dimensionless'), "Q is not dimensionless"
Q.ito('dimensionless')  # NOTE: Q.check raises a KeyError for some reason, doing this instead → create an Issue?
Q

# %% Fit: Q(E) – Effizienz in Abhängigkeit von der Energie


def fit_fn_Q(E, Q_max, exponent):
    # NOTE: E / (1 keV)
    if isinstance(E, ureg.Quantity):
        E = E.to('keV').magnitude
    return Q_max * E**exponent


Q_max, exponent = tools.pint_curve_fit(fit_fn_Q,
                                       lit_energies, Q,
                                       (ureg.dimensionless, ureg.dimensionless),
                                       p0=(15 * ureg.dimensionless, -1 * ureg.dimensionless),
                                       return_p0=True,  # TODO
                                       )
print(f"Q_max={Q_max:.2f}")
print(f"exponent={exponent:.2f}")

# %% Plot: Q(E) – Effizienz in Abhängigkeit von der Energie
energy_linspace = tools.linspace(*tools.bounds(lit_energies), 100)
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', "Energie", "Effizienz") as plt2:
    plt2.plot(lit_energies, Q, fmt='x', label="Messwerte")
    plt2.plot(energy_linspace, fit_fn_Q(energy_linspace, Q_max, exponent), label="Fit")
plt.legend()
plt.show()
# %%
