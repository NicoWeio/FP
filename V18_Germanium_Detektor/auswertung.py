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

import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.define('percent = 1 / 100 = %')

console = Console()


@ureg.wraps(ureg.s**(-1), [ureg.dimensionless, ureg.s])
def calc_I(N, T):
    """
    Zählrate aus Anzahl N und Zeit T.
    """
    I = unp.uarray(N / T, np.sqrt(N)/T) / ureg.s
    return I


def load_spe(filename):
    # COULDDO: move to tools.py
    N = np.genfromtxt(filename, skip_header=12, skip_footer=14)
    assert len(N) == 8192
    x = np.arange(0, len(N))  # * ureg.dimensionless

    # TODO: Irgendwo den Untergrund abziehen!

    N *= ureg.dimensionless

    return x, N


def load_lara(filename):
    # http://www.lnhb.fr/Laraweb/index.php
    df = pd.read_csv(filename, sep=" ; ", skiprows=13, skipfooter=1, engine='python', index_col=False)
    # Filtern nach Zerfallsart (TODO untested)
    decay_mode_mask = df['Type'] == 'g'
    print(f"Filtered {len(df)} → {len(df[decay_mode_mask])} rows by decay mode")
    df = df[decay_mode_mask]
    return df


# █ Energiekalibrierung
x_calib, N_calib = load_spe("data/2022-11-28/1_Eu.Spe")

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


def fit_fn_peak(x, a, x_0, sigma, N_0):
    return a * np.exp(-((x - x_0) ** 2) / (2 * sigma ** 2)) + N_0


def fit_peak(peak, x, N, plot=True):
    # TODO: Stelle sicher, dass x Kanäle, nicht Energien sind.
    # assert not isinstance(x, ureg.Quantity) or x.units == ureg.dimensionless

    # print(f"Peak bei ({peak}, {N[peak].m})")
    AREA_HALFWIDTH = 40  # Radius des Bereichs, in dem die Gauß-Kurve gefittet wird
    area_x = x[(peak - AREA_HALFWIDTH): (peak + AREA_HALFWIDTH)]
    area_N = N[(peak - AREA_HALFWIDTH): (peak + AREA_HALFWIDTH)]

    # ▒ Fitte Gauß-Kurve
    a, x_0, σ, N_0 = tools.curve_fit(fit_fn_peak,
                                       area_x, area_N.m,
                                       p0=[max(area_N), peak, 1, min(area_N)],
                                       )
    # print(f"→ a={a}, x_0={x_0}, sigma={sigma}, N_0={N_0}")

    if plot:
        # ▒ Plotte Peak
        plt.plot(area_x, area_N)
        plt.plot(area_x, fit_fn_peak(
            area_x, *[param.n for param in [a, x_0, σ, N_0]]
        ))
        plt.show()

    # return a, x_0, sigma, c
    # return x_0
    return {
        'a': a,
        'x_0': x_0,
        'σ': σ,
        'N_0': N_0,
    }


peak_fits = [fit_peak(peak, x_calib, N_calib, plot=False) for peak in peaks]
peak_channels = [fit['x_0'] for fit in peak_fits]


# %%
# Lade Emissionsspektrum aus der Literatur
df = load_lara("data/emissions/Eu-152.lara.txt")
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
mask = df['Energy (keV)'] > 100  # ignore that first peak… # TODO: Might be obsolete since we're ignoring non-gamma decays now
lit_energies = df['Energy (keV)'][mask].values[:len(peak_channels)]
lit_intensities = df['Intensity (%)'][mask].values[:len(peak_channels)]
# ↓ all values
lit_energies_all = df['Energy (keV)'].values
lit_intensities_all = df['Intensity (%)'].values

# sort by energy
lit_energies, lit_intensities = zip(*sorted(zip(lit_energies, lit_intensities)))

lit_energies *= ureg.keV
lit_energies_all *= ureg.keV
lit_intensities *= ureg.percent
lit_intensities_all *= ureg.percent

peak_channels_n = tools.nominal_values(peak_channels * ureg.dimensionless)
m, n = tools.linregress(peak_channels_n, lit_energies)
print(f"m={m}, n={n}")

# energy = m*channel + b
# channel = (energy - b)/m
lit_channels = tools.nominal_values((lit_energies - n) / m)
lit_channels_all = tools.nominal_values((lit_energies_all - n) / m)

# Plot: Energie(Kanal)
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'keV', r"\text{Kanal}", "E") as plt2:
# TODO: \text funzt nicht ohne BUILD…
with tools.plot_context(plt, 'dimensionless', 'keV', "x", "E") as plt2:
    plt2.plot(peak_channels_n, m * peak_channels_n + n, label="Regressionsgerade")
    plt2.plot(peak_channels_n, lit_energies, 'x', zorder=5, label="Literaturwerte")
plt.legend()
if tools.BUILD:
    plt.savefig("build/plt/energy_calibration.pdf")
# plt.show()
# %%

# ▒ Plotte Spektrum
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'dimensionless', r"\text{Kanal}", "N") as plt2:
# TODO: \text funzt nicht ohne BUILD…
with tools.plot_context(plt, 'dimensionless', 'dimensionless', "x", "N") as plt2:
    plt.plot(N_calib, label="Messwerte")
    plt.plot(peaks, N_calib[peaks], 'x', label="Peaks")
    # for lit_channel, lit_intensity in zip(lit_channels, lit_intensities):
    for lit_channel, lit_intensity in zip(lit_channels_all, lit_intensities_all):
        # TODO: Label
        plt.axvline(lit_channel, color='C2', alpha=min(1, lit_intensity.to('dimensionless').m*3))
# plt.xlim(right=4000)
plt.yscale('log')
plt.legend()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_152-Eu.pdf")
plt.show()

# %%
# NOTE: wird später noch gebraucht
a_np = np.array([fit['a'].n for fit in peak_fits]) * ureg.dimensionless  # Amplituden
x_0_np = np.array([fit['x_0'].n for fit in peak_fits]) * ureg.dimensionless  # Mittelwerte
σ_np = np.array([fit['σ'].n for fit in peak_fits]) * ureg.dimensionless  # Breiten (Standardabweichungen)
N_0_np = np.array([fit['N_0'].n for fit in peak_fits]) * ureg.dimensionless  # Hintergründe

# %% █ Tabelle generieren
# TODO: Unsicherheiten
generate_table.generate_table_pint(
    "build/tab/1_energiekalibrierung.tex",
    (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
    (r"x_0", ureg.dimensionless, x_0_np, 1),
    (r"a", ureg.dimensionless, a_np, 1),
    (r"\sigma", ureg.dimensionless, σ_np, 2),
    (r"N_0", ureg.dimensionless, N_0_np, 2),
)

# %% ███ Effizienz ███
console.rule("Effizienz")
# TODO: Aus .Spe-Datei lesen
t_mess = ureg('3388 s')  # Messzeit

r = ureg('45 mm') / 2  # Radius [Versuchsanleitung]
l = ureg('9.5 cm')  # Abstand Probe–Detektor # TODO: Quelle?

Ω = 2*np.pi*(1 - l/np.sqrt(l**2 + r**2))
# print(f"Ω={Ω.to('pi'):.4f}")
print(f"Ω = {Ω:.4f} = 4π · {(Ω / (4*np.pi)).to('dimensionless'):.4f}")


probe_creationdate = datetime.datetime(2000, 10, 1, 0, 0, 0)  # [Versuchsanleitung]
durchfuehrung_date = datetime.datetime(2022, 11, 28, 0, 0, 0)
age_probe = (durchfuehrung_date - probe_creationdate).total_seconds() * ureg.s
print(f"Alter Probe: {age_probe.to('s'):.2e} = {age_probe.to('a'):.2f}")

t_hw = ureg('4943 d')  # TODO: Unsicherheit
A_0 = ufloat(4130, 60) * ureg.Bq  # Aktivität am Tag der Herstellung [Versuchsanleitung]
A = A_0 * np.exp(-np.log(2) / t_hw * age_probe)  # Aktivität am Tag der Durchführung
print(f"A={A.to('Bq'):.2f}")

# %% Flächeninhalte der Gaußkurven
# TODO: Richtige Faktoren (sigma etc.)?
Z = a_np * np.sqrt(2*np.pi * σ_np**2)  # Flächeninhalte
Z

# %% Effizienz
Q = 4*np.pi*Z / (Ω * A * lit_intensities * t_mess)  # Effizienz
# assert Q.check('dimensionless'), "Q is not dimensionless"
Q.ito('dimensionless')  # NOTE: Q.check raises a KeyError for some reason, doing this instead → create an Issue?
Q

# %% Fit: Q(E) – Effizienz in Abhängigkeit von der Energie


def fit_fn_Q(E, Q_max, exponent):
    # NOTE: E / (1 keV)
    if isinstance(E, ureg.Quantity):
        E = E.to('keV').magnitude
    return Q_max * E**exponent


Q_max, exponent = tools.pint_curve_fit(
    fit_fn_Q,
                                       lit_energies, Q,
                                       (ureg.dimensionless, ureg.dimensionless),
    p0=(50 * ureg.dimensionless, -1 * ureg.dimensionless),
                                       )
print(f"Q_max={Q_max:.2f}")
print(f"exponent={exponent:.2f}")

# %% Plot: Q(E) – Effizienz in Abhängigkeit von der Energie
energy_linspace = tools.linspace(*tools.bounds(lit_energies), 100)
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', "E", "Q") as plt2:
    plt2.plot(lit_energies, Q, fmt='x', label="Messwerte")
    plt2.plot(energy_linspace, fit_fn_Q(energy_linspace, Q_max, exponent), label="Fit")
plt.legend()
if tools.BUILD:
    plt.savefig("build/plt/effizienz.pdf")
plt.show()

# %% Tabelle generieren
# TODO: Unsicherheiten
generate_table.generate_table_pint(
    "build/tab/2_effizienz.tex",
    (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
    (r"W", ureg.percent, lit_intensities),
    (r"x_0", ureg.dimensionless, x_0_np, 1),  # redundant s.o.
    (r"\sigma", ureg.dimensionless, σ_np, 2),  # redundant s.o.
    (r"Z", ureg.dimensionless, Z, 0),
    (r"Q", ureg.dimensionless, tools.nominal_values(Q), 4),
)


# %% ███ Spektrum von 137Cs ███
console.rule("Spektrum von 137Cs")
# TODO: In mehrere .py-Dateien aufteilen

x, N = load_spe("data/2022-11-28/2_Cs.Spe")

E = m * x + n

# Peaks
# peaks, _ = find_peaks(N, height=100, distance=100)
peaks, _ = find_peaks(N, height=50, distance=1000)

assert len(peaks) == 2, f"Es sollten 2 Peaks (Rückstreupeak und Photopeak) gefunden werden. Gefunden wurden {len(peaks)} Peaks."
rueckstreupeak, photopeak = peaks
E_rueckstreupeak, E_photopeak = E[peaks]

# TODO: Statt plump die Maxima zu nehmen, Gaußkurven fitten!

print(f"Rückstreupeak: {E_rueckstreupeak}")
print(f"Photopeak: {E_photopeak}")

# %% Plot: N(E) – Spektrum von 137Cs
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', r"E_\gamma", r"N") as plt2:
    plt2.plot(E, N, label="Messwerte")  # fmt='x'
    plt2.plot(E[peaks], N[peaks], fmt='x', label="Peaks")
    plt2.plot(E, np.convolve(N.m, np.ones(20)/20, mode='same'), fmt='-', label="Smoothed")  # TODO: Gut so?
plt.yscale('log')
plt.xlim(right=800)
plt.legend()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_137-Cs.pdf")
plt.show()

# %% Fit: Photopeak
photopeak_fit = fit_peak(photopeak, x, N, plot=True)

# %%
