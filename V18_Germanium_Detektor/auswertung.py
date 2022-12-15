# %%
import datetime
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pint
import uncertainties.unumpy as unp
from rich.console import Console
from scipy.signal import find_peaks
from uncertainties import ufloat, UFloat

import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.define('percent = 1 / 100 = %')

console = Console()


# @ureg.wraps(ureg.s**(-1), [ureg.dimensionless, ureg.s])
# def calc_I(N, T):
#     """
#     Zählrate aus Anzahl N und Zeit T.
#     """
#     I = unp.uarray(N / T, np.sqrt(N)/T) / ureg.s
#     return I


@ureg.check('[energy]')
def calc_ε(E):
    return E / (ureg.m_e * ureg.c**2)


def _load_spe(filename):
    # COULDDO: move to tools.py

    # Messzeit lesen
    lines = Path(filename).read_text().splitlines()
    assert lines[8] == "$MEAS_TIM:"
    T = int(lines[9].split(" ")[0]) * ureg.s

    N = np.genfromtxt(filename, skip_header=12, skip_footer=14)
    assert len(N) == 8192
    x = np.arange(0, len(N))  # * ureg.dimensionless

    N *= ureg.dimensionless

    return x, N, T


def load_spe(filename, subtract_background=True):
    """
    Lade eine Spe-Datei und subtrahiere den Untergrund.
    """
    x, N, T = _load_spe(filename)

    if subtract_background:
        x_bg, N_bg, T_bg = _load_spe("data/2022-11-28/5_Background.Spe")
        N -= N_bg * (T / T_bg)

        if not all(N >= 0):
            print("⚠ Durch Subtraktion des Untergrunds sind negative Werte entstanden!")

    return x, N, T


def load_lara(filename, parent=None):
    """
    Loads a LARAWEB file and returns the energy and intensity columns, limited to γ decays.
    Data source: http://www.lnhb.fr/Laraweb/index.php
    """
    assert filename.endswith('.lara.txt')
    # LARA files have a header (with varying length) and footer that we need to skip.
    # That is all this block does.
    contents = Path(filename).read_text()
    lines = np.array(contents.splitlines(), dtype=object)
    # find delimiter lines (starting with "---" or "===")
    delimiter_mask = np.array([
        line.startswith("---") or line.startswith("===")
        for line in lines
    ])
    # add the header to the mask
    delimiter_mask[:np.where(delimiter_mask)[0][0]] = True
    # get the lines without delimiters as a file-like object
    data_lines = lines[~delimiter_mask]
    data = StringIO("\n".join(data_lines))

    df = pd.read_csv(data, sep=" ; ", engine='python', index_col=False)

    # filter out non-gamma decays
    df = df[df['Type'] == 'g']

    if parent is not None:
        df = df[df['Parent'] == parent]

    assert len(df) > 0, f"No data found for {filename} (parent={parent})"

    lit_energies = df['Energy (keV)'].values * ureg.keV
    lit_intensities = df['Intensity (%)'].values * ureg.percent

    return lit_energies, lit_intensities


def n_most_intense(energies, intensities, n):
    # select the most intense energies, so that len(lit_energies) == n
    i_by_intensity = np.argsort(intensities)[::-1]  # sort by intensity (descending)
    filtered_energies = energies[i_by_intensity][:n]
    filtered_intensities = intensities[i_by_intensity][:n]

    # sort by energy again
    i_by_energy = np.argsort(filtered_energies)
    filtered_energies = filtered_energies[i_by_energy]
    filtered_intensities = filtered_intensities[i_by_energy]

    return filtered_energies, filtered_intensities


def intensity_to_alpha(intensity, exponent=0.25):
    assert isinstance(intensity, pint.Quantity)
    i = intensity.to('dimensionless').m
    return i**exponent


# %%
console.rule("0. Untergrund")
x_bg, N_bg, T_bg = load_spe("data/2022-11-28/5_Background.Spe", subtract_background=False)

# Plot: Untergrund
plt.figure()
with tools.plot_context(plt, 'dimensionless', 'dimensionless', "x", "N") as plt2:
    plt2.plot(x_bg, N_bg, label="Messwerte")
plt.yscale('log')
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/background.pdf")
# if tools.PLOTS:
#     plt.show()

# %% █ Energiekalibration
console.rule("1. Energiekalibration: Eu-152")
x_calib, N_calib, T_calib = load_spe("data/2022-11-28/1_Eu.Spe")

N_calib_smooth = np.convolve(N_calib.m, np.ones(50)/50, mode='same') * ureg.dimensionless

N_calib_for_peaks = N_calib_smooth.copy()
N_calib_for_peaks[:400] = 0  # ignore first channels

# finde Peaks
peaks, _ = find_peaks(
    # np.log(N_calib + 1),
    # prominence=3.2,
    # distance=100,

    np.log(N_calib_smooth + 1),
    # np.log(N_calib_for_peaks + 1),
    prominence=0.6,
    distance=100,
)
# COULDDO: mehr Peaks?

peaks_to_ignore = [
    79,
    5861,
    6287,
]
peaks = [p for p in peaks if p not in peaks_to_ignore]
peaks = list(sorted(peaks))

assert len(peaks) == 11, f"Expected 11 peaks, found {len(peaks)}: {peaks}"  # for testing (@Mampfzwerg)
peaks

# %%


def fit_fn_peak(x, a, x_0, σ, N_0):
    """
    x: Kanal
    ---
    a: Amplitude
    x_0: Mittelwert
    σ: Breite (Standardabweichung)
    N_0: Hintergrund
    """
    return a * np.exp(-((x - x_0)**2) / (2 * σ**2)) + N_0


def fit_peak(peak_seed, x, N, fit_radius=40, plot=True, plot_path=None, channel_to_E=None):
    """
    Helfer-Funktion, um eine Gauß-Kurve auf einen Peak zu fitten.

    peak: Kanal, in dem der Peak etwa liegt (gefunden z.B. mit scipy.signal.find_peaks)
    x: Kanäle
    N: Zählraten
    fit_radius: Radius des Bereichs, in dem die Gauß-Kurve gefittet wird (in Kanälen)
    """
    # TODO: Stelle sicher, dass x Kanäle, nicht Energien sind.
    # assert not isinstance(x, ureg.Quantity) or x.units == ureg.dimensionless

    # print(f"Peak bei ({peak}, {N[peak].m})")
    range_x = x[(peak_seed - fit_radius): (peak_seed + fit_radius)]
    range_N = N[(peak_seed - fit_radius): (peak_seed + fit_radius)]

    # wahres Maximum im Suchbereich
    true_peak = np.argmax(range_N) + peak_seed - fit_radius

    # ▒ Fitte Gauß-Kurve
    # COULDDO: tools.pint_curve_fit
    a, x_0, σ, N_0 = tools.curve_fit(
        fit_fn_peak,
        range_x, range_N.m,
        p0=[max(range_N), true_peak, 1, min(range_N)],
    )

    # Sanity checks
    try:
        assert a > 0, f"Amplitude a={a} should be positive"
        assert x_0 > 0, f"x_0={x_0} should be positive"
        assert σ > 0, f"σ={σ} should be positive"
        assert N_0 > 0, f"N_0={N_0} should be positive"
    except AssertionError as e:
        print(e)
        print("Continuing for now, but this should be fixed.")

    # FWHM / FWTM
    fwhm = 2 * σ * np.sqrt(2 * np.log(2))
    fwhm_height = a / 2 + N_0
    fwtm = 2 * σ * np.sqrt(2 * np.log(10))
    fwtm_height = a / 10 + N_0

    # (Flächen-)Inhalt
    Z = a * unp.sqrt(2*np.pi * σ**2)  # TODO: Richtige Faktoren (sigma etc.)?

    if plot:
        plt.figure()
        with tools.plot_context(plt, 'dimensionless', 'dimensionless', r"xTODO", r"N") as plt2:
            plt2.plot(range_x, range_N, label="Messwerte")
            plt2.plot(
                range_x,
                fit_fn_peak(range_x, *[param.n for param in [a, x_0, σ, N_0]]),
                label="Fit",
            )

            plt2.plot([x_0 - fwhm / 2, x_0 + fwhm / 2], [fwhm_height, fwhm_height], show_yerr=False, label="FWHM")
            plt2.plot([x_0 - fwtm / 2, x_0 + fwtm / 2], [fwtm_height, fwtm_height], show_yerr=False, label="FWTM")

        plt.xlim(range_x[0], range_x[-1])
        plt.legend()
        plt.tight_layout()
        if plot_path:
            plt.savefig(plot_path)
        if tools.PLOTS:
            plt.show()

    data = {
        # ↓ Fit-Parameter
        'a': a,
        'x_0': x_0,
        'σ': σ,
        'N_0': N_0,
        # ↓ berechnete Werte
        'fwhm': fwhm,
        'fwtm': fwtm,
        'Z': Z,
        # ↓ andere
        'true_peak': true_peak,
    }

    if channel_to_E:
        data['E'] = channel_to_E(x_0)

    return data


def peak_fits_to_arrays(peak_fits):
    keys = peak_fits[0].keys()
    return {
        key: (
            # für 'E':
            tools.pintify([tools.nominal_value(fit[key]) for fit in peak_fits])
            if isinstance(peak_fits[0][key], pint.Quantity)
            # für alle anderen:
            else (
                np.array([fit[key].n for fit in peak_fits]) * ureg.dimensionless
                if isinstance(peak_fits[0][key], UFloat)
                else np.array([fit[key] for fit in peak_fits]) * ureg.dimensionless
            )
        )
        for key in keys
    }


peak_fits = [fit_peak(peak, x_calib, N_calib, plot=False) for peak in peaks]
peak_fit_arrays = peak_fits_to_arrays(peak_fits)
peak_channels = [fit['x_0'] for fit in peak_fits]  # COULDDO: Brauche ich nicht wirklich, oder?
peak_channels_n = tools.nominal_values(peak_channels * ureg.dimensionless)


# %%
lit_energies_all, lit_intensities_all = load_lara("data/emissions/Eu-152.lara.txt")
# select the most intense energies, so that len(lit_energies) == len(peak_channels)
lit_energies, lit_intensities = n_most_intense(lit_energies_all, lit_intensities_all, n=len(peak_channels))

# %% Lineare Regression (Zuordnung Energie ↔ Kanal)

m, n = tools.linregress(peak_channels_n, lit_energies)
print(f"m={m}, n={n}")


def channel_to_E(x):
    # COULDDO: warn about negative energies
    return m * x + n


def E_to_channel(E):
    return (E - n) / m


lit_channels = tools.nominal_values(E_to_channel(lit_energies))
lit_channels_all = tools.nominal_values(E_to_channel(lit_energies_all))

# Plot: Energie(Kanal)
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'keV', r"\text{Kanal}", "E") as plt2:
# TODO: \text funzt nicht ohne BUILD…
with tools.plot_context(plt, 'dimensionless', 'keV', "x", "E") as plt2:
    plt2.plot(peak_channels_n, channel_to_E(peak_channels_n), show_yerr=False, label="Regressionsgerade")
    plt2.plot(peak_channels_n, lit_energies, 'x', zorder=5, label="Literaturwerte")
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/energy_calibration.pdf")
# if tools.PLOTS:
#     plt.show()
# %%

# ▒ Plotte Spektrum
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'dimensionless', r"\text{Kanal}", "N") as plt2:
# TODO: \text funzt nicht ohne BUILD…
with tools.plot_context(plt, 'dimensionless', 'dimensionless', "x", "N") as plt2:
    plt.plot(N_calib, label="Messwerte")
    # plt.plot(x_calib, N_calib_smooth, '-', label="Messwerte (geglättet)")

    plt.plot(peak_fit_arrays['true_peak'], N_calib[peak_fit_arrays['true_peak']], 'x', label="Peaks")

    for lit_energy, lit_intensity in zip(lit_channels_all, lit_intensities_all):
        plt.axvline(lit_energy, color='C2', alpha=intensity_to_alpha(lit_intensity))

# add labels for axvlines
handles, labels = plt.gca().get_legend_handles_labels()
handles.append(
    mpatches.Patch(color='C2', label="Eu-152 (Literaturwerte)")
)

# plt.xlim(right=4000)
plt.ylim(bottom=0.5)
plt.yscale('log')
plt.legend(handles=handles)
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_152-Eu.pdf")
if tools.PLOTS:
    plt.show()

# %% █ Tabelle generieren
# COULDDO: Unsicherheiten
generate_table.generate_table_pint(
    "build/tab/1_energiekalibrierung.tex",
    (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
    (r"x_0", ureg.dimensionless, peak_fit_arrays['x_0'], 1),
    (r"a", ureg.dimensionless, peak_fit_arrays['a'], 1),
    (r"\sigma", ureg.dimensionless, peak_fit_arrays['σ'], 2),
    (r"N_0", ureg.dimensionless, peak_fit_arrays['N_0'], 2),
)

# %% ███ Effizienz ███
console.rule("1. Vollenergienachweiswahrscheinlichkeit / Effizienz: Eu-152")

r = ureg('45 mm') / 2  # Radius [Versuchsanleitung]
# Abstand Probe–Detektor =
#   Abstand Probe–Schutzhaube [eigene Messung]
# + Abstand Schutzhaube–Detektor [Versuchsanleitung]
l = ureg('7.0 cm') + ureg('1.5 cm')

Ω = 2*np.pi*(1 - l/np.sqrt(l**2 + r**2))
# print(f"Ω={Ω.to('pi'):.4f}")
print(f"Ω = {Ω:.4f} = 4π · {(Ω / (4*np.pi)).to('dimensionless'):.4f}")


probe_creationdate = datetime.datetime(2000, 10, 1, 0, 0, 0)  # [Versuchsanleitung]
durchfuehrung_date = datetime.datetime(2022, 11, 28, 0, 0, 0)
age_probe = (durchfuehrung_date - probe_creationdate).total_seconds() * ureg.s
print(f"Alter Probe: {age_probe.to('s'):.2e} = {age_probe.to('a'):.2f}")

# t_hw = ureg('4943 d')  # @Mampfzwerg
t_hw = ufloat(13.522, 0.016) * ureg.a  # [lara]
A_0 = ufloat(4130, 60) * ureg.Bq  # Aktivität am Tag der Herstellung [Versuchsanleitung]
A = A_0 * unp.exp((-np.log(2) / t_hw * age_probe).to('dimensionless').m)  # Aktivität am Tag der Durchführung
print(f"A={A.to('Bq'):.2f}")

# %% Flächeninhalte der Gaußkurven
Z = peak_fit_arrays['Z']

# %% Effizienz
Q = 4*np.pi*Z / (Ω * A * lit_intensities * T_calib)  # Effizienz
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
print(f"Q_max (aka p): {Q_max:.2f}")
print(f"exponent (aka q): {exponent:.2f}")

# %% Plot: Q(E) – Effizienz in Abhängigkeit von der Energie
energy_linspace = tools.linspace(*tools.bounds(lit_energies), 100)
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', "E", "Q") as plt2:
    plt2.plot(lit_energies, Q, fmt='x', label="Messwerte")
    # COULDDO: Unsicherheiten als Schattierung
    plt2.plot(energy_linspace, fit_fn_Q(energy_linspace, Q_max, exponent), show_yerr=False, label="Fit")
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/effizienz.pdf")
if tools.PLOTS:
    plt.show()

# %% Tabelle generieren
# COULDDO: Unsicherheiten
generate_table.generate_table_pint(
    "build/tab/2_effizienz.tex",
    (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
    (r"W", ureg.percent, lit_intensities),
    (r"Z", ureg.dimensionless, Z, 0),
    (r"Q", ureg.dimensionless, tools.nominal_values(Q), 4),
)


# %% ███ Spektrum von Cs-137 ███
console.rule("2. Spektrum von Cs-137")
# COULDDO: In mehrere .py-Dateien aufteilen

x, N, T = load_spe("data/2022-11-28/2_Cs.Spe")
x_pint = x * ureg.dimensionless

E = channel_to_E(x)  # NOTE: This means E has uncertainties!

# Peaks
peaks, _ = find_peaks(N, height=50, distance=1000)

assert len(peaks) == 2, f"Es sollten 2 Peaks (Rückstreupeak und Photopeak) gefunden werden. Gefunden wurden {len(peaks)} Peaks."
rueckstreupeak_raw, photopeak_raw = peaks

# %% Plot: N(E) – Spektrum von 137Cs
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', r"E_\gamma", r"N") as plt2:
    plt2.plot(E, N, label="Messwerte")  # fmt='x'
    # plt2.plot(E[peaks], N[peaks], 'x', show_xerr=False, label="Peaks")
    plt2.plot(E[rueckstreupeak_raw], N[rueckstreupeak_raw], 'x', show_xerr=False, label="Rückstreupeak")
    plt2.plot(E[photopeak_raw], N[photopeak_raw], 'x', show_xerr=False, label="Photopeak")
plt.ylim(bottom=0.5)
plt.yscale('log')
plt.xlim(right=800)
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_137-Cs.pdf")
if tools.PLOTS:
    plt.show()

# %% Fit: Rückstreupeak & Photopeak
# peak_fits = [fit_peak(peak, x, N, plot=True, channel_to_E=channel_to_E) for peak in peaks]
peak_fits = [fit_peak(peak, x, N, plot=True, plot_path=f"build/plt/3_gauss_{i}.pdf",
                      channel_to_E=channel_to_E, fit_radius=([150, 40][i])) for i, peak in enumerate(peaks)]
peak_fit_arrays = peak_fits_to_arrays(peak_fits)
rueckstreupeak_fit, photopeak_fit = peak_fits

E_photopeak_lit = ureg('661.657 keV')  # TODO: Quelle, Unsicherheit
print(tools.fmt_compare_to_ref(photopeak_fit['E'], E_photopeak_lit, name="Photopeak (Fit vs. Literatur)"))

E_fwhm_fit = channel_to_E(photopeak_fit['fwhm'])
E_fwtm_fit = channel_to_E(photopeak_fit['fwtm'])

print(f"Photopeak → FWHM: {E_fwhm_fit:.2f}")
print(f"Photopeak → FWTM: {E_fwtm_fit:.2f}")
print(f"Photopeak → FWHM/FWTM: {(E_fwhm_fit/E_fwtm_fit):.2f}")

# %% █ Tabelle generieren
# COULDDO: Unsicherheiten
generate_table.generate_table_pint(
    "build/tab/3_Cs-137.tex",
    # ↓ Fit-Parameter
    (r"x_0", ureg.dimensionless, peak_fit_arrays['x_0'], 1),
    (r"a", ureg.dimensionless, peak_fit_arrays['a'], 1),
    (r"\sigma", ureg.dimensionless, peak_fit_arrays['σ'], 2),
    (r"N_0", ureg.dimensionless, peak_fit_arrays['N_0'], 2),
    # ↓ berechnete Werte
    (r"E(x_0)", ureg.kiloelectron_volt, peak_fit_arrays['E']),
    (r"Z(a, \sigma)", ureg.dimensionless, peak_fit_arrays['Z'], 0),
)

# %% Compton-Kante


@ureg.check('[energy]')
def calc_compton_edge(E_photopeak):
    ε = calc_ε(E_photopeak)
    return E_photopeak * 2*ε / (1 + 2*ε)


E_compton_lit = calc_compton_edge(E_photopeak_lit)

# Fit links-rechts
# TODO: subject to confirmation bias
mask_l = (E < E_compton_lit) & (E > (E_compton_lit - 100*ureg.keV))
mask_r = (E > E_compton_lit) & (E < (E_compton_lit + 20*ureg.keV))
mask_lr = mask_l | mask_r

m_l, n_l = tools.linregress(tools.nominal_values(E[mask_l]), N[mask_l])
m_r, n_r = tools.linregress(tools.nominal_values(E[mask_r]), N[mask_r])

# compton peak is at the intersection of the two lines
E_compton_fit = (n_r - n_l) / (m_l - m_r)

print(tools.fmt_compare_to_ref(E_compton_fit, E_compton_lit, name="Compton-Kante (Fit vs. Literatur)"))


# %% Plot: Compton-Kante
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', r"E_\gamma", r"N") as plt2:
    plt2.plot(E[mask_l], N[mask_l], 'x', show_xerr=False, color='C0', alpha=0.5, label="Messwerte links")
    plt2.plot(E[mask_r], N[mask_r], 'x', show_xerr=False, color='C1', alpha=0.5, label="Messwerte rechts")
    # plt2.plot(E[peaks], N[peaks], fmt='x', label="Peaks")
    plt2.plot(E[mask_lr], m_l*E[mask_lr] + n_l, show_yerr=False, color='C0', label="Fit links")
    plt2.plot(E[mask_lr], m_r*E[mask_lr] + n_r, show_yerr=False, color='C1', label="Fit rechts")

    plt.axvline(tools.nominal_value(E_compton_fit).to('keV').m, color='C2', label="Compton-Kante (Fit)")
    plt.axvline(E_compton_lit.to('keV').m, color='C3', label="Compton-Kante (Literatur)")
plt.ylim(top=max(N[mask_lr]))
plt.yscale('linear')
plt.legend(loc='lower left')  # TODO
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/compton-kante.pdf")
# %%


def fit_fn_klein_nishina(E, a, N_0):
    # E: Energie des gestoßenen Elektrons

    E_γ = tools.nominal_value(photopeak_fit['E'])  # Energie des einfallenden Gammaquants (=hν)
    ε = calc_ε(E_γ)
    m_0 = ureg.m_e  # ?

    if not isinstance(E, ureg.Quantity):
        E *= ureg.keV

    dσ_div_dE = (
        a *
        (
            (np.pi * ureg.r_e**2) / (m_0 * ureg.c**2 * ε**2)
        ) *
        (
            2 +
            (
                E / (E_γ - E)
            )**2 *
            (
                ε**-2 +
                (E_γ - E) / E_γ -
                (2 * (E_γ - E)) / (ε * E_γ)
            )
        )
        # + N_0
    )
    assert dσ_div_dE.check('m²/keV')
    return dσ_div_dE.to('nm²/keV').m * ureg.dimensionless


# def fit_fn_klein_nishina_2(E, a, N_0):
#     return fit_fn_klein_nishina(
#         E, a, N_0,
#         E_γ=E_photopeak_fit,
#     )

# mask_compton = (E > E_rueckstreupeak) & (E < E_compton_fit)
mask_compton_fit = (E > ureg('350 keV')) & (E < E_compton_fit)
mask_compton_plot = (E > ureg('200 keV')) & (E < ureg('500 keV'))

a, N_0 = tools.pint_curve_fit(
    fit_fn_klein_nishina,
    tools.nominal_values(E[mask_compton_fit]), N[mask_compton_fit],
    (ureg.dimensionless, ureg.dimensionless),
    # p0=(5E14 * ureg.dimensionless, 1 * ureg.dimensionless),
    p0=(5E10 * ureg.dimensionless, 1 * ureg.dimensionless),
    # return_p0=True,
)

print(f"a: {a:.2f}")
print(f"N_0: {N_0:.2f}")

N_fit = fit_fn_klein_nishina(E, a, N_0)
# Der Inhalt ergibt sich direkt aus Summation,
# weil wir die Fläche mit Kanälen (statt Energien) auf der x-Achse betrachten…
Z = sum(N_fit[E < E_compton_lit])

print(f"→ Z: {Z:.1f}")

# %% Plot: Klein-Nishina
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', r"E_\gamma", r"N") as plt2:
    mask_compton_plotonly = mask_compton_plot & ~mask_compton_fit
    plt2.plot(E[mask_compton_plot], N[mask_compton_plot], 'x', show_xerr=False, color='C0', alpha=0.5, label="Messwerte")
    plt2.plot(E[mask_compton_plotonly], N[mask_compton_plotonly], 'x', show_xerr=False,
              color='C1', alpha=0.5, label="Messwerte (nicht berücksichtigt)")
    plt2.plot(E[mask_compton_plot], fit_fn_klein_nishina(E[mask_compton_plot], a, N_0),
              show_xerr=False, show_yerr=False, color='C2', label="Fit")
plt.yscale('linear')
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/klein-nishina.pdf")

# %%
lit_energies_all_Ba, intensities_all_Ba = load_lara("data/emissions/Ba-133.lara.txt")
lit_energies_all_Sb, intensities_all_Sb = load_lara("data/emissions/Sb-125.lara.txt")
# energies, intensities = n_most_intense(energies_all, intensities_all, n=len(peak_channels))

# %%
console.rule("3. Aktivitätsbestimmung: Ba-133")

x_Ba, N_Ba, T_Ba = load_spe("data/2022-11-28/3_Ba.Spe")
E_Ba = channel_to_E(x_Ba)

# ▒ Plotte Spektrum
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'dimensionless', r"\text{Kanal}", "N") as plt2:
# TODO: \text funzt nicht ohne BUILD…
with tools.plot_context(plt, 'keV', 'dimensionless', "E", "N") as plt2:
    plt2.plot(E_Ba, N_Ba, label="Messwerte")
    # plt.plot(peaks, N_Ba[peaks], 'x', label="Peaks")  # TODO: Peaks tatsächlich suchen…
    for lit_energy, lit_intensity in zip(lit_energies_all_Ba, intensities_all_Ba):
        plt.axvline(lit_energy.to('keV').m, color='C2', alpha=intensity_to_alpha(lit_intensity))
    for lit_energy, lit_intensity in zip(lit_energies_all_Sb, intensities_all_Sb):
        plt.axvline(lit_energy.to('keV').m, color='C3', alpha=intensity_to_alpha(lit_intensity))

# add labels for axvlines
handles, labels = plt.gca().get_legend_handles_labels()
handles += [
    mpatches.Patch(color='C2', label="Ba-133 (Literaturwerte)"),
    mpatches.Patch(color='C3', label="Sb-125 (Literaturwerte)"),
]

plt.xlim(
    left=0,  # COULDDO: Why is this necessary!?
    right=1.1 * max(lit_energies_all_Ba.max(), lit_energies_all_Sb.max()).to('keV').m,
)
plt.ylim(bottom=0.5)
plt.yscale('log')
plt.legend(handles=handles)
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_133-Ba.pdf")
if tools.PLOTS:
    plt.show()

# %%
console.rule("4. Nuklididentifikation und Aktivitätsbestimmung: Uran & Zerfallsprodukte")

x_U, N_U, T_U = load_spe("data/2022-11-28/4_U.Spe")
E_U = channel_to_E(x_U)

PARENTS = [
    "U-238",
    "Th-234",
    "Pa-234",
    "U-234",
    "Th-230",
    "Ra-226",
    "Rn-222",
    "Po-218",
    # "At-218", # keine Gamma-Zerfälle
    "Pb-214",
    "Rn-218",
    "Bi-214",
    "Po-214",
    "Tl-210",
    "Pb-210",
    "Bi-210",
    "Hg-206",
    "Po-210",
    "Tl-206",
]

data_dict = {
    parent:
    load_lara("data/emissions/U-238.lara.txt", parent=parent)
    for parent in PARENTS
}

# %% rank parents by their maximum intensity
parent_to_max_intensity = [(parent, max(intensities)) for parent, (energies, intensities) in data_dict.items()]
parent_to_max_intensity.sort(key=lambda x: x[1], reverse=True)
for parent, max_intensity in parent_to_max_intensity:
    print(f"{parent}: {max_intensity:.2f}")

# filter out parents with too low (max) intensity
# data_dict_filtered = {
#     parent: (energies, intensities)
#     for parent, (energies, intensities) in data_dict.items()
#     if max(intensities) > 0.1
# }
data_dict_filtered = {
    parent: data_dict[parent]
    for parent, max_intensity in parent_to_max_intensity
    if max_intensity > 0.1
}

print(f"Filtered out {len(data_dict) - len(data_dict_filtered)} parents with too low intensity.")

# %% COULDDO[WIP]: Score nuclides by multiplying their emission spectrum with the measured spectrum
# for parent, (lit_energies, intensities) in data_dict.items():
#     E_parent = np.zeros_like(E_U)
#     for lit_energy, lit_intensity in zip(lit_energies, intensities):
#         E_parent += fit_fn_peak(x=E_U, a=lit_intensity, x_0=lit_energy, sigma=5, N_0=ureg('0 keV'))
#     score = E_parent * N_U
#     print(f"{parent}: {score.sum():.2f}")

# %% Peaks
peaks, _ = find_peaks(
    np.log10(N_U + 1),
    prominence=0.9,
    distance=100,
    # threshold=0.03,
)
peak_fits = [fit_peak(peak, x_U, N_U, plot=False) for peak in peaks]
peak_channels = [fit['x_0'] for fit in peak_fits] * ureg.dimensionless
peak_energies = channel_to_E(peak_channels)

print(f"Found {len(peaks)} peaks")

# %%


def find_nearest_lit_energy(E, data_dict):
    E_RADIUS = ureg('5 keV')  # TODO: test this
    possible_lit_energies = []  # (parent, lit_energy, intensity)
    for parent, (lit_energies, intensities) in data_dict.items():
        for lit_energy, intensity in zip(lit_energies, intensities):
            if abs(lit_energy - E) < E_RADIUS:
                possible_lit_energies.append((parent, lit_energy, intensity))

    # if len(possible_lit_energies) == 0:
        # return None
    assert len(possible_lit_energies) > 0, f"No lit energies found for {E}"

    possible_lit_energies.sort(key=lambda x: x[2], reverse=True)
    print(
        f"Found {len(possible_lit_energies)} possible lit energies with intensities {[f'{x[2].m:.1f}' for x in possible_lit_energies]}")
    return possible_lit_energies[0]

# for peak_energy in peak_energies:
#     print(f"█ Peak at {peak_energy:.2f}")
#     nearest_lit_energy = find_nearest_lit_energy(peak_energy, data_dict_filtered)
#     if nearest_lit_energy is not None:
#         parent, lit_energy, intensity = nearest_lit_energy
#         print(f"Nearest lit energy: {parent} at {lit_energy:.2f} with intensity {intensity:.2f}")


nearest_lit_energies = [find_nearest_lit_energy(peak_energy, data_dict) for peak_energy in peak_energies]
# COULDDO: again, filter out energies < 100 keV

# Tabelle generieren
generate_table.generate_table_pint(
    "build/tab/5_uranstein.tex",
    (r"E", ureg.kiloelectron_volt, tools.nominal_values(peak_energies)),
    (r"E_\text{lit}", ureg.kiloelectron_volt, [x[1] for x in nearest_lit_energies]),
    (r"(E - E_\text{lit})", ureg.kiloelectron_volt, [x[1] - tools.nominal_value(peak_energy)
     for x, peak_energy in zip(nearest_lit_energies, peak_energies)]),
    (r"W", ureg.percent, [x[2] for x in nearest_lit_energies]),
    (r"\text{Nuklid}", str, [r"\ce{^" + ''.join(x[0].split('-')[::-1]) + r"}" for x in nearest_lit_energies]),
)


# %% Plotte Spektrum
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'dimensionless', r"\text{Kanal}", "N") as plt2:
# TODO: \text funzt nicht ohne BUILD…
with tools.plot_context(plt, 'keV', 'dimensionless', "E", "N") as plt2:
    plt2.plot(E_U, N_U, label="Messwerte")
    # plt2.plot(E_U, np.convolve(N_U.m, np.ones(20)/20, mode='same'), fmt='-', label="Messwerte (geglättet)")  # TODO: Gut so?
    plt2.plot(E_U[peaks], N_U[peaks], 'x', show_xerr=False, label="Peaks")

    yheight = 1/len(data_dict_filtered)

    for i, (parent, (lit_energies, intensities)) in enumerate(list(data_dict_filtered.items())[::-1]):
        plt_text = plt.text(
            0.05,
            (i + 1/2)*yheight,
            f"{parent}",
            horizontalalignment='left',
            verticalalignment='center',
            transform=plt.gca().transAxes
        )
        plt_text.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))

        for lit_energy, lit_intensity in zip(lit_energies, intensities):
            if lit_energy.m > max(E_U.m):
                continue
            # plt.axvline(lit_energy.to('keV').m, color=f'C{i+2}', alpha=intensity_to_alpha(lit_intensity, exponent=0.4))
            # line that only occupies 1/len(data_dict_filtered) of the vertical space
            plt.axvline(
                lit_energy.to('keV').m,
                color=f'C{i+2}',
                alpha=intensity_to_alpha(lit_intensity, exponent=0.4),
                ymin=i*yheight,
                ymax=(i+1)*yheight,
                zorder=5,
            )

# add labels for axvlines
handles, labels = plt.gca().get_legend_handles_labels()
# handles += [
#     mpatches.Patch(color=f'C{i+2}', label=f"{parent} (Literaturwerte)")
#     for i, parent in enumerate(list(data_dict_filtered.keys()))
# ]

# plt.xlim(right=750)
plt.ylim(bottom=0.5)
plt.yscale('log')
plt.legend(handles=handles, loc='lower right')
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_uranstein.pdf")
if tools.PLOTS:
    plt.show()

# %%
