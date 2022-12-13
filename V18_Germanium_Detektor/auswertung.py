# %%
import datetime
from io import StringIO
from pathlib import Path

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
    """
    Loads a LARAWEB file and returns the energy and intensity columns, limited to γ decays.
    Data source: http://www.lnhb.fr/Laraweb/index.php
    """
    assert filename.endswith('.lara.txt')
    # LARA files have a header (with varying length) and footer that we need to skip.
    # That is all this block does.
    contents = Path(filename).read_text()
    lines = contents.splitlines()
    # find lines that only contain "-" or "="
    delimiter_mask = [all(c in "-=" for c in line) for line in lines]
    # find the two delimiter lines' indices
    delimiter_indices = np.where(delimiter_mask)[0]
    assert len(delimiter_indices) == 2
    # get the lines between the delimiters as a file-like object
    data_lines = lines[delimiter_indices[0] + 1:delimiter_indices[1]]
    data = StringIO("\n".join(data_lines))

    df = pd.read_csv(data, sep=" ; ", engine='python', index_col=False)

    # filter out non-gamma decays
    df = df[df['Type'] == 'g']

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
    58, 62, 150, 200, 226, 415, 536, 576, 764,
}
peaks = [peak for peak in peaks if peak not in peaks_to_ignore]

peaks.append(4200)  # TODO!
peaks = list(sorted(peaks))

assert len(peaks) == 11, f"Expected 11 peaks, found {len(peaks)}: {peaks}"  # for testing (@Mampfzwerg)
peaks

# %%


def fit_fn_peak(x, a, x_0, sigma, N_0):
    return a * np.exp(-((x - x_0) ** 2) / (2 * sigma ** 2)) + N_0


def fit_peak(peak, x, N, plot=True, channel_to_E=None):
    # TODO: Stelle sicher, dass x Kanäle, nicht Energien sind.
    # assert not isinstance(x, ureg.Quantity) or x.units == ureg.dimensionless

    # print(f"Peak bei ({peak}, {N[peak].m})")
    FIT_RADIUS = 40  # Radius des Bereichs, in dem die Gauß-Kurve gefittet wird
    range_x = x[(peak - FIT_RADIUS): (peak + FIT_RADIUS)]
    range_N = N[(peak - FIT_RADIUS): (peak + FIT_RADIUS)]

    # ▒ Fitte Gauß-Kurve
    a, x_0, σ, N_0 = tools.curve_fit(
        fit_fn_peak,
        range_x, range_N.m,
        p0=[max(range_N), peak, 1, min(range_N)],
    )

    fwhm = 2 * σ * np.sqrt(2 * np.log(2))
    fwhm_height = a / 2 + N_0
    fwtm = 2 * σ * np.sqrt(2 * np.log(10))
    fwtm_height = a / 10 + N_0

    if plot:
        with tools.plot_context(plt, 'dimensionless', 'dimensionless', r"xTODO", r"N") as plt2:
            plt2.plot(range_x, range_N, label="Messwerte")
            plt2.plot(
                range_x,
                fit_fn_peak(range_x, *[param.n for param in [a, x_0, σ, N_0]]),
                label="Fit",
            )

            plt2.plot([x_0 - fwhm / 2, x_0 + fwhm / 2], [fwhm_height, fwhm_height], label="FWHM")
            plt2.plot([x_0 - fwtm / 2, x_0 + fwtm / 2], [fwtm_height, fwtm_height], label="FWTM")

        plt.legend()
        plt.tight_layout()
        plt.show()

    return {
        'a': a,
        'x_0': x_0,
        'σ': σ,
        'N_0': N_0,
        # ---
        'fwhm': fwhm,
        'fwtm': fwtm,
    }


peak_fits = [fit_peak(peak, x_calib, N_calib, plot=False) for peak in peaks]
peak_channels = [fit['x_0'] for fit in peak_fits]
peak_channels_n = tools.nominal_values(peak_channels * ureg.dimensionless)


# %%
lit_energies_all, lit_intensities_all = load_lara("data/emissions/Eu-152.lara.txt")
# select the most intense energies, so that len(lit_energies) == len(peak_channels)
lit_energies, lit_intensities = n_most_intense(lit_energies_all, lit_intensities_all, n=len(peak_channels))

# %% Lineare Regression (Zuordnung Energie ↔ Kanal)

m, n = tools.linregress(peak_channels_n, lit_energies)
print(f"m={m}, n={n}")


def channel_to_E(x):
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
plt.tight_layout()
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

# t_hw = ureg('4943 d')  # @Mampfzwerg
t_hw = ufloat(13.522, 0.016) * ureg.a  # [lara]
A_0 = ufloat(4130, 60) * ureg.Bq  # Aktivität am Tag der Herstellung [Versuchsanleitung]
A = A_0 * unp.exp((-np.log(2) / t_hw * age_probe).to('dimensionless').m)  # Aktivität am Tag der Durchführung
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
plt.tight_layout()
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
x_pint = x * ureg.dimensionless

E = channel_to_E(x)  # NOTE: This means E has uncertainties!

# Peaks
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
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/spektrum_137-Cs.pdf")
plt.show()

# %% Fit: Photopeak
photopeak_fit = fit_peak(photopeak, x, N, plot=True)
E_photopeak_fit = channel_to_E(photopeak_fit['x_0'])

E_photopeak_lit = ureg('661.657 keV')  # TODO: Quelle, Unsicherheit
print(tools.fmt_compare_to_ref(E_photopeak_fit, E_photopeak_lit, name="Photopeak (Fit vs. Literatur)"))

E_fwhm_fit = channel_to_E(photopeak_fit['fwhm'])
E_fwtm_fit = channel_to_E(photopeak_fit['fwtm'])

print(f"FWHM: {E_fwhm_fit:.2f}")
print(f"FWTM: {E_fwtm_fit:.2f}")
print(f"FWHM/FWTM: {(E_fwhm_fit/E_fwtm_fit):.2f}")

# %% Compton-Kante


def calc_compton_edge(E_photopeak):
    ε = E_photopeak / (ureg.m_e * ureg.c**2)
    ε.ito('dimensionless')
    return E_photopeak * 2*ε / (1 + 2*ε)


E_compton_lit = calc_compton_edge(E_photopeak_lit)

# Fit links-rechts
# TODO: subject to confirmation bias
mask_l = (E < E_compton_lit) & (E > (E_compton_lit - 50*ureg.keV))
mask_r = (E > E_compton_lit) & (E < (E_compton_lit + 10*ureg.keV))
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

    plt.axvline(tools.nominal_value(E_compton_fit).m, color='C2', label="Compton-Kante (Fit)")
    plt.axvline(E_compton_lit.m, color='C3', label="Compton-Kante (Literatur)")
plt.yscale('linear')
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/compton-kante.pdf")
# %%


def fit_fn_klein_nishina(E, a, N_0):
    # E: Energie des gestoßenen Elektrons

    # ε = ureg.epsilon_0  # ?
    ε = tools.nominal_value(E_photopeak) / (ureg.m_e * ureg.c**2)  # muss dimensionslos sein; siehe "2 +"
    m_0 = ureg.m_e  # ?
    # E_γ = h*nu
    E_γ = tools.nominal_value(E_photopeak_fit)  # Energie des einfallenden Gammaquants (=hν)

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
mask_compton_plot = (E > ureg('300 keV')) & (E < ureg('500 keV'))

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

# %% Plot: Klein-Nishina
plt.figure()
with tools.plot_context(plt, 'keV', 'dimensionless', r"E_\gamma", r"N") as plt2:
    mask_compton_plotonly = mask_compton_plot & ~mask_compton_fit
    plt2.plot(E[mask_compton_plot], N[mask_compton_plot], 'x', show_xerr=False, color='C0', alpha=0.5, label="Messwerte")
    plt2.plot(E[mask_compton_plotonly], N[mask_compton_plotonly], 'x', show_xerr=False, color='C1', alpha=0.5, label="nur Plot")
    plt2.plot(E[mask_compton_plot], fit_fn_klein_nishina(E[mask_compton_plot], a, N_0), show_xerr=False, color='C2', label="Fit")
plt.yscale('linear')
plt.legend()
plt.tight_layout()
if tools.BUILD:
    plt.savefig("build/plt/klein-nishina.pdf")

# %%
