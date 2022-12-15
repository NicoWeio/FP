import matplotlib.patches as mpatches
import generate_table
import matplotlib.pyplot as plt
import numpy as np
from common import (fit_peak, load_lara, load_spe, n_most_intense,
                    peak_fits_to_arrays, ureg, plot_energyspectrum, intensity_to_alpha)
from scipy.signal import find_peaks

import tools


def main():
    x, N, T = load_spe("data/2022-11-28/1_Eu.Spe")

    N_smooth = np.convolve(N.m, np.ones(50)/50, mode='same') * ureg.dimensionless

    # finde Peaks
    peaks, _ = find_peaks(
        # np.log(N + 1),
        # prominence=3.2,
        # distance=100,

        np.log(N_smooth + 1),
        # np.log(N_for_peaks + 1),
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

    peak_fits = [fit_peak(peak, x, N, plot=False) for peak in peaks]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)

    lit_energies_all, lit_intensities_all = load_lara("data/emissions/Eu-152.lara.txt")
    # select the most intense energies, so that len(lit_energies) == len(peaks)
    lit_energies, lit_intensities = n_most_intense(lit_energies_all, lit_intensities_all, n=len(peaks))

    # Lineare Regression (Zuordnung Energie ↔ Kanal)
    m, n = tools.linregress(peak_fit_arrays['x_0'], lit_energies)
    print(f"m={m}, n={n}")

    def channel_to_E(x):
        # COULDDO: warn about negative energies
        return m * x + n

    def E_to_channel(E):
        return (E - n) / m

    lit_channels = tools.nominal_values(E_to_channel(lit_energies))
    lit_channels_all = tools.nominal_values(E_to_channel(lit_energies_all))

    # ▒ Plot: Energie(Kanal)
    plt.figure()
    with tools.plot_context(plt, 'dimensionless', 'keV', "x", "E") as plt2:
        plt2.plot(peak_fit_arrays['x_0'], channel_to_E(peak_fit_arrays['x_0']), show_yerr=False, label="Regressionsgerade")
        plt2.plot(peak_fit_arrays['x_0'], lit_energies, 'x', zorder=5, label="Literaturwerte")
    plt.legend()
    plt.tight_layout()
    if tools.BUILD:
        plt.savefig("build/plt/energy_calibration.pdf")
    if tools.PLOTS:
        plt.show()

    # ▒ Plot: Spektrum
    plot_energyspectrum(
        tools.nominal_values(channel_to_E(x)), N,
        lit_energies_dict={
            "Eu-152": (lit_energies_all, lit_intensities_all),
        },
        path="build/plt/spektrum_152-Eu.pdf",
        # smooth_over=20,
        # stack_lit_energies=True,
    )

    # ▒ Plot: Spektrum → deprecated
    plt.figure()
    with tools.plot_context(plt, 'dimensionless', 'dimensionless', "x", "N") as plt2:
        plt.plot(N, label="Messwerte")
        # plt.plot(x, N_smooth, '-', label="Messwerte (geglättet)")

        plt.plot(peak_fit_arrays['true_peak'], N[peak_fit_arrays['true_peak']], 'x', label="Peaks")

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

    # █ Tabelle generieren
    # COULDDO: Unsicherheiten
    generate_table.generate_table_pint(
        "build/tab/1_energiekalibrierung.tex",
        (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
        (r"x_0", ureg.dimensionless, peak_fit_arrays['x_0'], 1),
        (r"a", ureg.dimensionless, peak_fit_arrays['a'], 1),
        (r"\sigma", ureg.dimensionless, peak_fit_arrays['σ'], 2),
        (r"N_0", ureg.dimensionless, peak_fit_arrays['N_0'], 2),
    )

    return {
        'lit_energies': lit_energies,
        'lit_intensities': lit_intensities,
        'Z': peak_fit_arrays['Z'],
        'T': T,
        'channel_to_E': channel_to_E,
        'E_to_channel': E_to_channel,
    }
