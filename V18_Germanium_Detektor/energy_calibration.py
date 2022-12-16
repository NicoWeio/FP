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

    # ⚠ Wir suchen die Peaks in N_smooth. Daher später `peak_fit_arrays['true_peak']` verwenden!
    # COULDDO: mehr Peaks?
    raw_peaks, _ = find_peaks(
        np.log(N_smooth + 1),
        prominence=0.6,
        distance=100,
    )

    peaks_to_ignore = [
        79,
        5861,
        6287,
    ]
    raw_peaks = [p for p in raw_peaks if p not in peaks_to_ignore]
    raw_peaks = list(sorted(raw_peaks))

    assert len(raw_peaks) == 11, f"Expected 11 peaks, found {len(raw_peaks)}: {raw_peaks}"  # for testing (@Mampfzwerg)

    peak_fits = [fit_peak(peak, x, N, plot=False) for peak in raw_peaks]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)

    lit_energies_all, lit_intensities_all = load_lara("data/emissions/Eu-152.lara.txt")
    # select the most intense energies, so that len(lit_energies) == len(peaks)
    lit_energies, lit_intensities = n_most_intense(lit_energies_all, lit_intensities_all, n=len(raw_peaks))

    # Lineare Regression (Zuordnung Energie ↔ Kanal)
    m, n = tools.linregress(peak_fit_arrays['x_0'], lit_energies)
    print(f"m={m}, n={n}")

    def channel_to_E(x):
        # COULDDO: warn about negative energies
        return m * x + n

    def E_to_channel(E):
        return (E - n) / m

    # ▒ Plot: Energie(Kanal)
    plt.figure()
    with tools.plot_context(plt, 'dimensionless', 'keV', "x", "E") as plt2:
        plt2.plot(peak_fit_arrays['x_0'], channel_to_E(peak_fit_arrays['x_0']), show_yerr=False, label="Regressionsgerade")
        plt2.plot(peak_fit_arrays['x_0'], lit_energies, 'x', zorder=5, label="Literaturwerte")
    plt.legend()
    plt.tight_layout()
    if tools.BUILD:
        plt.savefig("build/plt/energy_calibration.pdf")
    # if tools.PLOTS:
    #     plt.show()

    # ▒ Plot: Spektrum
    plot_energyspectrum(
        tools.nominal_values(channel_to_E(x)), N,
        peak_indices=peak_fit_arrays['true_peak'],
        lit_energies_dict={
            "Eu-152": (lit_energies_all, lit_intensities_all),
        },
        path="build/plt/spektrum_152-Eu.pdf",
        # smooth_over=20,
        # stack_lit_energies=True,
    )

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
