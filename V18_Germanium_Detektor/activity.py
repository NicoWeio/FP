import generate_table
import numpy as np
from common import (fit_peak, load_lara, load_spe, peak_fits_to_arrays,
                    plot_energyspectrum, ureg)
from scipy.signal import find_peaks

import tools


def main(channel_to_E):
    x, N, T = load_spe("data/2022-11-28/3_Ba.Spe")
    E = channel_to_E(x)

    lit_energies_all_Ba, intensities_all_Ba = load_lara("data/emissions/Ba-133.lara.txt")
    lit_energies_all_Sb, intensities_all_Sb = load_lara("data/emissions/Sb-125.lara.txt")

    # Peaks
    N_smooth = np.convolve(N.m, np.ones(50)/50, mode='same') * ureg.dimensionless

    raw_peaks, _ = find_peaks(
        np.log(N_smooth + 1),
        prominence=0.6,
        distance=100,
    )

    peak_fits = [fit_peak(peak, x, N, channel_to_E=channel_to_E, plot=False) for peak in raw_peaks]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)

    # █ Tabelle generieren
    # COULDDO: Unsicherheiten
    generate_table.generate_table_pint(
        "build/tab/4_Ba-133.tex",
        # ↓ Fit-Parameter
        (r"x_0", ureg.dimensionless, peak_fit_arrays['x_0'], 1),
        (r"a", ureg.dimensionless, peak_fit_arrays['a'], 1),
        (r"\sigma", ureg.dimensionless, peak_fit_arrays['σ'], 2),
        (r"N_0", ureg.dimensionless, peak_fit_arrays['N_0'], 2),
        # ↓ berechnete Werte
        (r"E(x_0)", ureg.kiloelectron_volt, peak_fit_arrays['E']),
        (r"Z(a, \sigma)", ureg.dimensionless, peak_fit_arrays['Z'], 0),
    )

    # ▒ Plotte Spektrum
    plot_energyspectrum(
        tools.nominal_values(E), N,
        lit_energies_dict={
            "Ba-133": (lit_energies_all_Ba, intensities_all_Ba),
            "Sb-125": (lit_energies_all_Sb, intensities_all_Sb),
        },
        path="build/plt/spektrum_133-Ba.pdf",
        # xlim='lit_energies', # Sb hat einige schwache, hochenergetische Emissionen, daher funktioniert das hier nicht.
        xlim=(0, 750),
    )
