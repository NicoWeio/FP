import generate_table
import numpy as np
from common import (fit_peak, get_Ω, load_lara, load_spe,
                    n_most_intense, peak_fits_to_arrays, plot_energyspectrum,
                    ureg)
from scipy.signal import find_peaks

import tools


def main(channel_to_E, E_to_Q):
    x, N, T = load_spe("data/2022-11-28/3_Ba.Spe")
    E = channel_to_E(x)

    lit_energies_all, lit_intensities_all = load_lara("data/emissions/Ba-133.lara.txt")

    # Peaks
    N_smooth = np.convolve(N.m, np.ones(50)/50, mode='same') * ureg.dimensionless

    raw_peaks, _ = find_peaks(
        np.log(N_smooth + 1),
        prominence=1.2,
        distance=100,
    )
    assert len(raw_peaks) == 5, f"Expected 5 peaks, found {len(raw_peaks)}: {raw_peaks}"
    # raw_peaks = list(sorted(raw_peaks))

    peak_fits = [fit_peak(peak, x, N, channel_to_E=channel_to_E, plot=False) for peak in raw_peaks]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)

    # select the most intense energies, so that len(lit_energies) == len(peaks)
    lit_energies, lit_intensities = n_most_intense(lit_energies_all, lit_intensities_all, n=len(raw_peaks))

    # █ Aktivitäten berechnen
    Ω = get_Ω()
    Q = E_to_Q(peak_fit_arrays['E'])

    A = 4*np.pi*peak_fit_arrays['Z'] / (Ω * Q * lit_intensities * T)  # Aktivität
    assert A.check('Bq')

    # NOTE: Der erste Wert ist deutlich zu klein (auch in den Altprotokollen) und wird daher nicht berücksichtigt
    A_mean = tools.ufloat_from_list(tools.nominal_values(A[1:])).to('Bq')
    print(f"A_mean = {A_mean:.2f}")

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
        # ↓ Literaturwerte
        (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
        # (r"W", ureg.percent, lit_intensities),
    )

    generate_table.generate_table_pint(
        "build/tab/4_Ba-133_activities.tex",
        # ↓ Literaturwerte
        (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
        (r"W", ureg.percent, lit_intensities),
        # ↓ berechnete Werte
        (r"E(x_0)", ureg.kiloelectron_volt, peak_fit_arrays['E']),
        # (r"Z(a, \sigma)", ureg.dimensionless, peak_fit_arrays['Z'], 0),
        (r"Q", ureg.percent, tools.nominal_values(Q), 1),
        (r"A", ureg.Bq, tools.nominal_values(A), 0),
    )

    # ▒ Plotte Spektrum
    plot_energyspectrum(
        tools.nominal_values(E), N,
        peak_indices=peak_fit_arrays['true_peak'],
        lit_energies_dict={
            "Ba-133": (lit_energies_all, lit_intensities_all),
            # "Sb-125": load_lara("data/emissions/Sb-125.lara.txt"),
        },
        path="build/plt/spektrum_133-Ba.pdf",
        # xlim='lit_energies', # Sb hat einige schwache, hochenergetische Emissionen, daher funktioniert das hier nicht.
        xlim=(0, 750),
    )
