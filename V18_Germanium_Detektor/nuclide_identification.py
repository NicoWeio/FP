import generate_table
import matplotlib.pyplot as plt
import numpy as np
from common import (fit_peak, get_Ω, load_lara, load_spe, peak_fits_to_arrays,
                    plot_energyspectrum, ureg)
from scipy.signal import find_peaks

import tools


def main(channel_to_E, E_to_Q):
    x, N, T = load_spe("data/2022-11-28/4_U.Spe")
    E = channel_to_E(x)

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
        np.log10(N + 1),
        prominence=0.9,
        distance=100,
        # threshold=0.03,
    )
    peak_fits = [fit_peak(peak, x, N, channel_to_E=channel_to_E, plot=False) for peak in peaks]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)
    peak_energies = channel_to_E(peak_fit_arrays['x_0'])

    print(f"Found {len(peaks)} peaks")

    # %%

    def find_nearest_lit_energy(E, data_dict):
        E_RADIUS = ureg('5 keV')
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

    nearest_lit_energies_results = [find_nearest_lit_energy(peak_energy, data_dict) for peak_energy in peak_energies]
    nearest_parents, nearest_lit_energies, nearest_lit_intensities = zip(*nearest_lit_energies_results)
    nearest_parents = np.array(nearest_parents)
    # COULDDO: again, filter out energies < 100 keV

    # █ Aktivitäten berechnen
    Ω = get_Ω(with_spacer=False)
    Q = E_to_Q(peak_fit_arrays['E'])

    A = 4*np.pi*peak_fit_arrays['Z'] / (Ω * Q * nearest_lit_intensities * T)  # Aktivität
    assert A.check('Bq')

    # for each parent, get the cumulated activity
    parent_to_A = {
        # NOTE: Mittelwert war Quatsch
        # parent: np.mean(A[nearest_parents == parent])
        # parent: tools.ufloat_from_list(tools.nominal_values(A[nearest_parents == parent])).to('Bq')
        parent: np.sum(A[nearest_parents == parent])
        for parent in set(nearest_parents)
    }
    print("Mittlere Aktivitäten:")
    for parent, A_ in parent_to_A.items():
        print(f"{parent}: {A_:.2f}")

    # █ Tabellen generieren

    generate_table.generate_table_pint(
        "build/tab/5_uranstein.tex",
        # ↓ Fit-Parameter
        (r"x_0", ureg.dimensionless, peak_fit_arrays['x_0'], 1),
        (r"a", ureg.dimensionless, peak_fit_arrays['a'], 1),
        (r"\sigma", ureg.dimensionless, peak_fit_arrays['σ'], 2),
        (r"N_0", ureg.dimensionless, peak_fit_arrays['N_0'], 2),
        # ↓ berechnete Werte
        (r"E(x_0)", ureg.kiloelectron_volt, peak_fit_arrays['E']),
        (r"Z(a, \sigma)", ureg.dimensionless, peak_fit_arrays['Z'], 0),
        # ↓ Literaturwerte
        # (r"E_\text{lit}", ureg.kiloelectron_volt, [x[1] for x in nearest_lit_energies_results]),
        # (r"W", ureg.percent, lit_intensities),
    )

    generate_table.generate_table_pint(
        "build/tab/5_uranstein_activities.tex",
        #
        (r"\text{Nuklid}", str, [r"\ce{^" + ''.join(p.split('-')[::-1]) + r"}" for p in nearest_parents]),
        (r"E_\text{lit}", ureg.kiloelectron_volt, nearest_lit_energies),
        (r"W", ureg.percent, nearest_lit_intensities),
        #
        (r"E(x_0)", ureg.kiloelectron_volt, peak_fit_arrays['E']),
        (r"Z(a, \sigma)", ureg.dimensionless, peak_fit_arrays['Z'], 0),
        (r"Q", ureg.percent, tools.nominal_values(Q), 1),
        (r"A", ureg.Bq, tools.nominal_values(A), 0),
    )

    generate_table.generate_table_pint(
        "build/tab/5_uranstein_activities_mean.tex",
        (r"\text{Nuklid}", str, [r"\ce{^" + ''.join(p.split('-')[::-1]) + r"}" for p in parent_to_A.keys()]),
        (r"\bar{A}", ureg.Bq, tools.nominal_values(tools.pintify(list(parent_to_A.values()))), 0),
    )

    # generate_table.generate_table_pint(
    #     "build/tab/5_uranstein_foo.tex",
    #     (r"E", ureg.kiloelectron_volt, tools.nominal_values(peak_energies)),
    #     (r"E_\text{lit}", ureg.kiloelectron_volt, [x[1] for x in nearest_lit_energies]),
    #     (r"(E - E_\text{lit})", ureg.kiloelectron_volt, [x[1] - tools.nominal_value(peak_energy)
    #                                                      for x, peak_energy in zip(nearest_lit_energies, peak_energies)]),
    #     (r"W", ureg.percent, [x[2] for x in nearest_lit_energies]),
    #     (r"\text{Nuklid}", str, [r"\ce{^" + ''.join(x[0].split('-')[::-1]) + r"}" for x in nearest_lit_energies]),
    # )

    # %% Plotte Spektrum
    plot_energyspectrum(
        tools.nominal_values(E), N,
        peak_indices=peaks,
        lit_energies_dict=data_dict_filtered,  # COULDDO: Reihenfolge umkehren
        path="build/plt/spektrum_uranstein.pdf",
        stack_lit_energies=True,
    )
