import generate_table
import matplotlib.pyplot as plt
import numpy as np
from common import (fit_peak, load_lara, load_spe, peak_fits_to_arrays,
                    plot_energyspectrum, ureg)
from scipy.signal import find_peaks

import tools


def main(channel_to_E):
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
    peak_fits = [fit_peak(peak, x, N, plot=False) for peak in peaks]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)
    peak_energies = channel_to_E(peak_fit_arrays['x_0'])

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
    plot_energyspectrum(
        E, N,
        lit_energies_dict=data_dict_filtered,  # COULDDO: Reihenfolge umkehren
        path="build/plt/spektrum_uranstein.pdf",
        # smooth_over=20,
        stack_lit_energies=True,
    )

    # %% deprecated

    plt.figure()
    with tools.plot_context(plt, 'keV', 'dimensionless', "E", "N") as plt2:
        plt2.plot(E, N, label="Messwerte")
        plt2.plot(E[peaks], N[peaks], 'x', show_xerr=False, label="Peaks")

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
