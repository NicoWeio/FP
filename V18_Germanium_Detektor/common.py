from io import StringIO
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pint
import uncertainties.unumpy as unp
from rich.console import Console
from uncertainties import UFloat, ufloat

import tools

console = Console()

ureg = pint.UnitRegistry()
ureg.define('percent = 1 / 100 = %')

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


def plot_energyspectrum(
    E, N,
    lit_energies_dict=None,
    stack_lit_energies=False,
    path=None,
    smooth_over: int | None = None,
):
    """
    lit_energies_dict: dict of {parent: (energies, intensities)}
    path: if given (and tools.BUILD == True), save the plot to this path
    smooth_over: if given, add a smoothed version over this many points to the plot
    """
    yheight = 1/len(lit_energies_dict)
    E_bounds = tools.bounds(E)

    plt.figure()
    with tools.plot_context(plt, 'keV', 'dimensionless', "E", "N") as plt2:
        plt2.plot(E, N, label="Messwerte")

        if smooth_over:
            N_smooth = np.convolve(N.m, np.ones(smooth_over)/smooth_over, mode='same') * ureg.dimensionless
            plt2.plot(E, N_smooth, fmt='-', label="Messwerte (geglättet)")  # TODO: Gut so?

        # plt.plot(peaks, N[peaks], 'x', label="Peaks")  # TODO: Peaks tatsächlich suchen…

        for i, (parent, (lit_energies, intensities)) in enumerate(list(lit_energies_dict.items())):
            for lit_energy, lit_intensity in zip(lit_energies, intensities):
                if lit_energy < E_bounds[0] or lit_energy > E_bounds[1]:
                    continue

                plt.axvline(
                    lit_energy.to('keV').m,
                    color=f'C{i+2}',
                    alpha=intensity_to_alpha(lit_intensity),
                    zorder=5,
                    **({
                        'ymin': i*yheight,
                        'ymax': (i+1)*yheight,
                    }
                        if stack_lit_energies else {}
                    )
                )

        # for i, (parent, (lit_energies, intensities)) in enumerate(list(data_dict_filtered.items())[::-1]):
            if stack_lit_energies:
                plt.text(
                    0.05,
                    (i + 1/2)*yheight,
                    f"{parent}",
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=plt.gca().transAxes
                ).set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))

    # add labels for axvlines
    handles, labels = plt.gca().get_legend_handles_labels()
    if not stack_lit_energies:
        handles += [
            mpatches.Patch(color=f'C{i+2}', label=f"{parent} (Literaturwerte)")
            for i, parent in enumerate(lit_energies_dict.keys())
        ]

    # plt.xlim(
    #     left=0,  # COULDDO: Why is this necessary!?
    #     right=1.1 * max([lit_energies.max() for (lit_energies, intensities) in lit_energies_dict.values()]).to('keV').m,
    # )
    # plt.xlim(*E_bounds.to('keV').m)  # COULDDO
    plt.ylim(bottom=0.5)
    plt.yscale('log')
    plt.legend(handles=handles)
    plt.tight_layout()
    if tools.BUILD and path:
        plt.savefig(path)
    if tools.PLOTS:
        plt.show()
    # plt.close()  # COULDDO

    return plt


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
