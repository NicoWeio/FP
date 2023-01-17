import generate_table
import matplotlib.pyplot as plt
import numpy as np
from common import (calc_ε, fit_peak, fit_peak_linregress, load_spe,
                    peak_fits_to_arrays, plot_energyspectrum, ureg, console)
from scipy.signal import find_peaks

import tools


@ureg.check('[energy]')
def calc_compton_edge(E_photopeak):
    ε = calc_ε(E_photopeak)
    return E_photopeak * 2*ε / (1 + 2*ε)


def get_fit_fn_klein_nishina(E_γ):
    """
    E_γ: Energie des einfallenden Gammaquants (=hν)
    """
    def fit_fn_klein_nishina(E, a, N_0):
        """
        E: Energie des gestoßenen Elektrons
        """
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

    return fit_fn_klein_nishina


def main_comptonkante(E, N, photopeak_fit, E_photopeak_lit):
    E_compton_lit = calc_compton_edge(E_photopeak_lit)

    fit_center = ureg('470 keV')
    mask_l = (E < fit_center) & (E > (fit_center - ureg('100 keV')))
    mask_r = (E > fit_center) & (E < (fit_center + ureg('20 keV')))
    mask_lr = mask_l | mask_r

    # v1 ↓
    # N_smooth = np.convolve(N.m, np.ones(10)/10, mode='same') * ureg.dimensionless
    # peak_i = np.argmax(N_smooth[mask_lr])

    # v2 ↓
    # smoothe N by chunking it into bins of 16 and taking the 95% percentile
    N_smooth = np.percentile(N.m.reshape(-1, 16), 95, axis=1)
    # …but keep the original length
    N_smooth = np.repeat(N_smooth, 16) * ureg.dimensionless
    # peak_i = np.argmax(N_smooth[mask_lr])
    peak_i = np.argmin(np.diff(N_smooth, append=0)[mask_lr])

    E_compton_fit = E[mask_lr][peak_i]

    print(tools.fmt_compare_to_ref(E_compton_fit, E_compton_lit, "Compton-Kante", unit='keV'))

    # █ Plot: Compton-Kante
    plt.figure()
    with tools.plot_context(plt, 'keV', 'dimensionless', r"E", r"N") as plt2:
        plt2.plot(E[mask_lr], N[mask_lr], 'x', show_xerr=False, color='C0', alpha=0.5, label="Messwerte")
        # plt2.plot(E[mask_lr], N_smooth[mask_lr], '-', show_xerr=False, color='C1', label="Messwerte (geglättet)")
        # plt2.plot(E[peaks], N[peaks], fmt='x', label="Peaks")

        plt.axvline(tools.nominal_value(E_compton_fit).to('keV').m, color='C2', label="Compton-Kante")
        plt.axvline(E_compton_lit.to('keV').m, color='C3', label="Compton-Kante (Literatur)")
    plt.ylim(top=max(N[mask_lr]))
    plt.yscale('linear')
    plt.legend(loc='lower left')
    plt.tight_layout()
    if tools.BUILD:
        plt.savefig("build/plt/compton-kante.pdf")


    # █ Fit: Klein-Nishina
    mask_compton_fit = (E > ureg('350 keV')) & (E < E_compton_fit)
    mask_compton_plot = (E > ureg('200 keV')) & (E < ureg('500 keV'))

    fit_fn_klein_nishina = get_fit_fn_klein_nishina(E_γ=tools.nominal_value(photopeak_fit['E']))

    a, N_0 = tools.pint_curve_fit(
        fit_fn_klein_nishina,
        tools.nominal_values(E[mask_compton_fit]), N[mask_compton_fit],
        (ureg.dimensionless, ureg.dimensionless),
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

    # █ Plot: Klein-Nishina
    plt.figure()
    with tools.plot_context(plt, 'keV', 'dimensionless', "E", "N") as plt2:
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


def main(channel_to_E):
    x, N, T = load_spe("data/2022-11-28/2_Cs.Spe")
    E = channel_to_E(x)  # NOTE: This means E has uncertainties!

    # Peaks
    peaks, _ = find_peaks(N, height=50, distance=1000)

    assert len(peaks) == 2, f"{len(peaks)} statt 2 Peaks (Rückstreupeak und Photopeak) gefunden."
    rueckstreupeak_raw, photopeak_raw = peaks

    # Fit: Rückstreupeak & Photopeak
    # peak_fits = [fit_peak(peak, x, N, plot=True, channel_to_E=channel_to_E) for peak in peaks]
    peak_fits = [
        fit_peak(
            peak, x, N, plot=True,
            plot_path=f"build/plt/3_gauss_{['rueckstreu','photo'][i]}peak.pdf",
            channel_to_E=channel_to_E,
            fit_radius=([150, 40][i]),
            # fit_radius=([150, 20][i]),
        )
        for i, peak in enumerate(peaks)
    ]
    peak_fit_arrays = peak_fits_to_arrays(peak_fits)
    rueckstreupeak_fit, photopeak_fit = peak_fits

    E_photopeak_lit = ureg('661.657 keV')  # Quelle: LARA
    print(tools.fmt_compare_to_ref(photopeak_fit['E'], E_photopeak_lit, name="Photopeak (Fit vs. Literatur)"))

    # NOTE: E = E_right - E_left = (m·channel_right+b) - (m·channel_left+b) = m·(channel_right-channel_left)
    E_fwhm_fit = channel_to_E.m * photopeak_fit['fwhm']
    E_fwtm_fit = channel_to_E.m * photopeak_fit['fwtm']

    # NOTE: Nicht verwirren lassen. Das sind die aus dem Gauß-Fit berechneten Werte.
    print(f"Photopeak (via Gauß-Fit):")
    print(f"\tFWHM: {E_fwhm_fit:.2f}")
    print(f"\tFWTM: {E_fwtm_fit:.2f}")

    # █ FWHM und FWTM des Photopeaks (via Gauß-Fits)
    fwhm_height = tools.nominal_value((photopeak_fit['a'] / 2 + photopeak_fit['N_0']) * ureg.dimensionless)
    fwtm_height = tools.nominal_value((photopeak_fit['a'] / 10 + photopeak_fit['N_0']) * ureg.dimensionless)

    fwtm_div_fwhm = E_fwtm_fit / E_fwhm_fit
    assert fwtm_div_fwhm > 1
    fwtm_div_fwhm_optimal = 1.823  # TODO: Quelle!

    print(tools.fmt_compare_to_ref(fwtm_div_fwhm, fwtm_div_fwhm_optimal,
          name="Photopeak: FWTM/FWHM: gemessen vs. optimal gaußförmig"))

    # Plot: N(E) – Spektrum von 137Cs
    plot_energyspectrum(
        tools.nominal_values(channel_to_E(x)), N,
        path="build/plt/spektrum_137-Cs.pdf",
        plot_cb=lambda _, plt2: (
            plt2.plot(E[rueckstreupeak_raw], N[rueckstreupeak_raw], 'x', show_xerr=False, label="Rückstreupeak")
            and
            plt2.plot(E[photopeak_raw], N[photopeak_raw], 'x', show_xerr=False, label="Photopeak")
        ),
        xlim=(0, 800),
    )

    # █ Tabelle generieren
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

    # Compton-Kante
    console.rule("→ Compton-Kante", style="gray")
    main_comptonkante(
        E, N,
        photopeak_fit,
        E_photopeak_lit,
    )
