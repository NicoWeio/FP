"""
NOTE: Die Benennung der Schichtdicke ist uneinheitlich. Manchmal wird sie als z bezeichnet, manchmal als d.
"""

import matplotlib.pyplot as plt
import numpy as np

import tools


def α_to_q(α, λ):
    """
    q ist der Wellenvektorübertrag.
    """
    # TODO: Blind übernommen aus @Mampfzwerg
    # Die Faktoren sehen so aus, als wären sie nur für deg→rad 🤔
    q = 4 * np.pi / λ * np.sin(np.pi / 180 * α)
    return q


def calc_G(α, D, d_Strahl, α_g):
    """
    Berechnet den Geometriefaktor G ≤ 1.
    Dieser wird relevant, wenn die Probe nicht die gesamte Strahlbreite ausfüllt (für α < α_g).

    D: Probendurchmesser
    d_Strahl: Strahlbreite
    α_g: Geometriewinkel
    """
    # Quelle: Versuchsanleitung
    G = D * np.sin(α) / d_Strahl
    G[α > α_g] = 1
    assert all(G <= 1)
    return G


def calc_α_c(δ, ureg):
    # α_c_PS = λ * np.sqrt(litdata['PS']['r_e·ρ'] / np.pi)

    # https://github.com/NicoJG/Fortgeschrittenenpraktikum/blob/677f5868153db0d111c41329f5c517432f6487c9/V44_Reflektometrie/python/messung.py#L157-L158
    return np.sqrt(2 * δ) * ureg.rad  # !?


def calc_parratt(
    α,
    z,
    k,
    δ1, δ2,
    β1, β2,
    σ1, σ2,
    ureg,
    rauigkeit=False,
):
    """
    δ_i: Brechungsindex-Korrektur (n = 1 - δ_i + i·β_i)
    β_i: Brechungsindex-Korrektur (n = 1 - δ_i + i·β_i)
    σ_i: Rauigkeit der Grenzfläche
    """

    # https://de.wikipedia.org/wiki/Brechungsindex#Komplexer_Brechungsindex
    n1 = 1  # Luft
    n2 = 1 - δ1 + 1j*β1  # Polysterol
    n3 = 1 - δ2 + 1j*β2  # Silizium

    # ----------------------------

    kz1 = k * np.sqrt(n1**2 - np.cos(α)**2)  # removed abs(…)
    kz2 = k * np.sqrt(n2**2 - np.cos(α)**2)  # removed abs(…)
    kz3 = k * np.sqrt(n3**2 - np.cos(α)**2)  # removed abs(…)
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2

    if rauigkeit:
        r12 *= np.exp(-2 * kz1 * kz2 * σ1**2)
        r23 *= np.exp(-2 * kz2 * kz3 * σ2**2)
        r13 *= 0  # NOTE: Hierzu hat @Mampfzwerg keine Formel
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    R_parratt = np.abs(x1)**2

    # assert par.check('dimensionless') # does not work for whatever reason
    assert str(R_parratt.units) == 'dimensionless'
    return R_parratt, np.abs(r13)


def do_fit(
    *,  # force keyword arguments
    the_α,
    the_R,
    fit_mask,
    parratt_params,
    k,
    ureg,
):
    """
    Nimmt `parratt_params` als Startwerte entgegen.
    Falls nicht alle Parameter gefittet werden, enthält die Rückgabe dann auch diese.
    """
    def parrat_fitfn(α, *override_parratt_params_tuple):
        # if isinstance(σ1, ureg.Quantity):
        #     # assume all are quantities
        #     σ1 = σ1.to('m').m
        #     σ2 = σ2.to('m').m
        PASSED_PARAMS = ['δ1', 'δ2', 'σ1', 'σ2', 'z']
        UNITS = [ureg.dimensionless, ureg.dimensionless, ureg.m, ureg.m, ureg['Å']]
        override_parratt_params = dict(zip(PASSED_PARAMS, override_parratt_params_tuple, strict=True))
        for key, val in override_parratt_params.items():
            if not isinstance(val, ureg.Quantity):
                override_parratt_params[key] *= UNITS[PASSED_PARAMS.index(key)]

        return calc_parratt(
            α,
            **(parratt_params | override_parratt_params),
            # ↓ the rest
            k=k,
            ureg=ureg,
            rauigkeit=True,
        )[0]

    # def parrat_fitfn_nodim(α, δ1):
    #     return parrat_fitfn(α * ureg.rad, δ1 * ureg.dimensionless)

    # BOUND_POSITIVE = (ureg('0 dimensionless'), None)

    δ1_fit, δ2_fit, σ1_fit, σ2_fit, z_fit = tools.pint_curve_fit(
        parrat_fitfn,
        the_α.to('rad')[fit_mask],
        # tools.nominal_values(R_corr)[fit_mask],
        the_R[fit_mask],
        (ureg.dimensionless, ureg.dimensionless, ureg.m, ureg.m, ureg.angstrom),
        # (ureg.dimensionless, ureg.dimensionless, ureg.m, ureg.m),
        # p0=(1E-6 * ureg.dimensionless, 1E-6 * ureg.dimensionless, ureg('1E-10 m'), ureg('1E-10 m'), ureg('860 Å')),
        p0=tuple(parratt_params[key] for key in ['δ1', 'δ2', 'σ1', 'σ2', 'z']),
        # bounds=(None, None, None, None, None),
        # bounds=(BOUND_POSITIVE, BOUND_POSITIVE, None, None, None),
        # bounds=(BOUND_POSITIVE, BOUND_POSITIVE, (ureg('0 m'), None), (ureg('0 m'), None), None),
        # bounds=(None, None, None, None, (parratt_params['z'] * 0.8, parratt_params['z'] * 1.2)),
        # bounds=(None, None, None, None, (ureg('400 Å'), ureg('600 Å'))),
        # p0=(litdata['PS']['δ'], litdata['Si']['δ']),
        # p0=(litdata['PS']['δ'], litdata['Si']['δ'], 20E-10 * ureg.m, 7E-10 * ureg.m),
        maxfev=5000,
    )

    parratt_params_fit = parratt_params | {
        'δ1': tools.nominal_value(δ1_fit),
        'δ2': tools.nominal_value(δ2_fit),
        'σ1': tools.nominal_value(σ1_fit),
        'σ2': tools.nominal_value(σ2_fit),
        'z': tools.nominal_value(z_fit),
    }

    return parratt_params_fit


def main(
    name,
    mess_refl,
    mess_diff,
    ureg,
    d_Strahl,
    D,
    α_g,
    I_max,
    litdata,
    parratt_params_input,
    plot_configs,
):
    """
    d_Strahl: Strahlbreite (siehe Z-Scan)
    α_g: Geometriewinkel (siehe Rockingscan)
    I_max: Maximale Intensität aus dem Detektorscan
    """
    # █ Messwerte vorbereiten und korrigieren
    assert np.all(mess_refl[0] == mess_diff[0]), "α-Werte stimmen nicht überein"
    α, I_refl = mess_refl
    α, I_diff = mess_diff

    # ↓ Korrektur um diffusen Anteil
    I_corr_diff = I_refl - I_diff
    # ↓ Korrektur um Geometriefaktor
    G = calc_G(α, D=D, d_Strahl=d_Strahl, α_g=α_g)
    # G[0] = np.nan  # NOTE: Workaround for division by zero
    G[0] = G[1]  # NOTE: Workaround for division by zero; variant for fitting…
    I_corr_G = I_refl / G
    # ↓ Korrektur um beides
    I_corr = I_corr_diff / G

    R_corr_diff = I_corr_diff / I_max
    R_corr = I_corr / I_max

    λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
    k = 2*np.pi / λ  # Wellenvektor


    # █ Schichtdicke aus Peaks bestimmen
    # peaks, peak_props = sp.signal.find_peaks(tools.nominal_values(I_corr).to('1/s').m, height=(1E2, 1E4), prominence=5, distance=8)
    # COULDDO: Fast funktioniert es automatisch. Fast…
    peaks = [70, 81, 91, 100, 111, 119, (131), 144, 155]
    assert len(peaks) > 0, "Keine Peaks gefunden"
    print(f"Peak-Indizes: {peaks}")
    # COULDDO: add sem/ufloat
    Δα_mean = np.mean(np.diff(α[peaks].to('rad').m)) * ureg.rad
    d_from_peaks = λ / (2 * Δα_mean)
    print(f"Δα_mean = {Δα_mean}")
    print(f"d_from_peaks = {d_from_peaks.to('nm'):.2f} = {d_from_peaks.to('Å'):.1f}")


    # █ Skalierung der Messdaten (Anpassung an die Parratt-Theoriekurve)
    R_corr_plateau_mean = R_corr[(ureg('0.1°') < α) & (α < ureg('0.2°'))].mean()
    print(f"R_corr_plateau_mean = {R_corr_plateau_mean:.3f}")
    if True:  # Scaling of measured data (in-place)
        print("🛈 Measured data is scaled to match the plateau of the Parratt curve!")
        for var in [I_refl, I_diff, I_corr_diff, I_corr_G, I_corr, R_corr_diff, R_corr,]:
            var /= R_corr_plateau_mean


    # █ Fit der Parratt-Theoriekurve
    # fit_mask = (ureg('0.2°') < α) & (α < ureg('1.5°'))
    # fit_mask = (ureg('0.3°') < α) & (α < ureg('1.25°'))
    fit_mask = (α[peaks[0]] <= α) & (α <= α[peaks[-1]])
    parratt_params_fit = do_fit(
        the_α=α,
        the_R=R_corr,
        fit_mask=fit_mask,
        k=k,
        parratt_params=parratt_params_input,
        ureg=ureg,
    )

    # █ Berechnung der kritischen Winkel
    α_c_PS_fit = calc_α_c(parratt_params_fit['δ1'], ureg=ureg)
    α_c_Si_fit = calc_α_c(parratt_params_fit['δ2'], ureg=ureg)


    # █ Vergleich mit Literaturwerten
    COMPARISON_DATA = [
        {
            'key': 'δ1',
            'name': "δ1",
            # 'ours': parratt_params_fit['δ1'],
            'ref': litdata['PS']['δ'],
        },
        {
            'key': 'δ2',
            'name': "δ2",
            # 'ours': parratt_params_fit['δ2'],
            'ref': litdata['Si']['δ'],
        },
        {
            'key': 'z',
            'name': "z (vs. aus Peak-Abständen)",
            # 'ours': parratt_params_fit['z'],
            'ref': d_from_peaks,
        }
    ]
    for c in COMPARISON_DATA:
        print(tools.fmt_compare_to_ref(parratt_params_fit[c['key']], c['ref'], c['name']))
    # print(f"α_c_PS_fit = {α_c_PS_fit.to('°'):.2f}")
    # print(f"α_c_Si_fit = {α_c_Si_fit.to('°'):.2f}")
    print(tools.fmt_compare_to_ref(α_c_PS_fit, litdata['PS']['α_c'], "α_c_PS_fit", unit='°'))
    print(tools.fmt_compare_to_ref(α_c_Si_fit, litdata['Si']['α_c'], "α_c_Si_fit", unit='°'))


    # █ Berechnung der Parratt-Theoriekurven
    calc_parratt_common_paramters = {
        'α': α.to('rad').m,
        'k': k,
        'ureg': ureg,
        'rauigkeit': True,
    }

    par_input, r13_input = calc_parratt(
        **calc_parratt_common_paramters,
        **parratt_params_input,
    )

    par_fit, r13_fit = calc_parratt(
        **calc_parratt_common_paramters,
        **parratt_params_fit,
    )

    par_fit_glatt, r13_fit_glatt = calc_parratt(
        **(calc_parratt_common_paramters | dict(rauigkeit=False)),
        **parratt_params_fit,
    )


    # █ Plot 1: Messwerte und Korrekturen
    if tools.PLOTS:
        plt.figure()
        # COULDDO: Doppelachse mit Intensität und Reflektivität?
        with tools.plot_context(plt, '°', '1/s', "α", "I") as plt2:
            plt2.plot(α, tools.nominal_values(I_refl), '-', label="Messwerte")
            plt2.plot(α, tools.nominal_values(I_diff), '-', label="Messwerte (diffuse)")
            plt2.plot(α, tools.nominal_values(I_corr_diff), '-', zorder=5, label="korrigiert um diffuse")
            plt2.plot(α, tools.nominal_values(I_corr_G), '-', zorder=5, label="korrigiert um Geometriefaktor")
            plt2.plot(α, tools.nominal_values(I_corr), '-', zorder=5, label="korrigiert")
        plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}_messwerte.pdf")
        plt.show()

    if tools.PLOTS or True:  # TODO
        def plot_schichtdicke(config):
            # α_linspace = tools.linspace(*tools.bounds(α), 1000)

            # COULDDO: Doppelachse mit Intensität und Reflektivität?
            plt.clf()
            with tools.plot_context(plt, '°', 'dimensionless', "α", "R") as plt2:
                # TODO: R_corr_diff passt irgendwie viel besser als R_corr. Eigentlich sollte letzteres benuzt werden…
                if 'R_corr' in config['show']:
                    plt2.plot(α, R_corr, fmt='-', zorder=5, label="Messwerte (korrigiert)")
                if 'R_corr_diff' in config['show']:
                    plt2.plot(α, R_corr_diff, fmt='-', zorder=5, label="Messwerte (um diffuse korrigiert)")
                if 'R_corr[peaks]' in config['show']:
                    plt2.plot(α[peaks], R_corr[peaks], fmt='xk', zorder=5, label="Peaks")
                if 'R_corr_diff[peaks]' in config['show']:
                    plt2.plot(α[peaks], R_corr_diff[peaks], fmt='xk', zorder=5, label="Peaks")

                if 'par_input' in config['show']:
                    plt2.plot(α, par_input, '-', zorder=5, label="Theoriekurve (rau)")
                if 'r13_input' in config['show']:
                    plt2.plot(α, r13_input, '--', label="Theoriekurve (Fresnel)")
                if 'par_fit' in config['show']:
                    plt2.plot(α, par_fit, '-', zorder=5, label="Theoriekurve (rau)") # label → Fit
                if 'r13_fit' in config['show']:
                    plt2.plot(α, r13_fit, '--', label="Theoriekurve (Fresnel)") # label → Fit
                if 'r13_fit_glatt' in config['show']:
                    plt2.plot(α, r13_fit_glatt, '--', label="Fresnelreflektivität")
                if 'par_fit_glatt' in config['show']:
                    plt2.plot(α, par_fit_glatt, '-', label="Theoriekurve (glatt)")

                if 'α_g' in config['show']:
                    plt.axvline(α_g.to('°'), color='C2', linestyle='--', label="$α_g$")
                if 'α_c_PS_fit' in config['show']:
                    plt.axvline(α_c_PS_fit.to('°'), color='C3', linestyle='--',
                                label=r"$α_\mathrm{c, PS}$")  # TODO label=r"$α_\text{c, PS}$"
                if 'α_c_Si_fit' in config['show']:
                    plt.axvline(α_c_Si_fit.to('°'), color='C4', linestyle='--',
                                label=r"$α_\mathrm{c, Si}$")  # TODO label=r"$α_\text{c, Si}$"
                if 'fit_mask' in config['show']:
                    plt.axvspan(*tools.bounds(α[fit_mask]), color='C1', alpha=0.2, label="Fitbereich")

            if config.get('cut_plot') == "little":
                # cut a little
                plt.xlim(right=1.5)
                plt.ylim(bottom=1E-6)  # COULDDO: No idea why this is necessary
            if config.get('cut_plot') == "lot":
                # cut a lot
                # plt.xlim(0.1, 1.0)
                plt.xlim(0.0, 1.0)
                # plt.ylim(1E-5, 1E0)

            plt.yscale('log')
            plt.grid()
            # plt.legend(fontsize=8)
            plt.legend()
            plt.tight_layout()
            if tools.BUILD:
                plt.savefig(f"build/plt/{name}_{config['name']}.pdf")
            plt.show()

            return plt

        for plot_config in plot_configs:
            plot_schichtdicke(plot_config)
