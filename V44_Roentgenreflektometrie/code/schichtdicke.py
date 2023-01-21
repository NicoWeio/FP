# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# import generate_table
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
    parratt_params,
    plot_configs,
):
    """
    d_Strahl: Strahlbreite (siehe Z-Scan)
    α_g: Geometriewinkel (siehe Rockingscan)
    I_max: Maximale Intensität aus dem Detektorscan
    """
    assert np.all(mess_refl[0] == mess_diff[0]), "α-Werte stimmen nicht überein"
    α, I_refl = mess_refl
    α, I_diff = mess_diff

    # Korrektur um diffusen Anteil
    I_corr_diff = I_refl - I_diff
    # Korrektur um Geometriefaktor
    G = calc_G(α, D=D, d_Strahl=d_Strahl, α_g=α_g)
    G[0] = G[1]  # TODO: Workaround for division by zero
    I_corr_G = I_refl / G
    # Korrektur um beides
    I_corr = I_corr_diff / G

    R_corr_diff = I_corr_diff / I_max
    R_corr = I_corr / I_max

    λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
    k = 2*np.pi / λ  # Wellenvektor
    q = α_to_q(α, λ=λ)

    # █ Schichtdicke bestimmen (Peaks finden)
    # peaks, peak_props = sp.signal.find_peaks(tools.nominal_values(I_corr).to('1/s').m, height=(1E2, 1E4), prominence=5, distance=8)
    # TODO: Fast funktioniert es automatisch. Fast…
    peaks = [70, 81, 91, 100, 111, 119, (131), 144, 155]
    assert len(peaks) > 0, "Keine Peaks gefunden"
    print(f"Peak-Indizes: {peaks}")
    # TODO: add sem/ufloat
    Δα_mean = np.mean(np.diff(α[peaks].to('rad').m)) * ureg.rad
    Δq_mean = np.mean(np.diff(q[peaks].to('1/m').m)) * ureg['1/m']
    # d_estim_a = 2*np.pi / Δq_mean # anderes Resultat als bei @Mampfzwerg
    d_estim_b = λ / (2 * Δα_mean)  # korrekte Daten bei @Mampfzwerg
    print(f"Δα_mean = {Δα_mean}")
    print(f"Δq_mean = {Δq_mean}")
    # print(f"d_estim_a = {d_estim_a.to('m')}")
    print(f"d_estim_b = {d_estim_b.to('nm'):.2f} = {d_estim_b.to('Å'):.1f}")

    # berechenbar aus δ !? (@Mampfzwerg)
    # α_c_PS = ureg('0.068 °')
    # α_c_Si = ureg('0.210 °')

    plt.figure()

    # z = d_estim_b

    # α_c_Si = np.sqrt(2 * δ)

    # α_c_PS = λ * np.sqrt(litdata['PS']['r_e·ρ'] / np.pi)
    # α_c_Si = λ * np.sqrt(litdata['Si']['r_e·ρ'] / np.pi)
    #
    α_c_PS = np.sqrt(2 * parratt_params['δ1']) * ureg.rad  # !?
    α_c_Si = np.sqrt(2 * parratt_params['δ2']) * ureg.rad  # !?
    #
    # print(f"α_c_PS = {α_c_PS.to('°'):.2f}")
    # print(f"α_c_Si = {α_c_Si.to('°'):.2f}")
    #
    # print(tools.fmt_compare_to_ref(α_c_PS, litdata['PS']['α_c'], 'α_c_PS', unit='°'))
    # print(tools.fmt_compare_to_ref(α_c_Si, litdata['Si']['α_c'], 'α_c_Si', unit='°'))

    # █ Parameter
    # TODO: Move back here

    # print(tools.fmt_compare_to_ref(parratt_params['δ1'], litdata['PS']['δ'], "δ1"))
    # print(tools.fmt_compare_to_ref(parratt_params['δ2'], litdata['Si']['δ'], "δ2"))
    # print(tools.fmt_compare_to_ref(parratt_params['z'], d_estim_b, "Schichtdicke (Fit vs. Peak-Dist.)", unit='Å'))

    par, r13 = calc_parratt(
        α.to('rad').m,
        k=k,
        **parratt_params,
        ureg=ureg,
        rauigkeit=True,
    )

    # Scaling of Parratt
    R_corr_plateau_mean = R_corr[(ureg('0.1°') < α) & (α < ureg('0.2°'))].mean()
    print(f"R_corr_plateau_mean = {R_corr_plateau_mean:.3f}")
    # NOTE: par can be assumed to be 1 in this range
    par_scaled = par * tools.nominal_value(R_corr_plateau_mean / 1)


    # --- TEST (WIP): Fit ---
    # @ureg.wraps(ureg.dimensionless, (ureg.rad, ureg.dimensionless))
    # def parrat_fitfn(α, δ1):
    #     return calc_parratt(
    #         α.to('rad').m,
    #         # ↓ pass these first, so they can be overwritten
    #         **parratt_params,
    #         # ↓ overrides
    #         δ1=δ1,
    #         # ↓ the rest
    #         k=k,
    #         ureg=ureg,
    #         rauigkeit=True,
    #     )[0]

    # δ1_fit = tools.curve_fit(
    #     parrat_fitfn,
    #     α.to('rad').m,
    #     tools.nominal_values(R_corr),
    #     p0=litdata['PS']['δ'].m,
    # )
    # print(tools.fmt_compare_to_ref(δ1_fit, litdata['PS']['δ'], "δ1 (Fit)"))
    # --- TEST (WIP): Fit ---

    par_glatt, r13_glatt = calc_parratt(
        α.to('rad').m,
        k=k,
        **parratt_params,
        ureg=ureg,
        rauigkeit=False,
    )

    # █ Plot 1: Messwerte und Korrekturen
    if tools.PLOTS:
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

            if 'par' in config['show']:
                plt2.plot(α, par, '-', zorder=5, label="Theoriekurve (rau)")
            if 'par_scaled' in config['show']:
                plt2.plot(α, par_scaled, '-', zorder=5, label="Theoriekurve (rau)")
            if 'r13' in config['show']:
                plt2.plot(α, r13, '--', label="Theoriekurve (Fresnel)")
            if 'r13_glatt' in config['show']:
                plt2.plot(α, r13_glatt, '--', label="Fresnelreflektivität")
            if 'par_glatt' in config['show']:
                plt2.plot(α, par_glatt, '-', label="Theoriekurve (glatt)")

            if 'α_g' in config['show']:
                plt.axvline(α_g.to('°'), color='C2', linestyle='--', label="$α_g$")
            if 'α_c_PS' in config['show']:
                plt.axvline(α_c_PS.to('°'), color='C3', linestyle='--', label=r"$α_\text{c, PS}$")
            if 'α_c_Si' in config['show']:
                plt.axvline(α_c_Si.to('°'), color='C4', linestyle='--', label=r"$α_\text{c, Si}$")

        if config.get('cut_plot') == "little":
            # cut a little
            plt.xlim(right=1.5)
            plt.ylim(bottom=1E-6)  # COULDDO: No idea why this is necessary
        if config.get('cut_plot') == "lot":
            # cut a lot
            plt.xlim(0.1, 1.0)
            plt.ylim(1E-5, 1E0)

        plt.yscale('log')
        plt.grid()
        # plt.legend(fontsize=8)
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}_{config['name']}.pdf")
        plt.show()

    if tools.PLOTS or True: # TODO
        for plot_config in plot_configs:
            plot_schichtdicke(plot_config)
