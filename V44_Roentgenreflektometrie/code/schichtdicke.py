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
    k, α_c_Si,
    δ1, δ2,
    σ1, σ2,
    ureg,
    rauigkeit=False,
):
    """
    δ_i: Brechungsindex-Korrektur (n = 1 - δ_i)
    σ_i: Rauigkeit der Grenzfläche
    """

    n1 = 1  # Luft
    n2 = 1 - δ1  # Polysterol
    n3 = 1 - δ2  # Silizium

    # ----------------------------

    kz1 = k * np.sqrt(n1**2 - np.cos(α)**2) # removed abs(…)
    kz2 = k * np.sqrt(n2**2 - np.cos(α)**2) # removed abs(…)
    kz3 = k * np.sqrt(n3**2 - np.cos(α)**2) # removed abs(…)
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
    par = np.abs(x1)**2

    # Strecke vor Beginn der Oszillationen auf 1 setzen
    # par[α < α_c_Si] = 1
    # r13[α < α_c_Si] = 1
    return par, r13


def main(name, mess_refl, mess_diff, ureg, d_Strahl, α_g, litdata):
    """
    d_Strahl: Strahlbreite (siehe Z-Scan)
    α_g: Geometriewinkel (siehe Rockingscan)
    """
    assert np.all(mess_refl[0] == mess_diff[0]), "α-Werte stimmen nicht überein"
    α, I_refl = mess_refl
    α, I_diff = mess_diff

    # Korrektur um diffusen Anteil
    I_corr_diff = I_refl - I_diff
    # Korrektur um Geometriefaktor
    G = calc_G(α, D=ureg('20 mm'), d_Strahl=d_Strahl, α_g=α_g)
    G[0] = G[1]  # TODO: Workaround for division by zero
    I_corr_G = I_refl / G
    # Korrektur um beides
    I_corr = I_corr_diff / G

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
    α_c_PS = λ * np.sqrt(litdata['PS']['r_e·ρ'] / 2)
    α_c_Si = λ * np.sqrt(litdata['Si']['r_e·ρ'] / 2)
    # print(f"α_c_PS = {α_c_PS.to('°'):.2f}")
    # print(f"α_c_Si = {α_c_Si.to('°'):.2f}")
    print(tools.fmt_compare_to_ref(α_c_PS, litdata['PS']['α_c'], 'α_c_PS'))
    print(tools.fmt_compare_to_ref(α_c_Si, litdata['Si']['α_c'], 'α_c_Si'))

    # █ Parameter
    parrat_params = {
        'α_c_Si': α_c_Si,
        # Brechungsindizes
        # δ1 = litdata['PS']['δ']  # Polysterol → Amplitude vergrößert + negativer Offset
        # δ2 = litdata['Si']['δ']  # Silizium → Amplitude verkleinert + positiver Offset
        'δ1': 0.7e-6,  # Polysterol → Amplitude vergrößert + negativer Offset
        'δ2': 6.8e-6,  # Silizium → Amplitude verkleinert + positiver Offset
        #
        # Rauigkeit
        'σ1': 10e-10 * ureg.m,  # Polysterol → Amplitude verkleinert bei hinteren Oszillationen
        'σ2': 7e-10 * ureg.m,  # Silizium → Senkung des Kurvenendes und Amplitudenverkleinerung der Oszillationen
        #
        # Schichtdicke
        'z': ureg('860 Å'),  # Schichtdicke | verkleinert Oszillationswellenlänge
    }

    print(tools.fmt_compare_to_ref(parrat_params['δ1'], litdata['PS']['δ'], "δ1"))
    print(tools.fmt_compare_to_ref(parrat_params['δ2'], litdata['Si']['δ'], "δ2"))
    print(tools.fmt_compare_to_ref(parrat_params['z'], d_estim_b, "Schichtdicke (Fit vs. Peak-Dist.)", unit='Å'))

    par, r13 = calc_parratt(
        α.to('rad').m,
        k=k,
        **parrat_params,
        ureg=ureg,
        rauigkeit=True,
    )

    par_glatt, r13_glatt = calc_parratt(
        α.to('rad').m,
        k=k,
        **parrat_params,
        ureg=ureg,
        rauigkeit=False,
    )

    # passe Höhe der Theoriekurve an Messwerte an
    # TODO: poor man's fit
    # theory_correction_factor = np.mean(I / par)
    theory_correction_factor = I_corr[peaks[0]] / par[peaks[0]]
    # theory_correction_factor = I_corr[peaks[-3]] / par[peaks[-3]]
    # theory_correction_factor = np.mean(I_corr[peaks[-3:]] / par[peaks[-3:]])
    # print(f"theory_correction_factor = {theory_correction_factor}")
    # NOTE: Bewusst nicht *=, um …?
    par = par * theory_correction_factor
    par_glatt = par_glatt * theory_correction_factor
    r13 = r13 * theory_correction_factor
    r13_glatt = r13_glatt * theory_correction_factor
    assert par.check('1/s'), "par hat falsche Dimension"

    # if tools.PLOTS:
    if True:
        # █ Plot 1: Messwerte und Korrekturen
        # α_linspace = tools.linspace(*tools.bounds(α), 1000)

        # COULDDO: Doppelachse mit Intensität und Reflektivität?
        with tools.plot_context(plt, '°', '1/s', "α", "I") as plt2:
            plt2.plot(α, tools.nominal_values(I_refl), '-', label="reflektiert")
            plt2.plot(α, tools.nominal_values(I_diff), '-', label="diffuse")
            plt2.plot(α, tools.nominal_values(I_corr_diff), '-', zorder=5, label="korrigiert um diffuse")
            plt2.plot(α, tools.nominal_values(I_corr_G), '-', zorder=5, label="korrigiert um Geometriefaktor")
            plt2.plot(α, tools.nominal_values(I_corr), '-', zorder=5, label="korrigiert")

        plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}_a.pdf")
        plt.show()

    # if True:
    if False:
        plt.plot(α, G, '-', zorder=5, label="G-Faktor")
        # with tools.plot_context(plt, '°', '1/s', "α", "I") as plt2:
        #     plt2.plot(q[1:], tools.nominal_values(np.diff(I)), '-', zorder=5, label="Differenzen")
        # plt.yscale('symlog')
        # plt.xlim(right=2E7)
        # plt.savefig('foo.pdf')
        plt.show()

    # if tools.PLOTS:
    if True:  # TODO
        # █ Plot 2: Fit
        # COULDDO: Doppelachse mit Intensität und Reflektivität?
        plt.clf()
        with tools.plot_context(plt, '°', '1/s', "α", "I") as plt2:
            plt2.plot(α, I_corr, fmt='-', zorder=5, label="Messwerte (korrigiert)")
            plt2.plot(α[peaks], I_corr[peaks], fmt='x', zorder=5, label="Peaks")

            plt2.plot(α, par, fmt='-', zorder=5, label="Theoriekurve (rau)")
            # plt2.plot(α, r13, fmt='--', label="Theoriekurve (Fresnel)")
            plt2.plot(α, r13_glatt, fmt='--', label="Fresnelreflektivität")
            plt2.plot(α, par_glatt, fmt='-', label="Theoriekurve (glatt)")

            plt.axvline(α_g, color='C0', linestyle='--', label="$α_g$")
            plt.axvline(α_c_PS, color='C1', linestyle='--', label="α_c,PS") # TODO label=r"$α_\text{c, PS}$")
            plt.axvline(α_c_Si, color='C2', linestyle='--', label="α_c,Si") # TODO label=r"$α_\text{c, Si}$")

        plt.xlim(right=1.5)

        plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}_b.pdf")
        plt.show()
