# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# import generate_table
import tools


# # Versuchsgrößen
# k = 2*np.pi / λ  # Wellenvektor
# qz = 2*k * np.sin(α)  # Wellenvektorübertrag -> y-Werte der Theoriekurve

# # Parameter des Parratt-Algorithmus

# ----------------------------


def α_to_q(α, λ):
    # TODO: Blind übernommen aus @Mampfzwerg
    q = 4 * np.pi / λ * np.sin(np.pi / 180 * α)
    return q


def calc_G(α, D, d_Strahl, α_g):
    """
    Berechnet den Geometriefaktor G.
    Dieser wird relevant, wenn die Probe nicht die gesamte Strahlbreite ausfüllt (für α < α_g).

    d_Strahl: Strahlbreite
    α_g: Geometriewinkel
    """
    # α_g = ureg('0.28 mm')
    α_g_alt = np.arcsin(d_Strahl / D)
    G = D * np.sin(α) / d_Strahl
    G[α > α_g] = 1
    return G


def calc_parratt(
    z, α,
    k, α_c_Si,
    d1, d2,
    s1, s2,
    ureg,
    rauigkeit=False,
):

    n1 = 1  # Luft
    n2 = 1 - d1  # Polysterol
    n3 = 1 - d2  # Silizium

    # ----------------------------

    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(α)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(α)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(α)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2

    if rauigkeit:
        r12 *= np.exp(-2 * kz1 * kz2 * s1**2)
        r23 *= np.exp(-2 * kz2 * kz3 * s2**2)
        r13 *= 0  # NOTE: Hierzu hat @Mampfzwerg keine Formel
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2

    # Strecke vor Beginn der Oszillationen auf 1 setzen
    par[α < α_c_Si] = 1
    r13[α < α_c_Si] = 1
    return par, r13


def main(name, mess_refl, mess_diff, ureg, d_Strahl, α_g):
    """
    d_Strahl: Strahlbreite (siehe Z-Scan)
    α_g: Geometriewinkel (siehe Rockingscan)
    """
    assert np.all(mess_refl[0] == mess_diff[0]), "x-Werte stimmen nicht überein"
    α, I_refl = mess_refl
    α, I_diff = mess_diff

    # Ausgangswerte
    I = I_refl.copy()
    # Korrektur um diffusen Anteil
    I -= I_diff
    I_corr_diff = I.copy()
    # Korrektur um Geometriefaktor
    G = calc_G(α, D=ureg('20 mm'), d_Strahl=d_Strahl, α_g=α_g)  # TODO: /= oder *=?
    G[0] = G[1]  # TODO: Workaround for division by zero
    # print(f"G = {G}")
    I /= G
    I_corr_G = I_refl / G

    λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
    k = 2*np.pi / λ  # Wellenvektor
    q = α_to_q(α, λ=λ)

    # █ Schichtdicke bestimmen (Peaks finden)
    peaks, peak_props = sp.signal.find_peaks(tools.nominal_values(I).to('1/s').m, height=(1E2, 1E4), prominence=5)
    assert len(peaks) > 0, "Keine Peaks gefunden"
    # TODO: add sem/ufloat
    Δα_mean = np.mean(np.diff(α[peaks].to('rad').m)) * ureg.rad
    Δq_mean = np.mean(np.diff(q[peaks].to('1/m').m)) * ureg['1/m']
    # d_estim_a = 2*np.pi / Δq_mean # anderes Resultat als bei @Mampfzwerg
    d_estim_b = λ / (2 * Δα_mean)  # korrekte Daten bei @Mampfzwerg
    print(f"Δα_mean = {Δα_mean}")
    print(f"Δq_mean = {Δq_mean}")
    # print(f"d_estim_a = {d_estim_a.to('m')}")
    print(f"d_estim_b = {d_estim_b.to('m')}")

    # berechenbar aus δ !? (@Mampfzwerg)
    α_c_PS = ureg('0.068 °')
    α_c_Si = ureg('0.210 °')

    plt.figure()

    # z = d_estim_b

    # █ Parameter
    # Brechungsindizes
    d1 = 0.7e-6  # Polysterol -> Amplitude vergrößert + negativer Offset
    d2 = 6.7e-6  # Silizium -> Amplitude vergkleinert + positiver Offset
    #
    # Rauigkeit
    s1 = 7.9e-10 * ureg.m  # Polysterol -> Amplitude verkleinert bei hinteren Oszillationen
    s2 = 5.7e-10 * ureg.m  # Silizium -> Senkung des Kurvenendes und  Amplitudenverkleinerung der Oszillationen
    #
    # Schichtdicke
    z = ureg('855 Å')  # Parameter: Schichtdicke | verkleinert Oszillationswellenlänge

    par, r13 = calc_parratt(
        z, α.to('rad').m,
        k=k, α_c_Si=α_c_Si,
        d1=d1, d2=d2,
        s1=s1, s2=s2,
        ureg=ureg,
        rauigkeit=True,
    )

    # passe Höhe der Theoriekurve an Messwerte an
    # TODO: poor man's fit
    # theory_correction_factor = np.mean(I / par)
    # theory_correction_factor = I[peaks[-3]] / par[peaks[-3]]
    theory_correction_factor = np.mean(I[peaks[-3:]] / par[peaks[-3:]])
    print(f"theory_correction_factor = {theory_correction_factor}")
    par = par * theory_correction_factor
    r13 = r13 * theory_correction_factor
    assert par.check('1/s'), "par hat falsche Dimension"

    # if tools.PLOTS:
    if True:  # TODO
        # █ Plot 1: Messwerte und Korrekturen
        # α_linspace = tools.linspace(*tools.bounds(α), 1000)

        # TODO: Doppelachse mit Intensität und Reflektivität?
        plt.figure()
        with tools.plot_context(plt, '1/m', '1/s', "q", "I") as plt2:
            plt2.plot(q, tools.nominal_values(I_refl), 'x', label="Messwerte (reflektiert)")
            plt2.plot(q, tools.nominal_values(I_diff), '.', label="Messwerte (diffuse)")
            # plt2.plot(q, tools.nominal_values(I_corr_diff), '-', zorder=5, label="Messwerte (korrigiert um diffuse)")
            # plt2.plot(q, tools.nominal_values(I_corr_G), '-', zorder=5, label="Messwerte (korrigiert um Geometriefaktor)")
            plt2.plot(q, tools.nominal_values(I), '-', zorder=5, label="Messwerte (korrigiert)")

        plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}_a.pdf")
        plt.show()

        # █ Plot 2: Fit
        # TODO: Doppelachse mit Intensität und Reflektivität?
        plt.figure()
        with tools.plot_context(plt, '1/m', '1/s', "q", "I") as plt2:
            plt2.plot(q, I, fmt='.', zorder=5, label="Messwerte (korrigiert)")  # oder 'x--'?
            plt2.plot(q[peaks], I[peaks], fmt='x', zorder=5, label="Peaks")

            plt2.plot(q, par, fmt='-', label="Theoriekurve (Parratt)")
            plt2.plot(q, r13, fmt='--', label="Theoriekurve (Fresnel)")

            plt.axvline(α_to_q(α_g, λ).to('1/m'), color='C0', linestyle='--', label="$α_g$")
            plt.axvline(α_to_q(α_c_PS, λ).to('1/m'), color='C1', linestyle='--', label="$α_c$ (PS)")
            plt.axvline(α_to_q(α_c_Si, λ).to('1/m'), color='C2', linestyle='--', label="$α_c$ (Si)")

        plt.xlim(0, 2E7)
        # plt.ylim(tools.nominal_value(min(I)), 1E2)
        plt.ylim(bottom=1E2)

        plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}_b.pdf")
        plt.show()
