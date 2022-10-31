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

# Brechungsindizes
d1 = 0.7e-6  # Polysterol -> Amplitude vergrößert + negativer Offset
d2 = 6.7e-6  # Silizium -> Amplitude vergkleinert + positiver Offset
#
n1 = 1  # Luft
n2 = 1 - d1  # Polysterol
n3 = 1 - d2  # Silizium

# # Rauigkeit
# s1 = 7.9e-10 * ureg.m  # Polysterol -> Amplitude verkleinert bei hinteren Oszillationen
# s2 = 5.7e-10 * ureg.m  # Silizium -> Senkung des Kurvenendes und  Amplitudenverkleinerung der Oszillationen

# # Schichtdicke
# z = 855e-10 * ureg.m  # verkleinert Oszillationswellenlänge

# ----------------------------


def α_to_q(α, λ):
    # TODO: Blind übernommen aus @Mampfzwerg
    q = 4 * np.pi / λ * np.sin(np.pi / 180 * α)
    return q


def calc_G(α, D, d_Strahl, α_g):
    """
    Berechnet den Geometriefaktor G.
    Dieser wird relevant, wenn die Probe nicht die gesamte Strahlbreite ausfüllt (für α < α_g).
    """
    # α_g = ureg('0.28 mm')
    α_g_alt = np.arcsin(d_Strahl / D)
    G = D * np.sin(α) / d_Strahl
    G[α > α_g] = 1
    return G

# █ Theoriekurve: Fresnelreflektivität Si


def parratt2(z, ai, k):
    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(ai)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(ai)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(ai)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2
    # Strecke vor Beginn der Oszillationen auf 1 setzen

    cutoff_i = 100  # TODO: anpassen (@Mampfzwerg: 296)
    par[:cutoff_i] = 1
    r13[:cutoff_i] = 1
    return par, r13


def main(name, mess_refl, mess_diff, ureg, d_Strahl, α_g):
    assert np.all(mess_refl[0] == mess_diff[0]), "x-Werte stimmen nicht überein"
    α, I_refl = mess_refl
    α, I_diff = mess_diff

    # Ausgangswerte
    I = I_refl
    # Korrektur um diffusen Anteil
    I -= I_diff
    # Korrektur um Geometriefaktor
    G = calc_G(α, D=ureg('20 mm'), d_Strahl=d_Strahl, α_g=α_g)  # TODO: /= oder *=?
    G[0] = G[1]  # TODO: Workaround for division by zero
    # print(f"G = {G}")
    I /= G

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

    # z = 855e-10  # Parameter: Schichtdicke | verkleinert Oszillationswellenlänge
    z = d_estim_b
    par, r13 = parratt2(z, α.to('rad').m, k=k)

    # berechenbar aus δ !? (@Mampfzwerg)
    α_c_PS = ureg('0.068 °')
    α_c_Si = ureg('0.210 °')

    # █ Plot
    # α_linspace = tools.linspace(*tools.bounds(α), 1000)

    plt.figure()
    # TODO: Doppelachse mit Intensität und Reflektivität?
    with tools.plot_context(plt, '1/m', '1/s', "q", "I") as plt2:
        plt2.plot(q, I, fmt='.', zorder=5, label="Messwerte (korrigiert)")  # oder 'x--'?
        plt2.plot(q[peaks], I[peaks], fmt='x', zorder=5, label="Peaks")

        # plt2.plot(q, I_diff, fmt='-', label="Messwerte (diffuse)")
        # plt2.plot(q, I_refl, fmt='-', label="Messwerte (reflektiert; roh)")

        plt2.plot(q, par * max(I), fmt='-', label="Theoriekurve (Parratt)")
        plt2.plot(q, r13 * max(I), fmt='--', label="Theoriekurve (Fresnel)")

        plt.axvline(α_to_q(α_g, λ).to('1/m'), color='C0', linestyle='--', label="$α_g$")
        plt.axvline(α_to_q(α_c_PS, λ).to('1/m'), color='C1', linestyle='--', label="$α_c$ (PS)")
        plt.axvline(α_to_q(α_c_Si, λ).to('1/m'), color='C2', linestyle='--', label="$α_c$ (Si)")

    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/{name}.pdf")
    plt.show()
