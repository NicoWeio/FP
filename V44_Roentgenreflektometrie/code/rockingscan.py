import matplotlib.pyplot as plt
import numpy as np

import tools


def main(name, α, I, ureg):
    I_thresh = tools.nominal_value(20 * min(I))  # TODO: optimize
    # find the first index where I > I_thresh and the last index where I > I_thresh
    i_min = next(i for i, x in enumerate(I) if x > I_thresh)
    i_max = next(i for i, x in reversed(list(enumerate(I))) if x > I_thresh)

    α_g_mean = np.mean(np.abs(α[[i_min, i_max]]))
    print(f'α_g_mean = {α_g_mean}')

    # █ alternative Berechnung
    d_0 = ureg('0.28 mm')  # Strahlbreite (@Mampfzwerg)
    D = ureg('20 mm')  # Probendicke (@Mampfzwerg)
    α_g_alt = np.arcsin(d_0 / D)
    print(f'α_g_alt = {α_g_alt.to("°")}')

    print(tools.fmt_compare_to_ref(α_g_mean, α_g_alt, unit='°'))

    return α_g_mean

    # █ Plot
    # α_linspace = tools.linspace(*tools.bounds(α), 1000)

    plt.figure()
    with tools.plot_context(plt, '°', '1/s', r"\alpha", "I") as plt2:
        plt2.plot(α, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

        plt.axhline(I_thresh, color='C1', alpha=0.5, zorder=0, label="Schwellwert")
        plt.axvline(α[i_min], color='C2', alpha=0.5, zorder=0, label="Geometriewinkel")
        plt.axvline(α[i_max], color='C2', alpha=0.5, zorder=0)

    plt.yscale('log')  # TODO
    plt.grid()
    plt.legend()
    plt.tight_layout()
    if tools.BUILD:
        plt.savefig(f"build/plt/{name}.pdf")
    plt.show()
