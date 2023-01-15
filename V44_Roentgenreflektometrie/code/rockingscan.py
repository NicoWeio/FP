import matplotlib.pyplot as plt
import numpy as np

import tools


def main(name, α, I, ureg, d_Strahl, D):
    # I_thresh = tools.nominal_value(np.percentile(I, 30))  # NOTE: Nicht gerade robust
    I_thresh = tools.nominal_value(np.percentile(I, 50))  # TODO: Das 30-Perzentil liefert eine kleinere Abweichung von α_g_alt…
    # find the first index where I > I_thresh and the last index where I > I_thresh
    i_min = next(i for i, x in enumerate(I) if x > I_thresh)
    i_max = next(i for i, x in reversed(list(enumerate(I))) if x > I_thresh)
    print(f"I_thresh = {I_thresh}")

    α_g = α[[i_min, i_max]]
    α_g_mean = np.mean(np.abs(α_g))
    print(f"α_g_1 = {α_g[0]:.2f}")
    print(f"α_g_2 = {α_g[1]:.2f}")
    print(f"α_g_mean = {α_g_mean:.2f}")

    # █ alternative Berechnung
    α_g_alt = np.arcsin((d_Strahl / D).to('dimensionless'))
    print(f"α_g_alt = arcsin(d_Strahl / D) = {α_g_alt.to('°'):.2f}")

    print(tools.fmt_compare_to_ref(α_g_mean, α_g_alt, unit='°', name="α_g_mean vs. α_g_alt"))

    # █ Plot
    # α_linspace = tools.linspace(*tools.bounds(α), 1000)

    if tools.PLOTS:
        plt.figure()
        with tools.plot_context(plt, '°', '1/s', r"\alpha", "I") as plt2:
            plt2.plot(α, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?
            # plt.axhline(I_thresh, color='C2', label="Schwellwert")

            plt.axvline(α[i_min], color='C1', label=r"$\alpha_{\text{g}, i}$")
            plt.axvline(α[i_max], color='C1')

            plt.axvline(-α_g_alt, color='C2', label=r"$\pm \alpha_\text{g}'$")
            plt.axvline(+α_g_alt, color='C2')

        plt.grid()
        for yscale in ['linear', 'log']:
            plt.yscale(yscale)
            plt.legend()
            plt.tight_layout()
            if tools.BUILD:
                plt.savefig(f"build/plt/{name}_{yscale}.pdf")
            plt.show()

    return α_g_mean
