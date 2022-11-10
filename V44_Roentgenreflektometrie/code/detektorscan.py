import matplotlib.pyplot as plt
import numpy as np
import tools


def fit_fn(α, I_max, I_0, σ, α_0):
    # ↓ Hier entspricht σ nicht der gängigen Definition
    # return I_max * np.exp(-(σ * (α - α_0))**2) + I_0

    # ↓ Mampfzwerg
    # return (I_max / np.sqrt(2 * np.pi * σ**2)) * np.exp(-((α - α_0)**2) / (2 * σ**2)) + I_0

    # return I_max * np.exp(-((α/σ)**2)) + I_0 # kinda works
    return I_max * np.exp(-((α - α_0)**2) / (2 * σ**2)) + I_0  # works!!!

    # mit uncertainties-Dings
    # return (I_0 / np.sqrt(2 * np.pi * σ**2)) * np.exp((-((α - α_0)**2) / (2 * σ**2)).to('dimensionless')) + I_max


def main(name, α, I, ureg):
    I_max, I_0, σ, α_0 = tools.pint_curve_fit(
        fit_fn,
        # α, I, # TODO
        α, tools.nominal_values(I),
        (1 / ureg.s, 1 / ureg.s, ureg.deg, ureg.deg),
        p0=(tools.nominal_value(max(I)), tools.nominal_value(min(I)), ureg('0.05°'), ureg('0°')),
        # return_p0=True,  # TODO
    )
    print(f"I_max = {I_max}")
    print(f"→ vs. größter Messwert = {max(I)}")
    print(f"I_0 = {I_0}")
    print(f"σ = {σ}")
    print(f"α_0 = {α_0}")

    # Halbwertsbreite
    fwhm_width = 2 * np.sqrt(2 * np.log(2)) * σ
    print(f"fwhm_width = {fwhm_width}")

    # █ Plot
    if tools.PLOTS:
        α_linspace = tools.linspace(*tools.bounds(α), 1000)

        plt.figure()
        with tools.plot_context(plt, '°', '1/s', r"\alpha", "I") as plt2:
            plt2.plot(α, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

            # plt2.plot(α_linspace, fit_fn(α_linspace, *[I_max, I_0, σ, α_0]), label="Fit")
            plt2.plot(α_linspace, fit_fn(α_linspace, *map(tools.nominal_value, [I_max, I_0, σ, α_0])), label="Fit")

            plt2.plot(
                tools.nominal_values(
                    tools.pintify([α_0 - fwhm_width / 2, α_0 + fwhm_width / 2])
                ),
                tools.nominal_values(
                    tools.pintify([(I_max + I_0) / 2]*2)
                ),
                '-', lw=2, label="Halbwertsbreite",
            )

        plt.grid()
        plt.legend()
        plt.tight_layout()
        if tools.BUILD:
            plt.savefig(f"build/plt/{name}.pdf")
        plt.show()
