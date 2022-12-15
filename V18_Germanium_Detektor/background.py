import matplotlib.pyplot as plt
from common import load_spe

import tools


def main():
    x, N, T = load_spe("data/2022-11-28/5_Background.Spe", subtract_background=False)

    # Plot: Untergrund
    plt.figure()
    with tools.plot_context(plt, 'dimensionless', 'dimensionless', "x", "N") as plt2:
        plt2.plot(x, N, label="Messwerte")
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    if tools.BUILD:
        plt.savefig("build/plt/background.pdf")
    if tools.PLOTS:
        plt.show()
