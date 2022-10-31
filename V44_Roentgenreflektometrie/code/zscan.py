# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import uncertainties.unumpy as unp
from rich.console import Console

# import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()


def main(name, z, I, ureg):
    # TODO: Auto-detect?
    flank_bound_indices = (17, 23)
    flank_bounds = tuple(z[list(flank_bound_indices)])

    # █ Plot
    # z_linspace = tools.linspace(*tools.bounds(z), 1000)

    plt.figure()
    with tools.plot_context(plt, 'mm', '1/s', "z", "I") as plt2:
        plt2.plot(z, I, fmt='x', zorder=5, label="Messwerte")  # oder 'x--'?

        plt.axvspan(*flank_bounds, color='C1', alpha=0.5, zorder=0, label="Strahlbreite")

    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/{name}.pdf")
    plt.show()
