
import generate_table
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from common import (fit_peak, get_age_probe, get_Ω, intensity_to_alpha,
                    load_lara, load_spe, n_most_intense, peak_fits_to_arrays,
                    plot_energyspectrum, ureg)
from scipy.signal import find_peaks
from uncertainties import UFloat, ufloat

import tools


def main(lit_energies, lit_intensities, Z, T, **kwargs):
    """
    lit_energies: Literaturwerte für die Energie der Peaks
    lit_intensities: Literaturwerte für die Intensität der Peaks
    Z: Peakinhalte
    T: Messzeit
    """
    Ω = get_Ω()
    print(f"Ω = {Ω:.4f} = 4π · {(Ω / (4*np.pi)).to('dimensionless'):.4f}")

    age_probe = get_age_probe()
    print(f"Alter Probe: {age_probe.to('s'):.2e} = {age_probe.to('a'):.2f}")

    t_hw = ufloat(13.522, 0.016) * ureg.a  # [lara]
    A_0 = ufloat(4130, 60) * ureg.Bq  # Aktivität am Tag der Herstellung [Versuchsanleitung]
    A = A_0 * unp.exp((-np.log(2) / t_hw * age_probe).to('dimensionless').m)  # Aktivität am Tag der Durchführung
    print(f"A={A.to('Bq'):.2f}")

    # %% Effizienz
    Q = 4*np.pi*Z / (Ω * A * lit_intensities * T)  # Effizienz
    # assert Q.check('dimensionless'), "Q is not dimensionless"
    Q.ito('dimensionless')  # NOTE: Q.check raises a KeyError for some reason, doing this instead → create an Issue?

    # %% Fit: Q(E) – Effizienz in Abhängigkeit von der Energie

    def fit_fn_Q(E, Q_max, exponent):
        # NOTE: E / (1 keV)
        if isinstance(E, ureg.Quantity):
            E = E.to('keV').magnitude
        return Q_max * E**exponent

    Q_max, exponent = tools.pint_curve_fit(
        fit_fn_Q,
        lit_energies, Q,
        (ureg.dimensionless, ureg.dimensionless),
        p0=(50 * ureg.dimensionless, -1 * ureg.dimensionless),
    )
    print(f"Q_max (aka p): {Q_max:.2f}")
    print(f"exponent (aka q): {exponent:.2f}")

    def E_to_Q(E):
        return fit_fn_Q(E, Q_max, exponent)

    # %% Plot: Q(E) – Effizienz in Abhängigkeit von der Energie
    energy_linspace = tools.linspace(*tools.bounds(lit_energies), 100)
    plt.figure()
    with tools.plot_context(plt, 'keV', 'dimensionless', "E", "Q") as plt2:
        plt2.plot(lit_energies, Q, fmt='x', label="Messwerte")
        # COULDDO: Unsicherheiten als Schattierung
        plt2.plot(energy_linspace, fit_fn_Q(energy_linspace, Q_max, exponent), show_yerr=False, label="Fit")
    plt.legend()
    plt.tight_layout()
    if tools.BUILD:
        plt.savefig("build/plt/effizienz.pdf")
    if tools.PLOTS:
        plt.show()

    # %% Tabelle generieren
    # COULDDO: Unsicherheiten
    generate_table.generate_table_pint(
        "build/tab/2_effizienz.tex",
        (r"E_\text{lit}", ureg.kiloelectron_volt, lit_energies),
        (r"W", ureg.percent, lit_intensities),
        (r"Z", ureg.dimensionless, Z, 0),
        (r"Q", ureg.dimensionless, tools.nominal_values(Q), 4),
    )

    return {
        'E_to_Q': E_to_Q,
    }
