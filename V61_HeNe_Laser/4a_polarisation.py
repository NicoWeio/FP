import scipy as sp
import tools
import matplotlib.pyplot as plt
import numpy as np
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()


def calc_I(
    ɑ,  # Parameter
    ɑ0, I_max, I_0,  # Konstanten
):
    """TEM₀₀: Theoretische Intensität in Abhängigkeit vom Abstand r von der Modenmitte."""
    # return I0 * np.cos(ɑ + ɑ0)**2
    return I_max * np.sin(ɑ + ɑ0)**2 + I_0


ɑ, I = np.genfromtxt('dat/4_polarisation.csv', delimiter=',', skip_header=1, unpack=True)
ɑ = np.deg2rad(ɑ)
# ɑ *= ureg('°')
# I *= ureg.microwatt

# Fit berechnen
params, pcov = sp.optimize.curve_fit(
    calc_I, ɑ, I,
    # p0=(0, max(I), 0),
    # bounds=[
    #     (0, 2*np.pi),
    #     (0, 2*max(I)),
    #     (-2*max(I), 2*max(I))
    # ]
    # bounds=(
    #     [0, 0.5*max(I), -2*max(I)],
    #     [2*np.pi, 2*max(I), 2*max(I)]
    # )
)

ɑ_linspace = np.linspace(min(ɑ), max(ɑ), 500)


plt.plot(ɑ, I, 'x', zorder=5, label='Messwerte')
plt.plot(ɑ_linspace, calc_I(ɑ_linspace, *params), label='Fit')

# plt.gca().xaxis.set_ticks(np.arange(0, 360 * 9/8, 360/8))
plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig('build/plt/4_polarisation.pdf')
plt.show()
