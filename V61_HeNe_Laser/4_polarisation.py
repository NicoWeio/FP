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
ɑ *= ureg('°')
I *= ureg.microwatt

# Fit berechnen
# params = tools.pint_curve_fit(calc_I, ɑ, I, (ureg.deg, ureg.microwatt)) # zu niedrig…

params = tools.pint_curve_fit(calc_I, ɑ.to('rad'), I.to('µW'), (1*ureg.rad, 1*ureg.microwatt, 1*ureg.microwatt), p0=(0 * ureg.rad, max(I).to('µW'), 0 * ureg.microwatt))

print(f"params: {params}")
nominal_params = [tools.nominal_value(p) for p in params]
print(f"nominal_params: {nominal_params}")
ɑ_linspace = tools.linspace(ureg('0 °'), ureg('360 °'), 500)


with tools.plot_context(plt, '°', 'microwatt', r'\alpha', 'I') as plt2:
    plt2.plot(ɑ, I, 'x', zorder=5, label='Messwerte')
    plt2.plot(ɑ_linspace, calc_I(ɑ_linspace, *nominal_params), label='Fit')

plt.gca().xaxis.set_ticks(np.arange(0, 360 * 9/8, 360/8))
plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig('build/plt/4_polarisation.pdf')
plt.show()
