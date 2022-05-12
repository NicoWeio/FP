from generate_table import generate_table_pint
import matplotlib.pyplot as plt
import numpy as np
import pint
import tools
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()


def calc_I(
    ɑ,  # Parameter
    ɑ0, I_max, I_0,  # Konstanten
):
    """TEM₀₀: Theoretische Intensität in Abhängigkeit vom Abstand r von der Modenmitte."""
    return I_max * np.sin(ɑ + ɑ0)**2 + I_0


ɑ, I = np.genfromtxt('dat/4_polarisation.csv', delimiter=',', skip_header=1, unpack=True)
ɑ *= ureg('°')
I *= ureg.microwatt

# Tabelle erzeugen
generate_table_pint(f'build/tab/4_polarisation.tex', (r'\alpha', ureg.deg, ɑ, 0), ('I', ureg.microwatt, I, 1))

# Fit berechnen
# WICHTIG: ɑ muss in hier manuell rad umgerechnet werden,
# weil die Fit-Funktion das nicht implizit umwandeln kann.
params = tools.pint_curve_fit(
    calc_I, ɑ.to('rad'), I,
    (ureg.rad, ureg.microwatt, ureg.microwatt),
    # bounds=[
    #     # (ureg('0 °'), ureg('360 °')),
    #     # (ureg('-180 °'), ureg('180 °')),
    #     (ureg('-5 °'), ureg('5 °')),
    #     None,
    #     None
    # ],
    p0=[ureg('0 °'), max(I), ureg('0 µW')]
)
params[0].ito('°')
print(f"params: {params}")


nominal_params = [tools.nominal_value(p) for p in params]
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
