# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
from rich.console import Console
import tools
import uncertainties.unumpy as unp
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

# █ Daten einlesen
t, channel = np.genfromtxt("2_mca.csv", unpack=True, delimiter=",", skip_header=1)
t *= ureg.microsecond
channel *= ureg.dimensionless


# █ Tabelle generieren
# TODO

# █ lineare Regression
m, b = tools.linregress(channel, t)
print(f"{(m,b)=}")

# █ Plot
channel_linspace = tools.linspace(*tools.bounds(channel))
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'µs', r'\text{channel}', r'\mathrm{\Delta} t') as plt2:
# with tools.plot_context(plt, 'dimensionless', 'µs', r'channel', r't') as plt2: # TODO
with tools.plot_context(plt, 'dimensionless', 'µs', "channel") as plt2:  # TODO
    plt2.plot(channel, t, 'x', zorder=5, label=r"Messwerte zu $^{85}$Rb")
    plt2.plot(channel_linspace, tools.nominal_values(
        m*channel_linspace + b), label=r"Ausgleichsgerade zu $^{85}$Rb")

plt.legend()
plt.tight_layout()
plt.savefig("build/plt/2_mca.pdf")
plt.show()
