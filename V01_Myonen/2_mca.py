# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import pint
import uncertainties.unumpy as unp
from rich.console import Console

import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

# █ Daten einlesen
t, channel = np.genfromtxt("2_mca.csv", unpack=True, delimiter=",", skip_header=1)
sorting_i = t.argsort()
t = t[sorting_i]
channel = channel[sorting_i]

t *= ureg.microsecond
channel *= ureg.dimensionless


# █ Tabelle generieren
generate_table.generate_table_pint(
    'build/tab/2_mca.tex',
    ('t', ureg.microsecond, t),
    ('channel', ureg.dimensionless, channel, 0),
)

# █ lineare Regression
m, b = tools.linregress(channel, t)
print(f"{(m,b)=}")

t_per_channel = m

# █ Plot
channel_linspace = tools.linspace(*tools.bounds(channel))
plt.figure()
# with tools.plot_context(plt, 'dimensionless', 'µs', r'\text{channel}', r'\mathrm{\Delta} t') as plt2:
# with tools.plot_context(plt, 'dimensionless', 'µs', r'channel', r't') as plt2: # TODO
with tools.plot_context(plt, 'dimensionless', 'µs', "channel", "t") as plt2:  # TODO
    plt2.plot(channel, t, 'x', zorder=5, label=r"Messwerte zu $^{85}$Rb")
    plt2.plot(channel_linspace, tools.nominal_values(
        m*channel_linspace + b), label=r"Ausgleichsgerade zu $^{85}$Rb")

plt.legend()
plt.tight_layout()
plt.savefig("build/plt/2_mca.pdf")
plt.show()
