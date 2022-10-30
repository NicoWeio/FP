# import generate_table
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# import generate_table
import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
T = ureg('5 s')

# █ Daten einlesen
α, N = np.genfromtxt("data/4_reflektivitätsscan.txt", unpack=True)  # skip_header=1

# Poisson-Fehler
N = unp.uarray(N, np.sqrt(N))

α *= ureg.degree
N *= ureg.dimensionless

I = N / T


# █ Tabelle generieren
# generate_table.generate_table_pint(
#     'build/tab/1_koinzidenz.tex',
#     (r't_\text{diff}', ureg.degree, Δt),
#     ('I', ureg.second**-1, tools.nominal_values(I)),  # TODO: make ufloats work with tables (again)
# )

q = 4 * np.pi / λ * np.sin(np.pi / 180 * α)  # TODO: Blind übernommen aus @Mampfzwerg

peaks, peak_props = sp.signal.find_peaks(tools.nominal_values(I).to('1/s').m, height=(1E2, 1E4), prominence=5)
Δα_mean = np.mean(np.diff(α[peaks].to('rad').m)) * ureg.rad
Δq_mean = np.mean(np.diff(q[peaks].to('1/m').m)) * ureg['1/m']
# d_estim_a = 2*np.pi / Δq_mean # anderes Resultat als bei @Mampfzwerg
d_estim_b = λ / (2 * Δα_mean) # korrekte Daten bei @Mampfzwerg
print(f"Δα_mean = {Δα_mean}")
print(f"Δq_mean = {Δq_mean}")
# print(f"d_estim_a = {d_estim_a.to('m')}")
print(f"d_estim_b = {d_estim_b.to('m')}")

    # █ Plot
    # α_linspace = tools.linspace(*tools.bounds(α), 1000)

    plt.figure()
    # TODO: Doppelachse mit Intensität und Reflektivität?
    with tools.plot_context(plt, '1/m', '1/s', "q", "I") as plt2:
        plt2.plot(q, I_refl, fmt='-', zorder=5, label="Messwerte")  # oder 'x--'?
        plt2.plot(q, I_diff, fmt='-', zorder=5, label="Messwerte (diffuse)")
        plt2.plot(q, I, fmt='-', zorder=5, label="Messwerte (ohne diffuse)")

    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plt/{name}.pdf")
    plt.show()
