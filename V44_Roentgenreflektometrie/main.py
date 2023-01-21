# %%
# %load_ext autoreload
# %autoreload 2

# import generate_table
import code.detektorscan as detektorscan
import code.rockingscan as rockingscan
import code.schichtdicke as schichtdicke
import code.zscan as zscan
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pint
import scipy as sp
import uncertainties.unumpy as unp
from rich.console import Console
from yaml import scan

import tools

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
console = Console()

SCANS = {
    '1_detektor': {
        'x_units': ureg.degree,
        # 'T': ureg('1 s'),
    },
    '2_z1': {
        'x_units': ureg.mm,
        # 'T': ureg('1 s'),
    },
    '3_x': {
        'x_units': ureg.mm,
        # 'T': ureg('1 s'),
    },
    '4_rocking1': {
        'x_units': ureg.degree,
        # 'T': ureg('1 s'),
    },
    '5_z2': {
        'x_units': ureg.mm,
        # 'T': ureg('1 s'),
    },
    '6_rocking2': {
        'x_units': ureg.degree,
        # 'T': ureg('1 s'),
    },
    '7_z3': {
        'x_units': ureg.mm,
        # 'T': ureg('1 s'),
    },
    '8_reflektivitaet': {
        'x_units': ureg.degree,
        # 'T': ureg('5 s'),
    },
    '9_diffus': {
        'x_units': ureg.degree,
        # 'T': ureg('5 s'),
    },
}


def load_scan(name, x_units):
    # data_path = Path(f"data_Mampfzwerg/{name}.txt") # skip_header=1
    data_path = Path(f"data/{name}.UXD")

    raw_data = data_path.read_text()
    T_line = raw_data.splitlines()[28]
    assert T_line.startswith('_STEPTIME=')
    T = float(T_line.split('=')[1].strip()) * ureg.s
    # print(f"T: {T}")

    x, N = np.genfromtxt(str(data_path), unpack=True, skip_header=56)

    # Poisson-Fehler
    N = unp.uarray(N, np.sqrt(N))

    x *= x_units
    N *= ureg.dimensionless

    I = N / T

    return x, I


scan_name = '1_detektor'
print(f"█ {scan_name}")
I_max = detektorscan.main(scan_name, *load_scan(scan_name, **SCANS[scan_name]), ureg=ureg)

scan_name = '2_z1'
print(f"█ {scan_name}")
d_Strahl = zscan.main(scan_name, *load_scan(scan_name, **SCANS[scan_name]), ureg=ureg)

D = ureg('20 mm')  # Probendurchmesser [Versuchsanleitung]

scan_name = '4_rocking1'
print(f"█ {scan_name}")
α_g = rockingscan.main(scan_name, *load_scan(scan_name, **SCANS[scan_name]), ureg=ureg, d_Strahl=d_Strahl, D=D)

LITDATA = {
    # Polysterol:
    'PS': {
        'r_e·ρ': ureg('9.5E10 / cm²'),
        'δ': 3.5E-6 * ureg.dimensionless,
        'μ': ureg('4 / cm'),
        'α_c': ureg('0.153 °'),
    },
    # Silizium:
    'Si': {
        'r_e·ρ': ureg('20E10 / cm²'),
        'δ': 7.6E-6 * ureg.dimensionless,
        'μ': ureg('86 / cm'),
        'α_c': ureg('0.174 °'),
    },
}

# %%
δ1 = 0.9e-6 + 0.0e-6
δ2 = 6.8e-6 - 1.6e-6
# δ1 = LITDATA['PS']['δ']
# δ2 = LITDATA['Si']['δ']

PARRATT_PARAMS = {
    # █ Brechungsindizes:
    # ▒ Korrekturterm
    'δ1': δ1,  # Polysterol → untere Einhüllende / Amplitude vergrößert + negativer Offset
    'δ2': δ2,  # Silizium → Amplitude verkleinert + positiver Offset
    # ▒ Absorption
    # 'β1': 0E-6,  # Polysterol → Amplitude vergrößert
    # 'β2': 0E-6,  # Silizium → Amplitude verkleinert + positiver Offset
    # ↓ vgl. https://raw.githubusercontent.com/bwgro/FP/main/V44/main.pdf bzw. [Parratt]
    'β1': δ1 / 200,
    'β2': δ2 / 40,
    #
    # █ Rauigkeit
    'σ1': (20.0e-10 + 0.0e-10) * ureg.m,  # Polysterol → Amplitude verkleinert bei größeren α
    'σ2': (7.3e-10 + 0.0e-10) * ureg.m,  # Silizium → Senkung des Kurvenendes und Amplitudenverkleinerung der Oszillationen
    #
    # █ Schichtdicke
    'z': ureg('860 Å'),  # Schichtdicke → verkleinert Oszillationswellenlänge
}

PLOT_CONFIGS = [
    {
        # █ Plot 2: Theoriekurven
        'name': 'theoriekurven',
        'show': {
            'par',
            # 'r13',
            'r13_glatt',
            'par_glatt',
        },
        'cut_plot': 'little',
    },
    {
        # █ Plot 3: Fit
        'name': 'fit',
        'show': {
            # 'R_corr_diff', # COULDDO
            'R_corr',
            'R_corr[peaks]',
            'par_scaled',
            # ---
            'α_g',
            # 'α_c_PS', # mit cut_plot='lot' nicht sichtbar
            'α_c_Si',
        },
        # 'cut_plot': 'little',
        'cut_plot': 'lot',
    },
]

print(f"█ Schichtdicke")
print({k: f"{v:.3E}" for k, v in PARRATT_PARAMS.items()})
schichtdicke.main(
    "schichtdicke",
    load_scan('8_reflektivitaet', **SCANS['8_reflektivitaet']),
    load_scan('9_diffus', **SCANS['9_diffus']),
    ureg=ureg,
    d_Strahl=d_Strahl,
    D=D,
    α_g=α_g,
    I_max=I_max,
    litdata=LITDATA,
    parratt_params=PARRATT_PARAMS,
    plot_configs=PLOT_CONFIGS,
)

# %%
