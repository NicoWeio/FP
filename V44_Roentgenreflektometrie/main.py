# import generate_table
import code.detektorscan as detektorscan
import code.rockingscan as rockingscan
import code.schichtdicke as schichtdicke
import code.zscan as zscan

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
        'T': ureg('1 s'),
    },
    '2_z1': {
        'x_units': ureg.mm,
        'T': ureg('1 s'),
    },
    '3_x': {
        'x_units': ureg.mm,
        'T': ureg('1 s'),
    },
    '4_rocking1': {
        'x_units': ureg.degree,
        'T': ureg('1 s'),
    },
    '5_z2': {
        'x_units': ureg.mm,
        'T': ureg('1 s'),
    },
    '6_rocking2': {
        'x_units': ureg.degree,
        'T': ureg('1 s'),
    },
    '7_z3': {
        'x_units': ureg.mm,
        'T': ureg('1 s'),
    },
    '8_reflektivitaet': {
        'x_units': ureg.degree,
        'T': ureg('5 s'),
    },
    '9_diffus': {
        'x_units': ureg.degree,
        'T': ureg('5 s'),
    },
}


def load_scan(name, x_units, T):
    x, N = np.genfromtxt(f"data/{name}.txt", unpack=True)  # skip_header=1

    # Poisson-Fehler
    N = unp.uarray(N, np.sqrt(N))

    x *= x_units
    N *= ureg.dimensionless

    I = N / T

    return x, I


# scan_name = '1_detektor'
# detektorscan.main(scan_name, *load_scan(scan_name, **SCANS[scan_name]), ureg=ureg)

# scan_name = '2_z1'
# zscan.main(scan_name, *load_scan(scan_name, **SCANS[scan_name]), ureg=ureg)

schichtdicke.main(
    "schichtdicke",
    load_scan('8_reflektivitaet', **SCANS['8_reflektivitaet']),
    load_scan('9_diffus', **SCANS['9_diffus']),
    ureg=ureg,
)
