# %%
from __future__ import print_function

import random

# append parent directory to path
import sys
sys.path.append('..')

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pint
from ipywidgets import fixed, interact, interact_manual, interactive

import schichtdicke


ureg = pint.UnitRegistry()

# %%


def calc_parratt(
    α,
    z,
    k, α_c_Si,
    δ1, δ2,
    σ1, σ2,
    ureg,
    rauigkeit=False,
):

    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(α)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(α)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(α)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2

    if rauigkeit:
        r12 *= np.exp(-2 * kz1 * kz2 * σ1**2)
        r23 *= np.exp(-2 * kz2 * kz3 * σ2**2)
        r13 *= 0  # NOTE: Hierzu hat @Mampfzwerg keine Formel
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2

    # Strecke vor Beginn der Oszillationen auf 1 setzen
    par[α < α_c_Si] = 1
    r13[α < α_c_Si] = 1
    return par, r13


# %%

def series(
    z_input, α_c_Si_input,
    δ1_input, δ2_input,
    σ1_input, σ2_input,
):
    # z = ureg('855 Å')  # Parameter: Schichtdicke | verkleinert Oszillationswellenlänge
    z = z_input * ureg['Å']
    # α_c_Si = ureg('0.210 °')
    α_c_Si = α_c_Si_input * ureg['°']

    α = np.linspace(0, 2.5, 500) * ureg.deg
    λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
    k = 2*np.pi / λ  # Wellenvektor
    q = schichtdicke.α_to_q(α, λ=λ)
    # q = np.linspace(0.0E9, 2.0E9) * ureg['1/m']

    par, r13 = schichtdicke.calc_parratt(
    # par = calc_parratt(
        α.to('rad').m,
        z=z,
        k=k, α_c_Si=α_c_Si,
        δ1=δ1_input * 1E-6, δ2=δ2_input * 1E-6,
        σ1=σ1_input * ureg['Å'], σ2=σ2_input * ureg['Å'],
        ureg=ureg,
        # rauigkeit=False,
        rauigkeit=True,
    )
    plt.plot(q, par)
    plt.plot(q, r13)
    plt.yscale('log')
    pass


interact(series,
         z_input=(500, 1000), α_c_Si_input=(0.1, 0.3, 0.01),
         δ1_input=(0.1, 1.0, 0.1), δ2_input=(0.1, 10, 0.1),
         σ1_input=(1, 10), σ2_input=(1, 10),
         );

# %%
