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


def series(
    z_input, α_c_Si_input,
    d1_input, d2_input,
    s1_input, s2_input,
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
        z, α.to('rad').m,
        k=k, α_c_Si=α_c_Si,
        d1=d1_input * 1E-6, d2=d2_input * 1E-6,
        s1=s1_input * ureg['Å'], s2=s2_input * ureg['Å'],
        ureg=ureg,
        # rauigkeit=False,
        rauigkeit=True,
    )
    plt.plot(q, par)
    plt.plot(q, r13)
    plt.yscale('log')
    return ()


interact(series,
         z_input=(500, 1000), α_c_Si_input=(0.1, 0.3, 0.01),
         d1_input=(0.1, 1.0, 0.1), d2_input=(0.1, 10, 0.1),
         s1_input=(1, 10), s2_input=(1, 10),
         );
