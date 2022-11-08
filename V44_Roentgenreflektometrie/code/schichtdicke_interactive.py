# %%
from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import matplotlib.pyplot as plt, random

# append parent directory to path
import sys
sys.path.append('..')

import schichtdicke

import numpy as np
import pint
ureg = pint.UnitRegistry()

def series(dots, colr):
    # z = ureg('855 Å')  # Parameter: Schichtdicke | verkleinert Oszillationswellenlänge
    z = ureg('855 Å')  # Parameter: Schichtdicke | verkleinert Oszillationswellenlänge
    α_c_Si = ureg('0.210 °')

    α = np.linspace(0, 2.5, 500) * ureg.deg
    λ = ureg('1.54 Å')  # ? (@Mampfzwerg)
    k = 2*np.pi / λ  # Wellenvektor
    q = schichtdicke.α_to_q(α, λ=λ)
    # q = np.linspace(0.0E9, 2.0E9) * ureg['1/m']


    par, r13 = schichtdicke.calc_parratt(z, α.to('rad').m, k=k, α_c_Si=α_c_Si, ureg=ureg, rauigkeit=True)
    plt.plot(q, par)
    return();
interact(series, dots=(1,100,1), colr=["red","orange","brown"]);
