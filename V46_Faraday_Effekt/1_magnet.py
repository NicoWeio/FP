import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

z, B = np.genfromtxt('data/magnet.csv', unpack=True)

plt.plot(z, B, 'x', label='Messwerte')
plt.axhline(y=405, color='orange', label= r'Maximale Magnetfeldstärke $B = \SI{405}{\milli\tesla}$')
plt.xlabel(r'$z/\si{\milli\meter}$')
plt.ylabel(r'$B/\si{\milli\tesla}$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("build/plt/magnet.pdf")