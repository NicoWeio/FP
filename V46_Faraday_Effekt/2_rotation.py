import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

# Messwerte der höher-dotierten Probe
lamb, grad1, gradminute1, grad2, gradminute2 = np.genfromtxt('data/probe_1.csv', unpack=True)

theta_a1 = grad1 + 1/60 * gradminute1
print(f'Theta1 Probe 1: {np.round(theta_a1,2)}')
theta_a2 = grad2 + 1/60 * gradminute2
print(f'Theta2 Probe 1: {np.round(theta_a2,2)}')

theta_a = 1/2 * np.abs(theta_a1 - theta_a2)
theta_a = theta_a * np.pi/180
print(f'Theta Probe 1: {np.round(theta_a,2)}')

L1 = 1.296e-3
theta_a_normiert = theta_a/L1
print(f'Theta Probe 1 normiert: {np.round(theta_a_normiert,2)}')

# Messwerte der weniger-dotierten Probe
lamb, grad1, gradminute1, grad2, gradminute2 = np.genfromtxt('data/probe_2.csv', unpack=True)

theta_b1 = grad1 + 1/60 * gradminute1
print(f'Theta1 Probe 2: {np.round(theta_b1,2)}')
theta_b2 = grad2 + 1/60 * gradminute2
print(f'Theta2 Probe 2: {np.round(theta_b2,2)}')

theta_b = 1/2 * np.abs(theta_b1 - theta_b2)
theta_b = theta_b * np.pi/180
print(f'Theta Probe 2: {np.round(theta_b,2)}')

L2 = 1.36e-3
theta_b_normiert = theta_b/L2
print(f'Theta Probe 2 normiert: {np.round(theta_b_normiert,2)}')

# Messwerte der reinen Probe
lamb, grad1, gradminute1, grad2, gradminute2 = np.genfromtxt('data/probe_3.csv', unpack=True)

theta_c1 = grad1 + 1/60 * gradminute1
print(f'Theta1 Probe 3: {np.round(theta_c1,2)}')
theta_c2 = grad2 + 1/60 * gradminute2
print(f'Theta2 Probe 3: {np.round(theta_c2,2)}')

theta_c = 1/2 * np.abs(theta_c1 - theta_c2)
theta_c = theta_c * np.pi/180
print(f'Theta Probe 3: {np.round(theta_c,2)}')

L3 = 5.11e-3
theta_c_normiert = theta_c/L3
print(f'Theta Probe 3 normiert: {np.round(theta_c_normiert,2)}')

#Plotten der Messwerte aller drei Proben
plt.plot(lamb**2, theta_a, 'x', label=r'Messwerte $N=\SI[per-mode=reciprocal]{2.8e18}{\per\cubic\centi\meter}$')
plt.plot(lamb**2, theta_b, 'x', label=r'Messwerte $N=\SI[per-mode=reciprocal]{1.2e18}{\per\cubic\centi\meter}$')
plt.plot(lamb**2, theta_c, 'x', label=r'Messwerte rein')
plt.xlabel(r'$\lambda^2/\si{\micro\meter}$')
plt.ylabel(r'$\theta/\si{\radian}$')
plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig('build/plt/rotation.pdf')

# Bestimmung der effektiven Masse

theta_norm_diff1 = theta_a_normiert - theta_c_normiert
theta_norm_diff2 = theta_b_normiert - theta_c_normiert

def theta_fit(a, x):
    return a * x

# Probe 1
#param, cov1 = curve_fit(theta_fit, lamb**2, theta_norm_diff1)
#err1 = np.sqrt(np.diag(cov1))
#
#print('a_1 = {:.3} $\pm$ {:.4}'.format(param[0], err1[0]))
