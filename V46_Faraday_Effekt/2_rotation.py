import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy.constants as const

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
#lamb = lamb * 10**(-6)

plt.plot(lamb**2, theta_a, 'x', label=r'Messwerte $N=\SI[per-mode=reciprocal]{2.8e18}{\per\cubic\centi\meter}$', color='blue')
plt.plot(lamb**2, theta_b, 'x', label=r'Messwerte $N=\SI[per-mode=reciprocal]{1.2e18}{\per\cubic\centi\meter}$', color='orange')
plt.plot(lamb**2, theta_c, 'x', label=r'Messwerte rein', color='limegreen')
plt.xlabel(r'$\lambda^2\mathbin{/}\si{\micro\meter\squared}$')
plt.ylabel(r'$\theta\mathbin{/}\si{\radian}$')
plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig('build/plt/rotation.pdf')
plt.close()

# Fit der Theta-Differenz

theta_norm_diff1 = theta_a_normiert - theta_c_normiert
theta_norm_diff2 = theta_b_normiert - theta_c_normiert

#u = np.linspace(np.min(df2["Filter [mikro m]"] ** 2), np.max(df2["Filter [mikro m]"] ** 2))

def theta_fit(a, x):
    return a * x

# Probe 1
param1, cov1 = curve_fit(theta_fit, lamb**2, theta_norm_diff1)
err1 = np.sqrt(np.diag(cov1))
# Probe 2
param2, cov2 = curve_fit(theta_fit, lamb**2, theta_norm_diff2)
err2 = np.sqrt(np.diag(cov2))

print('a_1 = {:.3} $\pm$ {:.4}'.format(param1[0], err1[0]))
print('a_2 = {:.3} $\pm$ {:.4}'.format(param2[0], err2[0]))

plt.plot(lamb**2, theta_norm_diff1, 'x', label=r'Probe $N=\SI[per-mode=reciprocal]{2.8e18}{\per\cubic\centi\meter}$', color='blue')
plt.plot(lamb**2, theta_norm_diff2, 'x', label=r'Probe $N=\SI[per-mode=reciprocal]{1.2e18}{\per\cubic\centi\meter}$', color='orange')
plt.plot(lamb**2, theta_fit(param1[0], lamb**2), label=r'Fit Differenz Probe $N=\SI[per-mode=reciprocal]{2.8e18}{\per\cubic\centi\meter}$', color='blue')
plt.plot(lamb**2, theta_fit(param2[0], lamb**2), label=r'Fit Differenz Probe $N=\SI[per-mode=reciprocal]{2.8e18}{\per\cubic\centi\meter}$', color='orange')
plt.xlabel(r'$\lambda^2\mathbin{/}\si{\micro\meter\squared}$')
plt.ylabel(r'$\symup{\Delta}\theta\mathbin{/}\si{\radian\per\meter}')
plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig('build/plt/fit_rotation.pdf')

# Berechnung der effektiven Masse

eps_0 = const.epsilon_0
N1 = 2.8e18#*(1/100)**(-3)
N2 = 1.2e18#*(1/100)**(-3)
e0 = const.elementary_charge
e_m = const.electron_mass
B = 405e-3
c = const.speed_of_light
a_1 = ufloat(param1[0], cov1)*1e6
a_2 = ufloat(param2[0], cov2)*1e6

lamb_mean = np.mean(lamb)
print(f'lambda Mittelwert={lamb_mean}')

# lamb_mean = 1.904
n = 3.4724 #bei lambda = 1.078 micrometer

#m1 = unp.sqrt(e0**3 * B * N1 / (8 * np.pi**2 * eps_0 * c**3 * a_1 * n)) / e_m
m_1 = unp.sqrt((e0**3 / (8 * np.pi**2 * eps_0 * c**3)) * 1/a_1 * (N1*B/n)) / e_m
m_2 = unp.sqrt((e0**3 / (8 * np.pi**2 * eps_0 * c**3)) * 1/a_2 * (N2*B/n)) / e_m
print(f'Effektive Masse Probe 1: {m_1}')
print(f'Effektive Masse Probe 2: {m_2}')

# Abweichung berechnen

m_theorie = 0.067

abw_1 = np.abs(m_1 - m_theorie) / m_theorie * 100
abw_2 = np.abs(m_2 - m_theorie) / m_theorie * 100

print(f'Abweichung Probe 1 in Prozent: {abw_1}')
print(f'Abweichung Probe 2 in Prozent: {abw_2}')