# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.optimize import curve_fit


# In[2]:


def linear(x, m, b):
    return m*x+b


def b_helm(r, N, I):
    return 8 * const.mu_0 * I * N / (np.sqrt(125) * r)


def g(m):
    return 4 * np.pi * const.m_e / (m * const.e)


def kernspin(g_f, J):
    return J * (2/g_f - 1)


def quad_zeemann(g_f, B, E_HFS, m_f):
    return g_f**2 * (const.e * const.h/(4 * np.pi * const.m_e))**2 * B**2 * (1 - 2 * m_f)/E_HFS


# Werte auslesen. rf_freq zunächst in kHz | U1_hor, U1_sweep, U2_hor, U2_sweep alles in V und hor hat 13.8V offset

# In[3]:


rf_freq, U1_hor, U1_sweep, U2_hor, U2_sweep = np.genfromtxt("data.txt", unpack=True)
rf_freq = rf_freq * 1000  # Nun von kHz in Hz umgerechnet
U1_hor = unp.uarray(U1_hor, 0.02)
U2_hor = unp.uarray(U2_hor, 0.02)
U1_sweep = unp.uarray(U1_sweep, 0.01)
U2_sweep = unp.uarray(U2_sweep, 0.01)

#print(U1_hor[1])
#print(noms(U1_hor)[1])
#print(stds(U1_hor)[1])


# In[4]:


#Horizontalspule:
r_hor = 0.1579
N_hor = 154
R_hor = 10/3
#Sweepspule:
r_sweep = 0.1639
N_sweep = 11
R_sweep = 10
#Vertikalspule:
r_ver = 0.11735
N_ver = 20


# In[5]:


#Spannungen in Volt
U1_hor = U1_hor - 13.8
U2_hor = U2_hor - 13.8
#Stromstärken in Ampere
I1_hor = U1_hor/R_hor
I2_hor = U2_hor/R_hor
I1_sweep = U1_sweep/R_sweep
I2_sweep = U2_sweep/R_sweep


# In[6]:


B1_hor = b_helm(r_hor, N_hor, I1_hor)
B1_sweep = b_helm(r_sweep, N_sweep, I1_sweep)
B1_ges = B1_hor + B1_sweep

B2_hor = b_helm(r_hor, N_hor, I2_hor)
B2_sweep = b_helm(r_sweep, N_sweep, I2_sweep)
B2_ges = B2_hor + B2_sweep
print(B2_ges)
#print("Gesamte horizontale Magnetfeldstärke in Tesla")
print(f"Dip1:{max(B1_ges), np.argmax(B1_ges)}, Dip2:{max(B2_ges), np.argmax(B2_ges)}")


# In[7]:


params_1, cov_1 = curve_fit(linear, rf_freq, noms(B1_ges), sigma=stds(B1_ges))
errors1 = np.sqrt(np.diag(cov_1))
params_2, cov_2 = curve_fit(linear, rf_freq, noms(B2_ges), sigma=stds(B2_ges))
errors2 = np.sqrt(np.diag(cov_2))

params_85 = unp.uarray(params_2, errors2)
params_87 = unp.uarray(params_1, errors1)

print(f"m, b von Rb-85: {params_85}")
print(f"m, b von Rb-87: {params_87}")

# Kernspin aus Steigung berechnen

# In[8]:


x = np.linspace(0, 1000000)
plt.plot(x/1000, linear(x, *params_2)*1000000, label=r"${}^{85}$Rb")
plt.errorbar(rf_freq/1000, noms(B2_ges)*1000000, yerr=stds(B2_ges)*1000000, fmt='rx')
plt.plot(x/1000, linear(x, *params_1)*1000000, label=r"${}^{87}$Rb")
plt.errorbar(rf_freq/1000, noms(B1_ges)*1000000, yerr=stds(B1_ges)*1000000, fmt='rx')
plt.ylabel("B-Feld in $\mu$T")
plt.xlabel("RF-Freq in kHz")
plt.legend(loc="best")
plt.savefig("Rb_85_87.pdf")


# In[9]:


#x = np.linspace(0,1000000)
#plt.plot(x/1000, linear(x, *params_1)*1000000, label=r"${}^{87}$Rb")
#plt.errorbar(rf_freq/1000, noms(B1_ges)*1000000,yerr=stds(B1_ges)*1000000, fmt='rx')
#plt.ylabel("B-Feld in $\mu$T")
#plt.xlabel("RF-Freq in kHz")
#plt.legend(loc="best")
#plt.savefig("Rb_87.pdf")

# g-Faktoren berechnen

# In[10]:


g_f1 = g(unp.uarray(params_1[0], errors1[0]))
g_f2 = g(unp.uarray(params_2[0], errors2[0]))
print(f"Rb-85: g_f={g_f2}")
print(f"Rb-87: g_f={g_f1}")


# In[11]:


I_1 = kernspin(g_f1, 1/2)
I_2 = kernspin(g_f2, 1/2)


# In[12]:

print(f"Rb-85: I={I_2}")
print(f"Rb-87: I={I_1}")


# Erdmagnetfeld in vertikaler und horizontaler Richtung bestimmen
#

# In[13]:


#Spannung in Volt
U_vert = unp.uarray(2.26, 0.01)
I_vert = U_vert/R_sweep
B_vert = b_helm(r_ver, N_ver, I_vert)
print("Vertikale Magnetfeldkomponente in muT aus Spannung der vertikalen Spule")
print(B_vert * 1000000)


# In[14]:


B_hor = (unp.uarray(params_1[1], errors1[1]) + unp.uarray(params_2[1], errors2[1]))/2
print("Horizontale Magnetfeldkomponente in muT aus y-Achsenabschnitt des lin. Zusammenhangs zwischen B-Feld unf RF-Frequenz")
print(B_hor * 1000000)


# Vergleichswerte:
# horizontale Magntefledkomponente ~20 Mikrotesla; vertikale Magnetfeldkomponente ~40 Mikrotesla

# Das Verhältnis der beiden Isotope lässt sich anhand des Amplitudenverhältnis der beiden Dips (1. Dip Rb-87 und 2. Dip Rb-85) ablesen. Es beträgt Rb-87 = 0.5 * Rb-85

# Abschätzung des quadratischen Zeemann-Effekts:
#

# In[19]:


E_HFS_87 = 4.53 * 10**(-24)
E_HFS_85 = 2.01 * 10**(-24)
E_QZ_87 = quad_zeemann(g_f1, B1_ges[9], E_HFS_87, 1)
E_QZ_85 = quad_zeemann(g_f2, B2_ges[9], E_HFS_85, 1)

#print("Aufspaltung durch den quadratischen Zeemann-Effekt bei höheren Magnetfeldstärken")
print(f"Rb-87:{abs(E_QZ_87)} J bei {max(B1_ges)*1000000} Mikrotesla, Rb-85: {abs(E_QZ_85)} J bei {max(B2_ges)*1000000} Mikrotesla")


# In[ ]:
