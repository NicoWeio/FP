import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from scipy import signal as sig
from uncertainties import ufloat
import uncertainties.unumpy as unp

# Versuchsgrößen
l = 1.54e-10  # Wellenlänge
ai = np.pi/180 * np.arange(6e-2, 1.505, 5e-4)  # Winkel -> x−Werte der Theoriekurve
k = 2*np.pi / l  # Wellenvektor
qz = 2*k * np.sin(ai)  # Wellenvektorübertrag -> y-Werte der Theoriekurve

# Parameter des Parratt-Algorithmus

# Brechungsindizes
d1 = 0.7e-6  # Polysterol -> Amplitude vergrößert + negativer Offset
d2 = 6.7e-6  # Silizium -> Amplitude vergkleinert + positiver Offset
#
n1 = 1  # Luft
n2 = 1 - d1  # Polysterol
n3 = 1 - d2  # Silizium

# Rauigkeit
s1 = 7.9e-10  # Polysterol -> Amplitude verkleinert bei hinteren Oszillationen
s2 = 5.7e-10  # Silizium -> Senkung des Kurvenendes und  Amplitudenverkleinerung der Oszillationen

# Schichtdicke
z = 855e-10  # verkleinert Oszillationswellenlänge

#################################################################################################################

# Messdaten
x, y = np.genfromtxt('data/4_reflektivitätsscan.txt', unpack=True)
y = y / (5 * 1636650)  # Normierung auf die 5fache maximale Intensität des Z-Scans


def plotten(phi, I):
    fig, ax1 = plt.subplots()
    #
    q = 2*k * np.sin(np.pi/180 * phi)  # Winkel -> Wellenvektorübertrag
    phi_d, I_d = np.genfromtxt('data/4_difscan.txt', unpack=True)  # Messdaten des diffusen Scans
    I_d = I_d / 1636650  # Normierung
    q_d = 2*k * np.sin(np.pi/180 * phi_d)
    #
    I_d_k = diff_k(I, I_d)  # Korrektur der Messdaten um Diffusen Scan
    I_res = geofaktor()  # Korrektur der Messdaten um Geometriefaktor
    #
    q_min, I_min, peak_array = min(I_res, q)  # findet Oszillationsminima
    d_layer = layer(q, peak_array)  # Bestimmt Schichtdicke
    print(d_layer)
    #
    u = 40  # Ende des Plateaus
    beg = 11  # Anfang sinnvoller Messdaten
    end = 300  # Ende sinnvoller Messdaten
    #
    unkorr = I[beg:end] * 1e-1  # Übersichtshalber werden die Daten nach unten verschoben
    diff_data = I_d[beg:end] * 1e-1
    # Plots
    ax1.plot(q[beg:end], unkorr, linestyle='dashed', linewidth=0.7, color='#FB9D02', label='Messdaten x 0,1')
    ax1.plot(q_d[beg:end], diff_data, linestyle='dashed', linewidth=0.7, color='#4DD30A', label='Diffuse Daten x 0,1')
    ax1.plot(q_d[beg:end], I_d_k[beg:end], linewidth=0.7, color='#AD0EBD', label='Daten um diffuse Daten korrigiert x 0,1')
    ax1.plot(q[beg:end], I_res[beg:end], linewidth=0.7, color='#045DF9', label='Daten um Geometriefaktor korrigiert')
    ax1.plot(q_min, I_min, 'x', color='#06189A', label='Schichtdickenminima')
    #
    phicrit_PS, qcrit_PS = phi_crit(d1)  # Bestimmung der kritischen Winkel
    ax1.axvline(x=qcrit_PS, ymin=0, ymax=1, linestyle='dashed', linewidth=0.7,
                color='#ABB2B9', label=r'$\alpha_c$ für PS')  # Vertikales Plotten des Winkels
    phicrit_Si, qcrit_Si = phi_crit(d2)
    ax1.axvline(x=qcrit_Si, ymin=0, ymax=1, linestyle='dashed', linewidth=0.7, color='#566573', label=r'$\alpha_c$ für Si')
    #
    par = parratt(z)  # Aufruf des Parrat-Algorithmus
    par2, r13 = parratt2(z)
    plt.plot(qz, par, linewidth=0.7, color='#F90421', label='Theoriekurve (raues Si)')
    plt.plot(qz, par2, linewidth=0.7, color='#FA58AC', label='Theoriekurve (glattes Si)', alpha=0.5)
    plt.plot(qz, r13, linewidth=0.7, color='#FA58AC', linestyle='dashed', label='Fresnelrefektivität Si', alpha=0.5)
    #
    #
    plt.yscale('log')
    plt.grid(alpha=0.2)
    ax1.legend(fancybox=True, ncol=1, loc='upper right', fontsize='small')
    ax1.set_xlabel('q (1/m)')
    ax1.set_ylabel('Reflektivität')
    fig.tight_layout()
    # plt.savefig('plot4.pdf')
    plt.show()


def diff_k(I, I_d):
    I_d_k = np.zeros(np.size(I))
    for i in np.arange(np.size(I)):
        I_d_k[i] = 1e-1 * (I[i] - I_d[i])
    return I_d_k


def geofaktor():
    ag = 0.74
    I_res = np.zeros(np.size(x))
    for i in np.arange(np.size(x)):
        if (x[i] <= ag and x[i] > 0):
            I_res[i] = y[i] * np.sin(ag) / np.sin(x[i])
        else:
            I_res[i] = y[i]
    return I_res


def min(I_res, q):
    I_log = np.zeros(np.size(I_res))  # logarithmierte Daten
    for i in np.arange(np.size(I_res)):
        if (I_res[i] != 0):
            I_log[i] = np.log10(I_res[i])
        else:
            I_log[i] = 0
    #
    I_log_n = -1 * I_log  # Negatives Vorzeichen -> Minima werden zu Maxima -> Zur Arbeit mit find_peaks
    peakfinder = sig.find_peaks(I_log_n, prominence=7e-2)
    #
    peak_array = [78, 88, 98, 108, 118, 128, 139, 149, 160, 170, 181, 192, 203, 226]  # ... Manuell angepasst
    #
    peak_q = []
    peak_I = []
    #
    for i in peak_array:
        peak_q.append(q[i])
        peak_I.append(I_res[i])
    #
    return peak_q, peak_I, peak_array


def layer(q, peak_array):
    dist_peaks = []  # Peakabstand
    for i in peak_array:
        if (i < 277):
            dist_peaks.append(q[i+1] - q[i])
        else:
            pass
    mean = np.mean(dist_peaks)
    std = np.std(dist_peaks)
    d = 2*np.pi / ufloat(mean, std)
    return d


def phi_crit(delta):
    phicrit = 180/np.pi * np.sqrt(2 * delta)
    qcrit = 2*k * np.sin(np.pi/180 * phicrit)
    return (phicrit, qcrit)


def parratt(z):
    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(ai)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(ai)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(ai)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * kz1 * kz2 * s1**2)
    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * kz2 * kz3 * s2**2)
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2
    # Strecke vor Beginn der Oszillationen auf 1 setzen
    for i in np.arange(np.size(par)):
        if (i <= 296):  # 296 manuell angepasst
            par[i] = 1
        else:
            pass
    return par


def parratt2(z):
    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(ai)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(ai)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(ai)**2))
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2
    # Strecke vor Beginn der Oszillationen auf 1 setzen
    for i in np.arange(np.size(par)):
        if (i <= 296):  # 296 manuell angepasst
            par[i] = 1
            r13[i] = 1
        else:
            pass
    return par, r13


plotten(x, y)
