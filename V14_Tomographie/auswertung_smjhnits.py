import numpy as np
from scipy.stats import sem
from uncertainties import ufloat
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit
from numpy.linalg import inv

A = np.matrix([[0, 0, 0, 0, 0, np.sqrt(2), 0, np.sqrt(2), 0],
               [0, 0, np.sqrt(2), 0, np.sqrt(2), 0, np.sqrt(2), 0, 0],
               [0, np.sqrt(2), 0, np.sqrt(2), 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 1, 1, 1],
               [0, 0, 0, 1, 1, 1, 0, 0, 0],
               [1, 1, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, np.sqrt(2), 0, 0, 0, np.sqrt(2), 0],
               [np.sqrt(2), 0, 0, 0, np.sqrt(2), 0, 0, 0, np.sqrt(2)],
               [0, np.sqrt(2), 0, 0, 0, np.sqrt(2), 0, 0, 0],
               [1, 0, 0, 1, 0, 0, 1, 0, 0],
               [0, 1, 0, 0, 1, 0, 0, 1, 0],
               [0, 0, 1, 0, 0, 1, 0, 0, 1]])

A_t = np.transpose(A)

Alu_lit = 0.202
Blei_lit = 1.25

#   Aluminium

Alu_Werte = np.array([8206, 2268, 2135, 2336, 2837, 2186, 2138, 2204, 3057, 2504, 2228, 2446])
Alu_Zeiten = np.array([84, 30, 26, 26, 30, 24, 24, 29, 34, 26, 24, 23])
Alu_Fehler = np.sqrt(Alu_Werte)
Alu_Rate_2 = Alu_Werte/Alu_Zeiten
Alu_Rate_Fehler = Alu_Fehler/Alu_Zeiten

Alu_Real = np.array([ufloat(n, Alu_Fehler[i]) for i, n in enumerate(Alu_Werte)])
Alu_Rate = Alu_Real/Alu_Zeiten

#print(Alu_Zeiten, '\n', Alu_Werte, '\n', Alu_Fehler, '\n')
#print("I:", Alu_Rate, '\n')

#   Blei

Blei_Werte = np.array([1218, 1197, 1203, 1198, 1205, 1206, 1212, 1197, 1240, 1255, 1207, 1207])
Blei_Zeiten = np.array([110, 375, 203, 161, 183, 187, 99, 424, 182, 124, 192, 195])
Blei_Fehler = np.sqrt(Blei_Werte)
Blei_Rate_2 = Blei_Werte/Blei_Zeiten
Blei_Rate_Fehler = Blei_Fehler/Blei_Zeiten

Blei_Real = np.array([ufloat(n, Blei_Fehler[i]) for i, n in enumerate(Blei_Werte)])
Blei_Rate = Blei_Real/Blei_Zeiten

#print(Blei_Zeiten, '\n', Blei_Werte, '\n', Blei_Fehler, '\n')
#print("I:", Blei_Rate, '\n')

#  Unbekannt

Unb_Werte = np.array([12575, 12758, 12622, 12916, 12965, 13018, 12905, 12840, 12898, 12888, 12510, 12536])
Unb_Zeiten = np.array([822, 1503, 690, 648, 700, 683, 341, 1345, 1193, 125, 1407, 799])
Unb_Fehler = np.sqrt(Unb_Werte)
Unb_Rate_2 = Unb_Werte/Unb_Zeiten
Unb_Rate_Fehler = Unb_Fehler/Unb_Zeiten

Unb_Real = np.array([ufloat(n, Unb_Fehler[i]) for i, n in enumerate(Unb_Werte)])
Unb_Rate = Unb_Real/Unb_Zeiten

#print(Unb_Zeiten, '\n', Unb_Werte, '\n', Unb_Fehler, '\n')
#print("I:", Unb_Rate, '\n')

#   Luft

Luft_Werte = np.array([16115, 16042, 16115, 15341, 15341, 15341, 16115, 16042, 16115, 15341, 15341, 15341])
Luft_Zeiten = np.array([87, 87, 87, 83, 83, 83, 87, 87, 87, 83, 83, 83])
Luft_Fehler = np.sqrt(Luft_Werte)
Luft_Rate_2 = Luft_Werte/Luft_Zeiten
Luft_Rate_Fehler = Luft_Fehler/Luft_Zeiten

Luft_Real = np.array([ufloat(n, Luft_Fehler[i]) for i, n in enumerate(Luft_Werte)])
Luft_Rate = Luft_Real/Luft_Zeiten
#print("I0: ", Luft_Rate, '\n')

# Bestimmung der Varianzmatrizen

delta_y_Alu = np.array([np.sqrt((1/Luft_Rate_2[i])**2*n**2 + (1/Alu_Rate_2[i])**2*Alu_Rate_Fehler[i])
                       for i, n in enumerate(Luft_Rate_Fehler)])
delta_y_Pb = np.array([np.sqrt((1/Luft_Rate_2[i])**2*n**2 + (1/Blei_Rate_2[i])**2*Blei_Rate_Fehler[i])
                      for i, n in enumerate(Luft_Rate_Fehler)])
delta_y_Unb = np.array([np.sqrt((1/Luft_Rate_2[i])**2*n**2 + (1/Unb_Rate_2[i])**2*Unb_Rate_Fehler[i])
                       for i, n in enumerate(Luft_Rate_Fehler)])

#print("delta y Al: ", delta_y_Alu, '\n')
#print("delta y Pb: ", delta_y_Pb, '\n')
#print("delta y Unb: ", delta_y_Unb, '\n')

V_Al_N = np.diag(delta_y_Alu**2)
V_Pb_N = np.diag(delta_y_Pb**2)
V_Unb_N = np.diag(delta_y_Unb**2)

#  Bestimmung der Intensitäten

# Fehlerfortpflanzung

I_Alu = np.array([unp.log(n/Alu_Rate[i]) for i, n in enumerate(Luft_Rate)])
I_Blei = np.array([unp.log(n/Blei_Rate[i]) for i, n in enumerate(Luft_Rate)])
I_Unb = np.array([unp.log(n/Unb_Rate[i]) for i, n in enumerate(Luft_Rate)])

Alu_t = np.transpose(np.array([I_Alu]))
Blei_t = np.transpose(np.array([I_Blei]))
Unb_t = np.transpose(np.array([I_Unb]))

# Varianzmatrizen

I_Alu_2 = np.array([unp.log(n/Alu_Rate_2[i]) for i, n in enumerate(Luft_Rate_2)])
I_Blei_2 = np.array([unp.log(n/Blei_Rate_2[i]) for i, n in enumerate(Luft_Rate_2)])
I_Unb_2 = np.array([unp.log(n/Unb_Rate_2[i]) for i, n in enumerate(Luft_Rate_2)])

#print("y Alu: ", I_Alu_2, '\n')
#print("y Blei: ", I_Blei_2, '\n')
#print("y Unb: ", I_Unb_2, '\n')

Alu_t_2 = np.transpose(np.array([I_Alu_2]))
Blei_t_2 = np.transpose(np.array([I_Blei_2]))
Unb_t_2 = np.transpose(np.array([I_Unb_2]))

#Bestimmung der Koeffizienten

# Varianzmethode

V_Al_mu = inv(A_t*inv(V_Al_N)*A)
Al_mu_2 = inv(A_t*inv(V_Al_N)*A) * A_t*inv(V_Al_N)*Alu_t_2
mu_Alu_2 = np.array([ufloat(Al_mu_2[i], n) for i, n in enumerate(np.sqrt(np.diag(V_Al_mu)))])

V_Pb_mu = inv(A_t*inv(V_Pb_N)*A)
Pb_mu_2 = inv(A_t*inv(V_Pb_N)*A) * A_t*inv(V_Pb_N)*Blei_t_2
mu_Pb_2 = np.array([ufloat(Pb_mu_2[i], n) for i, n in enumerate(np.sqrt(np.diag(V_Pb_mu)))])

V_Unb_mu = inv(A_t*inv(V_Unb_N)*A)
Unb_mu_2 = inv(A_t*inv(V_Unb_N)*A) * A_t*inv(V_Unb_N)*Unb_t_2
mu_Unb_2 = np.array([ufloat(Unb_mu_2[i], n) for i, n in enumerate(np.sqrt(np.diag(V_Unb_mu)))])

print("Koeffizienten durch Varianz:", '\n')
#print("Alu: ", np.mean(mu_Alu_2),'\n')
#print("Blei: ", np.mean(mu_Pb_2),'\n')
print("Unbekannt: ", '\n', mu_Unb_2, '\n')

#print("Abweichung Alu:", np.mean(mu_Alu_2)-Alu_lit, '\n')
#print("Abweichung Blei:", np.mean(mu_Pb_2)-Blei_lit, '\n')

print("Abweichungen:", '\n')
print("Alu_lit:", '\n', mu_Unb_2-Alu_lit, '\n')
print("Alu_exp:", '\n', mu_Unb_2-np.mean(mu_Alu_2), '\n')
print("Blei_lit:", '\n', mu_Unb_2-Blei_lit, '\n')
print("Blei_exp:", '\n', mu_Unb_2-np.mean(mu_Pb_2), '\n')

print("----------------------------------------------------------------------", '\n')

# Durch Fehlerrechnung

mu_Alu = (inv(A_t*A)*A_t)*Alu_t
mu_Blei = (inv(A_t*A)*A_t)*Blei_t
mu_Unb = (inv(A_t*A)*A_t)*Unb_t

#print("Koeffizienten durch Fehlerrechnung:", '\n')
#print("Alu: ", np.mean(mu_Alu),'\n')
#print("Blei: ", np.mean(mu_Blei),'\n')
#print("Unbekannt: ", mu_Unb, '\n')
