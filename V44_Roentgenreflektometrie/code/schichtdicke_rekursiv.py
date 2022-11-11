import numpy as np

k_zj = k * np.sqrt(np.abs(n_j**2 - np.cos(α)**2))
r_jjplus1 = (k_zj - k_zjplus1) / (k_zj + k_zjplus1)
X_j = np.exp(-2j * k_zj * z_j) * r23
