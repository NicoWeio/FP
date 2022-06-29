import matplotlib.pyplot as plt
import numpy as np

name = 'Würfel1_456'
path = f'2022-06-20_Messdaten/{name}.Spe'

# Die Dateien beinhalten zusätzliche Zeilen an Anfang und Ende.
# Wir kümmern uns selbst darum…
with open(path, 'r') as f:
    lines = f.readlines()

# discard all lines that don't start with a space; then strip it away
lines = [line.strip() for line in lines if line.startswith(' ')]

# convert the lines to a NumPy array
N = np.array(list(map(int, lines)))

# cut off after the last non-zero element
# N = N[:N.nonzero()[0][-1]]

# cut off after 99% of the events
N_cumsum = np.cumsum(N)
N_cumsum_max = N_cumsum[-1] * 0.99
N_cumsum_max_index = np.where(N_cumsum > N_cumsum_max)[0][0]
N = N[:N_cumsum_max_index]


plt.bar(np.arange(len(N)), N, zorder=5)
plt.axvline(21, color='C1', label='Rückstrahlpeak')
plt.axvline(44, color='C2', label='Compton-Kante')
plt.axvline(np.argmax(N), color='C3', label='Photopeak')
# plt.yscale('log')
plt.xlabel('Channel')
plt.ylabel('$N$')
plt.legend()
plt.savefig('build/plt/spektrum.pdf')
# plt.show()
