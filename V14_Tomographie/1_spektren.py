import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from numpy.linalg import inv
import pandas as pd
import pint
import rich
from rich.console import Console
import tools
ureg = pint.UnitRegistry()
console = Console()

path = '2022-06-20_Messdaten/' 'Würfel1_456.Spe'

# N = np.genfromtxt(path, dtype='int64', skip_header=1)

with open(path, 'r') as f:
    lines = f.readlines()

# discard all lines that don't start with a space
lines = [line for line in lines if line[0] == ' ']

# convert the lines to a numpy array (int64)
N = np.array(list(map(int, [line.split()[0] for line in lines])))

# cut off after the last non-zero element
# N = N[:N.nonzero()[0][-1]]

# cut off after 99% of the events
N_cumsum = np.cumsum(N)
N_cumsum_max = N_cumsum[-1] * 0.99
N_cumsum_max_index = np.where(N_cumsum > N_cumsum_max)[0][0]
N = N[:N_cumsum_max_index]


# print the last 15 lines
print(N[-15:])

plt.bar(np.arange(len(N)), N)
# plt.yscale('log')
plt.ylabel('N')
plt.show()
