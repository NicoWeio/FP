import numpy as np
import pandas as pd

dat = pd.read_csv('dat2/Unb.csv')
# dat.rename(str.strip)
# dat.rename_axis()



def A_row_from_indices(indices):
    print('→', indices)
    indices = indices.strip() # Zwischenlösung
    baserow = np.zeros(9)  # 9 = Anz. Würfel
    if '/' in indices:
        assert '|' not in indices
        # → diagonal
        # baserow[np.array(map(int,indices.split('/')))] = np.sqrt(2)
        for i in map(int, indices.split('/')):
            baserow[i-1] = np.sqrt(2)
    elif '|' in indices:
        assert '/' not in indices
        # → parallel
        # baserow[np.array(map(int,indices.split('|')))] = 1
        for i in map(int, indices.split('|')):
            baserow[i-1] = 1
    return baserow

def A_from_indices(all_indices):
    return np.row_stack(list(map(A_row_from_indices, all_indices)))

A = A_from_indices(dat['indices'])
print(A.round(1))
