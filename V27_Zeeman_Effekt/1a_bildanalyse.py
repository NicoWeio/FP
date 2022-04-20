from matplotlib.image import imread
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal

img = imread('IMG_0021_rotated.jpg')

# print("img", img)
# sum colors
sums1 = img.sum(axis=2)
# print("sums1", sums1)
# sum along its vertical axis
sums = sums1.sum(axis=0)
# print(sums)

min_height = 0.5 * max(sums)
# min_height = 0
peaks, _ = sp.signal.find_peaks(sums, height=min_height)
# print(peaks)

peak_dists = np.diff(peaks)
print(peak_dists)
print(peak_dists.mean())

fig, ax = plt.subplots()
ax.imshow(img)
# ax.imshow(img, extent=[0, 4000, 0, 2248])
x = np.array(range(len(sums)))
displaysums = sums / max(sums) * sums1.shape[0]
ax.plot(x, displaysums)
ax.plot(peaks, displaysums[peaks], 'x')
plt.show()
