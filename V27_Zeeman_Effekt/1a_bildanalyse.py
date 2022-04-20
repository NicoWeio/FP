from matplotlib.image import imread
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal

img = imread('IMG_0021_rotated.jpg')

# display the image
imgplot = plt.imshow(img)
# plt.show()

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

plt.figure()
x = np.array(range(len(sums)))
plt.plot(x, sums)
plt.plot(peaks, sums[peaks], 'x')
plt.show()
