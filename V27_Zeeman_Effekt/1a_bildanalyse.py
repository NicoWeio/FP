import itertools
from matplotlib.image import imread
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal


def analyze_img(img, min_height=0.5, show=False):
    """
    min_height: Diesen Anteil vom Maximum muss ein Peak haben, damit er als solcher gewertet wird.
    """
    # summiere Farbkanäle
    sums1 = img.sum(axis=2)
    # summiere entlang der vertikalen Achse
    sums = sums1.sum(axis=0)

    min_height = 0.5 * max(sums)
    peaks, _ = sp.signal.find_peaks(
        sums,
        height=min_height,
        distance=70,
    )

    peak_dists = np.diff(peaks)
    print(peak_dists)
    print(peak_dists.mean())

    # plt.figure()
    fig, ax = plt.subplots()
    # ax.imshow(img)
    ax.imshow(img, extent=[0, 4000, 0, 2248])
    # origin='lower' dreht Achsen und Bild
    x = np.array(range(len(sums)))
    displaysums = sums / max(sums) * sums1.shape[0]
    ax.plot(x, displaysums)
    # ax.plot(peaks, displaysums[peaks], 'x')
    for peak in peaks:
        ax.axvline(x=peak, color='r', alpha=0.25)
    plt.axis('off')
    plt.tight_layout()
    # plt.savefig("build/plt/bildanalyse_1.jpg", bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()

    return peaks
    # return peak_dists


img1 = imread('img/rot_0A.jpg')
img2 = imread('img/rot_8A.jpg')


peaks1 = analyze_img(img1, show=False)
peaks2 = analyze_img(img2, min_height=0.4, show=False)

Δs = np.diff(peaks1).mean()
Δs /= 2  # !?

# group peaks2 into pairs of 2

pairs = list(itertools.pairwise(peaks2))[::2]
diffs = [b - a for a, b in pairs]
δs = np.array(diffs).mean()
print(peaks2)
print(pairs)
print(diffs)
print(δs)

