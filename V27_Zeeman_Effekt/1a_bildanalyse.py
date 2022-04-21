import itertools
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
from PIL import Image

timestampImage = Image.open('img/timestamp.png')


def preprocess_image(path, rotate_deg=0):
    """Entferne den Zeitstempel und rotiere das Bild, sodass die Linien vertikal stehen."""
    colorImage = Image.open(path)
    # Zeitstempel entfernen
    colorImage.paste(timestampImage, (0, 0), timestampImage)
    # rotieren
    rotated = colorImage.rotate(rotate_deg)

    # rotated.show()
    return np.array(rotated)


def get_peaks(img, min_height=0.5, show=False):
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
    ax.plot(x, displaysums, alpha=0.25)
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


def get_Δs(img):
    """Analyse ohne B-Feld/Aufspaltung"""
    peaks = get_peaks(img, show=False)

    Δs = np.diff(peaks).mean()
    # Δs /= 2  # !?
    return Δs


def get_δs(img):
    """Analyse mit B-Feld/Aufspaltung"""
    peaks = get_peaks(img, min_height=0.4, show=False)

    # group peaks2 into pairs of 2
    pairs = list(itertools.pairwise(peaks))[::2]
    diffs = [b - a for a, b in pairs]
    δs = np.array(diffs).mean()
    return δs
