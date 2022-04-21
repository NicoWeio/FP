import itertools
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
from PIL import Image, ImageEnhance

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


def display_image(img):
    """Passt das Bild zur Darstellung an."""
    pil_img_in = Image.fromarray(img)
    pil_img_out = ImageEnhance.Contrast(pil_img_in).enhance(2.5)
    return np.array(pil_img_out)


def get_peaks(img, min_distance=1, min_height=0.0, prominence=0, show=False):
    """
    min_height: Diesen Anteil vom Maximum muss ein Peak haben, damit er als solcher gewertet wird.
    """
    # summiere Farbkanäle
    sums1 = img.sum(axis=2)
    # summiere entlang der vertikalen Achse
    sums = sums1.sum(axis=0)
    # normalisiere zwischen 0 und 1
    # sums /= sums.max()
    sums = sums / max(sums)

    peaks, _ = sp.signal.find_peaks(
        sums,
        distance=min_distance,
        height=min_height,
        prominence=prominence,
    )

    peak_dists = np.diff(peaks)
    print(peak_dists)
    print(peak_dists.mean())

    img_height = sums1.shape[0]
    # plt.figure()
    fig, ax = plt.subplots()
    ax.imshow(display_image(img), extent=[0, 4000, 0, 2248])
    # origin='lower' dreht Achsen und Bild
    x = np.array(range(len(sums)))
    displaysums = sums * img_height
    ax.plot(x, displaysums, color='g', alpha=0.25)
    ax.plot(peaks, displaysums[peaks], 'x', alpha=0.5)
    if min_height:
        ax.axhline(min_height * img_height, color='gray')
        # ax.axhline(1600, color='gray')
    for peak in peaks:
        ax.axvline(x=peak, color='r', alpha=0.5)
    plt.axis('off')
    plt.tight_layout()
    # plt.savefig("build/plt/bildanalyse_1.jpg", bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()

    return peaks
    # return peak_dists


def get_Δs(img, **kwargs):
    """Analyse ohne B-Feld/Aufspaltung"""
    peaks = get_peaks(img, **kwargs)

    Δs = np.diff(peaks).mean()
    # Δs /= 2  # !?
    return Δs


def get_δs(img, **kwargs):
    """Analyse mit B-Feld/Aufspaltung"""
    peaks = get_peaks(img, **kwargs)

    # group peaks2 into pairs of 2
    pairs = list(itertools.pairwise(peaks))[::2]

    diffs = [b - a for a, b in pairs]

    # ↓ Ähm, mache ich das wirklich nicht schon in pairs?
    # use only inner distances
    # print("pre", diffs)
    # diffs = diffs[::2]
    # print("post", diffs)

    δs = np.array(diffs).mean()
    return δs
