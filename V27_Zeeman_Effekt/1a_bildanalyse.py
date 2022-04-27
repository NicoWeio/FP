import itertools
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
from PIL import Image, ImageEnhance

timestampImage = Image.open('img/timestamp.png')


def preprocess_image(path, rotate_deg=0):
    """Entferne den Zeitstempel und rotiere das Bild, sodass die Linien vertikal stehen."""
    # Bild einlesen
    colorImage = Image.open(path)

    # Zeitstempel entfernen
    # Variante 1:
    # colorImage.paste(timestampImage, (0, 0), timestampImage)
    # Variante 2:
    colorImage = Image.alpha_composite(colorImage.convert('RGBA'), timestampImage).convert('RGB')

    # rotieren
    rotated = colorImage.rotate(rotate_deg)

    return np.array(rotated)


def display_image(img):
    """Passt das Bild zur Darstellung an."""
    pil_img_in = Image.fromarray(img)
    pil_img_out = ImageEnhance.Brightness(pil_img_in).enhance(2.5)
    # pil_img_out = ImageEnhance.Contrast(pil_img_in).enhance(2.5)
    return np.array(pil_img_out)


def get_peaks(img, min_distance=100, min_height=0.4, prominence=0.2, name=None, show=False):
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

    # █ Plot
    img_height = img.shape[0]
    # plt.figure()
    fig, ax = plt.subplots()
    ax.imshow(display_image(img), extent=[0, 4000, 0, 2248])
    # origin='lower' dreht Achsen und Bild
    x = np.array(range(len(sums)))
    displaysums = sums * img_height
    ax.plot(x, displaysums, color='g', alpha=0.8)
    ax.plot(peaks, displaysums[peaks], 'x', alpha=0.5)
    if min_height:
        ax.axhline(min_height * img_height, color='gray')
    for peak in peaks:
        ax.axvline(x=peak, color='r', alpha=0.5)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f"build/plt/{name}.pdf", bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()

    return peaks
