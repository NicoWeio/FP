import bildanalyse
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
from PIL import Image, ImageEnhance


img = bildanalyse.preprocess_image('Bilder/blau/5_8A_pi/IMG_0034.JPG', rotate_deg=-2.3)


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

    peaks, peak_props = sp.signal.find_peaks(
        sums,
        distance=1,
        height=min_height,
        prominence=prominence,
        width=(1, 50),
    )

    # █ Plot
    img_height = img.shape[0]
    # plt.figure()
    fig, ax = plt.subplots()
    ax.imshow(display_image(img), extent=[0, 4000, 0, 2248])
    # origin='lower' dreht Achsen und Bild
    x = np.array(range(len(sums)))
    displaysums = sums * img_height
    ax.plot(x, displaysums, color='g', alpha=0.8, label='Signal')
    ax.plot(peaks, displaysums[peaks], 'x', alpha=0.5, label='Maxima')

    # if min_height:
    #     ax.axhline(min_height * img_height, color='gray')

    peak_widths = sp.signal.peak_widths(sums, peaks, rel_height=0.4) # wichtiger Parameter rel_height, war 0.5!
    # plt.hlines(*peak_widths[1:], color='r', alpha=0.5)
    plt.hlines(peak_widths[1] * img_height, peak_widths[2], peak_widths[3], color='r', alpha=0.5)

    for pw in zip(*peak_widths):
        ax.axvline(x=pw[2], color='r', alpha=0.5)
        ax.axvline(x=pw[3], color='r', alpha=0.5)

    # for peak in peaks:
    #     ax.axvline(x=peak, color='r', alpha=0.5)

    plt.axis('off')
    plt.legend()
    plt.tight_layout()
    # plt.savefig(f"build/plt/{name}.pdf", bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()

    # return peaks
    widths = peak_widths[0]
    return widths.mean()

w = get_peaks(
    img,
    show=True,
)
print(w)
