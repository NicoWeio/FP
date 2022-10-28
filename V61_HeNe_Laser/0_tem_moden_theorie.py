from matplotlib import pyplot as plt
import numpy as np
from scipy.special import hermite


def draw_tem(ax, l, m):
    def E(x, y):
        # https://en.wikipedia.org/wiki/Gaussian_beam#Hermite-Gaussian_modes
        # ohne Phaseninformation und Skalierung
        return hermite(n=l)(x) * hermite(n=m)(y) * np.exp(-(x**2 + y**2))

    size = 2
    x = np.linspace(-size, size, 1000)
    y = np.linspace(-size, size, 1000)
    X, Y = np.meshgrid(x, y)
    I = abs(E(X, Y))**2
    ax.imshow(I, cmap='gray')
    ax.set_xticks([])
    ax.set_yticks([])


def draw_tem_grid(gridsize):
    # Sollte aussehen wie https://en.wikipedia.org/wiki/Gaussian_beam#/media/File:Hermite-gaussian.png
    fig, axes = plt.subplots(gridsize, gridsize)

    for i in range(gridsize):
        for j in range(gridsize):
            draw_tem(axes[i, j], i, j)

    plt.show()


TEMs = [
    (0, 0),
    (1, 0),
]


for l, m in TEMs:
    # Bilder ohne Rand zu exportieren ist nicht gerade trivial: https://stackoverflow.com/a/8218887
    fig = plt.figure(frameon=False)
    fig.set_size_inches(1, 1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    draw_tem(fig.gca(), l, m)
    fig.savefig(f'build/plt/tem_{l}{m}.png', dpi=1000)
    # plt.show()
