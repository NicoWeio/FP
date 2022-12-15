from common import load_lara, load_spe, plot_energyspectrum


def main(channel_to_E):
    x, N, T = load_spe("data/2022-11-28/3_Ba.Spe")
    E = channel_to_E(x)

    lit_energies_all_Ba, intensities_all_Ba = load_lara("data/emissions/Ba-133.lara.txt")
    lit_energies_all_Sb, intensities_all_Sb = load_lara("data/emissions/Sb-125.lara.txt")

    # ▒ Plotte Spektrum
    plot_energyspectrum(
        E, N,
        lit_energies_dict={
            "Ba-133": (lit_energies_all_Ba, intensities_all_Ba),
            "Sb-125": (lit_energies_all_Sb, intensities_all_Sb),
        },
        path="build/plt/spektrum_133-Ba.pdf",
        # smooth_over=20,
        # stack_lit_energies=True,
    )
