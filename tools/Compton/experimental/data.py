# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    omega2012 = 60  # in MeV
    data2012 = np.array(
        [
            [40.0, 196.0, 2.0 * np.sqrt(34)],
            [55.0, 176.0, 10.0],
            [75.0, 125.0, np.sqrt(65)],
            [90.0, 135.0, np.sqrt(65)],
            [110.0, 138.0, np.sqrt(89)],
            [125.0, 177.0, 3.0 * np.sqrt(13)],
            [145.0, 178.0, 3.0 * np.sqrt(13)],
            [159.0, 193.0, 10.0],
        ],
        dtype=float,
    ).T

    omega2014 = 86  # in MeV
    data2014 = np.array(
        [
            [40.0, 203.0, 2.0 * np.sqrt(106)],
            [55.0, 147.0, np.sqrt(170)],
            [75.0, 140.0, np.sqrt(193)],
            [90.0, 159.0, np.sqrt(233)],
            [110.0, 146.0, np.sqrt(449)],
            [145.0, 167.0, 2.0 * np.sqrt(97)],
            [159.0, 172.0, 15.0],
        ],
        dtype=float,
    ).T
    omegas = [omega2012, omega2014]
    # labels = [f"{omega2012}, 2012", f"{omega2014}, 2014 data"]
    data = [data2012, data2014]
    fig, axs = plt.subplots(2, 1, figsize=(10, 10), layout="tight")
    for i, dataset in enumerate(data):
        axs[i].errorbar(
            dataset[0],
            dataset[1],
            yerr=dataset[2],
            # label=labels[i],
            linestyle=None,
            fmt="o",
            capsize=5,
        )
        # axs[i].legend()

        title = (
            r"Compton ${}^{6} \mathrm{Li}$ data at $"
            + str(omegas[i])
            + r"\mathrm{MeV}$"
        )
        axs[i].set_title(title)
        axs[i].set_ylabel(r"${d\sigma}/{d\Omega}\; [\mathrm{\mu barn} ]$", fontsize=12)
        axs[i].set_xlabel(r"$\theta$", fontsize=12)
    plt.show()


if __name__ == "__main__":
    main()
