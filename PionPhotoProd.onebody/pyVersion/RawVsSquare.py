# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import PionPhotoLib as ppl
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    poleData = ppl.SaidPoles()
    xs = np.arange(-1, 1, 0.1)
    # ys=np.zeros(len(xs))
    nucs = ["pp0", "nn0"]
    sqrtS = 1300  # Mandelstam S
    fig, axs = plt.subplots(len(nucs), layout="tight")
    # for i, x in enumerate(xs):
    S = sqrtS**2
    for i, nuc in enumerate(nucs):
        RawYs = np.array([ppl.calcCCFromRaw(S, x, nuc, poleData) for x in xs])
        SquareYs = np.array([ppl.calcCrossSection(S, x, nuc, poleData) for x in xs])
        axs[i].plot(xs, RawYs, label="Raw")
        axs[i].plot(xs, SquareYs, label="Square")
        print("RawYs-SquareYs=\n", RawYs - SquareYs)
    axs[1].legend()
    plt.show()


if __name__ == "__main__":
    main()
