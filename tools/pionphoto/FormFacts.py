# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import parseThreshOneBody as pt
from CrossSection import meanAndSpread, TwoBod_AveAndSpread


def main():
    # onebod()
    twobod()


def twobod():
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/"
    means, spread = TwoBod_AveAndSpread(folder, verbose=False)
    means *= 2
    spread *= 2
    print("means=", means)
    print("spread=", spread)


def onebod():
    f1 = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/132MeV/onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=fac31e8f17bc2b1a6d84411b408e05cca71aae63c04be8f6fbf3ba311acb45a3.v2.0.dat"
    f2 = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/132MeV/onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=3bd9b5d55f401de84f369724cdeb02eb176fda6928c442e1611c2ba14efdbc8b.v2.0.dat"
    files = [f1, f2]
    val1 = []
    val2 = []
    for f in files:
        d = pt.He3Parse(f)
        val1.append(d[1])
        val2.append(d[3])
        # for k, v in d.items():
        #     print(f"{k}: {v}")
    val1 = np.vstack((val1[0], val1[1])).T
    val2 = np.vstack((val2[0], val2[1])).T
    Transverse_minus = meanAndSpread(val1[0])
    Transverse_plus = meanAndSpread(val1[1])

    Long_minus = meanAndSpread(val2[0])
    Long_plus = meanAndSpread(val2[1])

    results = {
        "F_T^{S-V}": Transverse_minus,
        "F_T^{S+V}": Transverse_plus,
        "F_L^{S-V}": Long_minus,
        "F_L^{S+V}": Long_plus,
    }
    for label, (mean, spread) in results.items():
        print(f"{label} = {mean:.3f} +/- {spread:.3f}")


if __name__ == "__main__":
    main()
