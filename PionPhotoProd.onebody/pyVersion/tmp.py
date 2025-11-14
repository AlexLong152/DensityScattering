# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import PionPhotoLib as ppl


def main():
    omegaLab = 160
    #   thetas = np.arange(0, 181, 10)
    theta = 0
    thetas = np.array([0])
    # print(thetas)
    # print(len(thetas))
    nucs = ["pp0", "nn0"]
    poleData = ppl.SaidPoles()
    for nuc in nucs:
        if nuc == "pp0":
            mNucl = 938.272
        else:
            assert nuc == "nn0"
            mNucl = 939.565

        sqrtS = np.sqrt(2 * omegaLab * mNucl + mNucl * mNucl)

        print("nuc=", nuc)
        print("sqrtS=", sqrtS)
        print("")

        ccsCalc = []
        for i, theta in enumerate(thetas):
            S = sqrtS**2
            x = np.cos(theta * np.pi / 180)
            # crossSec = ppl.calcCrossSectionFromRaw(S, x, nuc, poleData)
            crossSec = ppl.calcCrossSection(S, x, nuc, poleData)
            print("crossSec=", crossSec)
            ccsCalc.append(crossSec)

        print("\n\n")


if __name__ == "__main__":
    main()
