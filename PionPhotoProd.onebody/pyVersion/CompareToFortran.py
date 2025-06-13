# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import PionPhotoLib as ppl


def main():
    sqrtS = 1086.555
    thetas = np.arange(0, 181, 10)
    print("sqrtS=", sqrtS)
    nuc = "pp0"
    poleData = ppl.SaidPoles()

    for theta in thetas:
        S = sqrtS**2
        x = np.cos(theta * np.pi / 180)
        cc = ppl.calcCrossSectionFromRaw(S, x, nuc, poleData)
        print(cc, "---", theta)


if __name__ == "__main__":
    main()
