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
    omegaLab = 160
    thetas = np.arange(0, 181, 10)
    nucs = ["pp0", "nn0"]
    poleData = ppl.SaidPoles()

    fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # Values of "ccs" come from experiment
    for j, nuc in enumerate(nucs):
        if nuc == "pp0":
            mNucl = 938.272

            ccs = np.array(
                [
                    0.1007e00,
                    0.1027e00,
                    0.1084e00,
                    0.1175e00,
                    0.1293e00,
                    0.1432e00,
                    0.1581e00,
                    0.1729e00,
                    0.1866e00,
                    0.1982e00,
                    0.2070e00,
                    0.2126e00,
                    0.2151e00,
                    0.2151e00,
                    0.2133e00,
                    0.2106e00,
                    0.2080e00,
                    0.2061e00,
                    0.2054e00,
                ]
            )
        else:
            mNucl = 939.565
            assert nuc == "nn0"

            ccs = np.array(
                [
                    0.2351e00,
                    0.2362e00,
                    0.2391e00,
                    0.2430e00,
                    0.2464e00,
                    0.2477e00,
                    0.2456e00,
                    0.2387e00,
                    0.2265e00,
                    0.2091e00,
                    0.1871e00,
                    0.1618e00,
                    0.1348e00,
                    0.1078e00,
                    0.8286e-01,
                    0.6152e-01,
                    0.4522e-01,
                    0.3501e-01,
                    0.3154e-01,
                ]
            )

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
            ccsCalc.append(crossSec)
            ccExp = ccs[i]

            diff = abs(ccExp - crossSec)
            pct_err = 100.0 * diff / ccExp
            print(
                "theta={:4d}, cross sec={:8.6f} ---  Experimental Val={:8.6f}, Diff={:8.6f} --> {:6.3f}% Error".format(
                    theta, crossSec, ccExp, diff, pct_err
                )
            )
        ccCalc = np.array(ccsCalc)
        titles = {0: "Proton, Neutral Pion", 1: "Neutron, Neutral Pion"}
        axs[j].set_title(titles[j])
        axs[j].plot(thetas, ccs, c="b", label="Experiment")
        axs[j].plot(thetas, ccCalc, c="r", label="Calculated")
        # plt.legend()
        axs[j].legend()
        if j == 0:
            print("\n" + 100 * "#")
    fig.supxlabel(r"$\theta$")
    fig.supylabel(r"$\mathrm{d} \sigma/\mathrm{d}\Omega \;[\mu \mathrm{bn}]$")
    plt.show()


def resultOnly():
    thetas = np.arange(0, 181, 10)
    # ccs1 = np.array(
    #     [
    #         0.1007e00,
    #         0.1027e00,
    #         0.1084e00,
    #         0.1175e00,
    #         0.1293e00,
    #         0.1432e00,
    #         0.1581e00,
    #         0.1729e00,
    #         0.1866e00,
    #         0.1982e00,
    #         0.2070e00,
    #         0.2126e00,
    #         0.2151e00,
    #         0.2151e00,
    #         0.2133e00,
    #         0.2106e00,
    #         0.2080e00,
    #         0.2061e00,
    #         0.2054e00,
    #     ]
    # )
    #
    # ccs2 = np.array(
    #     [
    #         0.2351e00,
    #         0.2362e00,
    #         0.2391e00,
    #         0.2430e00,
    #         0.2464e00,
    #         0.2477e00,
    #         0.2456e00,
    #         0.2387e00,
    #         0.2265e00,
    #         0.2091e00,
    #         0.1871e00,
    #         0.1618e00,
    #         0.1348e00,
    #         0.1078e00,
    #         0.8286e-01,
    #         0.6152e-01,
    #         0.4522e-01,
    #         0.3501e-01,
    #         0.3154e-01,
    #     ]
    # )
    names = ["pp0", "nn0"]

    omegaLab = 160
    poleData = ppl.SaidPoles()
    print("\n")
    print("Now printing results we have calculated from poles directly")
    for nuc in names:
        print(nuc)
        print("Theta                Cross Section ")
        if nuc == "pp0":
            mNucl = 938.272
        else:
            mNucl = 939.565
        sqrtS = np.sqrt(2 * omegaLab * mNucl + mNucl * mNucl)
        for i, theta in enumerate(thetas):
            x = np.cos(theta * np.pi / 180)
            crossSec = ppl.calcCrossSection(sqrtS**2, x, nuc, poleData)
            print(f"{theta:3d}    --------    {crossSec:6.6f} [\u03bc barns]")
        print("\n" + 50 * "#")


if __name__ == "__main__":
    main()
    resultOnly()
