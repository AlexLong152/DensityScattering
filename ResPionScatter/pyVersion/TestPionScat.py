# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
from PionScatLib import Tlab2Wcm, getGH, getcsGH


def main():
    Tlab = 100
    Wcm = Tlab2Wcm(Tlab)
    sqrtS = Wcm
    print("Wcm=", Wcm, "  Tlab=", Tlab)
    print("g and h values are weird here, and do not refer to Pions and Nuceli labels")
    thetas = np.arange(0, 181, 30)
    topStr = (
        "  Theta     DSG             f                  g/sin(theta)            f ratio"
    )
    topStr += "             g ratio         abs(f ratio)      abs(g ratio)"
    reactions = ["pi+p", "pi-p"]
    print("f and g results are in Fm, DSG(mb)=10*(f**2+g**2)")
    fsPlusp = np.array(
        [
            0.600989 + 0.353415j,
            0.409825 + 0.317779j,
            0.183982 + 0.196838j,
            -0.196938 + 0.029738j,
            -0.588941 - 0.137558j,
            -0.881934 - 0.260161j,
            -0.990733 - 0.305077j,
        ]
    )

    gsPlusp = np.array(
        [
            0.454827 + 0.159417j,
            0.460794 + 0.165636j,
            0.459258 + 0.165683j,
            0.463200 + 0.165800j,
            0.470091 + 0.165982j,
            0.476839 + 0.166165j,
            0.479505 + 0.166244j,
        ]
    )

    fsMinusp = np.array(
        [
            0.308959 + 0.155428j,
            0.352815 + 0.133858j,
            0.192784 + 0.090937j,
            0.061296 + 0.030082j,
            -0.049513 - 0.031133j,
            -0.122461 - 0.076067j,
            -0.147733 - 0.092535j,
        ]
    )
    gsMinusp = np.array(
        [
            0.154558 + 0.061401j,
            0.148789 + 0.059283j,
            0.155114 + 0.059225j,
            0.158785 + 0.059167j,
            0.161073 + 0.059127j,
            0.162326 + 0.059107j,
            0.162921 + 0.059101j,
        ]
    )

    gsDict = {}
    fsDict = {}
    ReactNum = {"pi+p": 1, "pi-p": 2}

    isospinVal = {"pi+p": 1, "pi-p": 1}
    piCharge = {"pi+p": 1, "pi-p": -1}
    gsDict["pi+p"] = gsPlusp
    gsDict["pi-p"] = gsMinusp

    fsDict["pi+p"] = fsPlusp
    fsDict["pi-p"] = fsMinusp
    for reaction in reactions:
        print("")
        print(133 * "#")
        print(133 * "#")
        print("\n")
        print(
            "Reaction "
            + str(reaction)
            + ", in Ron Workman's code this is reaction #"
            + str(ReactNum[reaction])
        )
        fs = fsDict[reaction]
        gs = gsDict[reaction]
        print(topStr)

        for i, theta in enumerate(thetas):
            x = np.cos(theta * np.pi / 180)
            f, g = getGH(
                sqrtS, x, isospin=isospinVal[reaction], piCharge=piCharge[reaction]
            )
            DSG = abs(f) ** 2 + abs(g * np.sin(theta * np.pi / 180)) ** 2
            DSG = DSG * 10
            # DSG = getcsGH(sqrtS, x, isospinVal[reaction], piCharge[reaction])
            fRatio = fs[i] / f
            gRatio = gs[i] / g
            fRatioAbs = abs(fRatio)
            gRatioAbs = abs(gRatio)
            print(
                f"{theta:7.2f}  {DSG:9.6f}   "
                f"{f.real:9.5f}{f.imag:+9.5f}j  "
                f"{g.real:9.5f}{g.imag:+9.5f}j  -- "
                f"{fRatio.real:8.4f}{fRatio.imag:+8.4f}j   "
                f"{gRatio.real:8.4f}{gRatio.imag:+8.4f}j     "
                f"{fRatioAbs:8.4}          "
                f"{gRatioAbs:8.4} "
            )

        print("\n\nFortran result from Ron Workman, ratios are to these values")

        print("  Theta      DSG             f                 g/sin(theta)")
        for i, theta in enumerate(thetas):
            DSG = abs(fs[i]) ** 2 + abs(gs[i] * np.sin(theta * np.pi / 180)) ** 2
            DSG = DSG * 10
            f = fs[i]
            g = gs[i]
            print(
                f"{theta:7.2f}  {DSG:9.6f}   "
                f"{f.real:9.5f}{f.imag:+9.5f}j  "
                f"{g.real:9.5f}{g.imag:+9.5f}j "
            )


if __name__ == "__main__":
    main()
