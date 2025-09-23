# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import CrossSection as cc
from matplotlib import pyplot as plt
from matplotlib import rcParams
import os

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


energy = 86
Odeltaonebod = "Odelta3"
Odeltatwobod = "Odelta4"
saveFolder = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/results/FinalCrossSections-And-Uncertainties/"
fileName = f"Compton-6Li-{energy}MeV-1bod={Odeltaonebod}-2bod={Odeltatwobod}-output.txt"
outfile = saveFolder + fileName
if os.path.exists(outfile):
    os.remove(outfile)


def main():
    Ntot = 14

    twobody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/2bod/{energy}MeV/"
    onebody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/1bod/{energy}MeV/"

    lambdaSRGs = np.array([1.880, 2.236])
    lambdaCuts = np.array([450, 500])
    omegaHs = np.array([16, 18])
    # omegaHs = np.array([18])
    thetas = np.array([1, 40, 55, 75, 90, 110, 125, 145, 159, 180])
    out = np.zeros((len(thetas), 2))

    counts = []
    for i, theta in enumerate(thetas):
        values = []
        numVal = 0
        for lambdaSRG in lambdaSRGs:
            for lambdaCut in lambdaCuts:
                for omegaH in omegaHs:
                    try:
                        crossSec = cc.ccForDict(
                            onebody_dir,
                            twobody_dir,
                            Odeltaonebod=Odeltaonebod,
                            Odeltatwobod=Odeltatwobod,
                            returnFull=True,
                            energy=energy,
                            angle=theta,
                            lambdaCut=lambdaCut,
                            lambdaSRG=lambdaSRG,
                            Ntotmax=Ntot,
                            omegaH=omegaH,
                        )["cc"]
                        values.append(crossSec)
                        numVal += 1

                    except ValueError:
                        pass
        counts.append(numVal)
        out[i] = np.mean(values), maxDev(values)
    # loop over mean values in a for loop and sort them into a new array

    plt.figure(figsize=(8, 6))
    plt.errorbar(
        thetas,
        out[:, 0],
        yerr=out[:, 1],
        linestyle="none",
        fmt="o",
        markersize=4,
        capsize=5,  # length of horizontal bars at the top/bottom
        capthick=1.5,  # thickness of those bars
    )
    plt.title(
        r"${}^6\mathrm{Li}$ Compton Scattering at $"
        + str(energy)
        + r"\, \mathrm{MeV}$",
        fontsize=15,
    )

    plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\Omega\; [\mathrm{nb}]$", fontsize=14)
    plt.xlabel(r"$\theta$ (Degrees) ", fontsize=14)

    wprint(f"{'theta':>6} {'Cross Section [nb]':>20} {'Uncertainty':>12} {'Points':>8}")
    for i, theta in enumerate(thetas):
        wprint(f"{theta:6.0f} {out[i, 0]:20.6f} {out[i, 1]:12.6f} {counts[i]:8d}")
    wprint("Polarisability values - deviations from Fortran code")
    wprint("cc.alphap=", cc.alphap)
    wprint("cc.alphan=", cc.alphan)
    wprint("cc.betap=", cc.betap)
    wprint("cc.betan=", cc.betan)
    wprint("Odeltaonebod=", Odeltaonebod)
    wprint("Odeltatwobod=", Odeltatwobod)
    figureSave = saveFolder + "plot" + fileName.replace(".txt", ".pdf")
    plt.tight_layout()
    plt.savefig(figureSave)
    plt.show()


def maxDev(ar):
    minVal = np.min(ar)
    maxVal = np.max(ar)
    return (maxVal - minVal) / 2


def wprint(*args, sep=" ", end="\n", file=None, flush=False):
    """
    Works like print(), but also writes to a global file `outfile`.
    """
    # First, print to the original destination (usually stdout)
    print(*args, sep=sep, end=end, file=file, flush=flush)

    # Then, append the same text to the outfile
    with open(outfile, "a", encoding="utf-8") as f:
        print(*args, sep=sep, end=end, file=f)


if __name__ == "__main__":
    main()
    print("outfile=\n", outfile)
