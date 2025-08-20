# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import CrossSection as cc

# from os import listdir
# from os.path import isfile, join
import sys

sys.path.insert(1, "..")
import readDensity as rd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    base = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/"

    energy = 60
    # angle = 180
    lambdaSRG = 1.880
    Ntotmax = 14
    lambdaCuts = [450, 500]
    lambdaSRGs = [1.880, 2.236]
    omegaH = 18
    plotting = "theta"

    thetas = np.array([40, 55, 75, 90, 110, 125, 145, 159, 180, 1])
    xs = thetas
    ccs = []
    ccs0 = {}
    xsPlot = []
    dString = ""
    MeVtmp = 0

    fig, ax_scatter = plt.subplots(figsize=(10, 7))
    markers = np.array(["1", "^", "+", "x"])
    colors = np.array(["b", "g", "r", "k"])

    resultSave = base + "results/"
    onebody_dir = base + r"1bod/ENERGY!!!MeV/"
    twobody_dir = base + r"2bod/ENERGY!!!MeV/"
    onebody_dir = onebody_dir.replace("ENERGY!!!", str(energy))
    twobody_dir = twobody_dir.replace("ENERGY!!!", str(energy))
    # print("onebody_dir=", onebody_dir)
    # print("twobody_dir=", twobody_dir)
    outDict = {}
    #
    for lambdaSRG in lambdaSRGs:
        for lambdaCut in lambdaCuts:
            outDict[(lambdaSRG, lambdaCut)] = []
    for x in xs:
        # kwargsfind[plotting] = x

        # print("kwargsfind=\n", kwargsfind)
        i = 0
        for lambdaSRG in lambdaSRGs:
            for lambdaCut in lambdaCuts:
                kwargsfind = {
                    "energy": energy,
                    "angle": x,
                    "lambdaCut": lambdaCut,
                    "lambdaSRG": lambdaSRG,
                    "Ntotmax": Ntotmax,
                    "omegaH": omegaH,
                }
                onebod, twobod, kwargs = cc.getMatchingFiles(
                    onebody_dir, twobody_dir, **kwargsfind
                )

                if onebod is not None and twobod is not None:
                    yVal = cc.crossSection(onebody_dir + onebod, twobody_dir + twobod)[
                        "cc"
                    ]
                    ccs.append(yVal)
                    if x == 40:
                        ccs0[(lambdaSRG, lambdaCut)] = yVal
                    xsPlot.append(x)
                    if x == xs[0]:
                        ax_scatter.scatter(
                            x,
                            ccs0[(lambdaSRG, lambdaCut)] / yVal,
                            marker=markers[i],
                            c=colors[i],
                            label=rf"$\Lambda_{{SRG}}={lambdaSRG}, \Lambda_{{NN}}={lambdaCut}$ at 40 deg/value at $\theta$",
                        )
                    else:
                        ax_scatter.scatter(
                            x, ccs[0] / yVal, marker=markers[i], c=colors[i]
                        )
                i += 1
                # print("xsPlot=", xsPlot)
                # print("ccs=", ccs)
                # else:
                #     print(f"For x={x}:\n   onebod={onebod}\n   twobod={twobod}\n\n")

        dString = cc.dictString(kwargs)
        MeVtmp = kwargs["energy"]
    # make a figure with 1 row, 2 columns

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
    )
    data = data2012
    ccAt40 = data[0][1]
    errAt40 = data[0][2]
    for i, dataset in enumerate(data):
        upperBound = ccAt40 / (dataset[1] + dataset[2])
        yValue = ccAt40 / dataset[1]

        yerr = np.array(
            [
                [((ccAt40 + errAt40) / dataset[1])],  # lower
                [((ccAt40 - errAt40) / dataset[1])],  # upper
            ]
        )
        yerr = abs(yerr - yValue)
        if i == 0:
            ax_scatter.errorbar(
                dataset[0],
                yValue,
                yerr=yerr,
                # label=labels[i],
                linestyle=None,
                fmt="o",
                capsize=5,
                c="C0",
                label="2012 Experiment (Gerry) at 40 deg/value at $\\theta$",
            )
        else:
            ax_scatter.errorbar(
                dataset[0],
                yValue,
                yerr=yerr,
                linestyle=None,
                c="C0",
                fmt="o",
                capsize=5,
            )

    ax_scatter.set_ylabel(r"$\mathrm{d} \sigma /\mathrm{d} \Omega$ ", fontsize=14)
    ax_scatter.set_xlabel("$\\theta$", fontsize=14)
    # ax_scatter.set_ylim([0, 1.2 * np.max(ccs)])
    ax_scatter.set_title(
        r"$\mathrm{d} \sigma /\mathrm{d} \Omega$ vs "
        + plotting
        + r" for Compton Scattering on ${}^6\mathrm{Li}$"
        + f" at ${MeVtmp} \\mathrm{{MeV}}$",
        fontsize=15,
    )

    plt.tight_layout()
    plt.legend()
    # saveString = "plot_" + plotting
    # for key, value in kwargsfind.items():
    #     if key != plotting:
    #         saveString += "-" + key + "=" + str(value)
    saveString = (
        "6Li"
        + "-"
        + str(energy)
        + "MeV-omegaH="
        + str(kwargsfind["omegaH"])
        + "-Ntotmax="
        + str(kwargsfind["Ntotmax"])
    )
    # if factor != 1.0:
    #     saveString = "factor=" + str(factor) + "-" + saveString
    saveString = "Ratio-" + saveString + ".pdf"
    saveString = resultSave + saveString

    plt.show()

    for lambdaSRG in lambdaSRGs:
        for lambdaCut in lambdaCuts:
            print("lambdaSRG=", lambdaSRG, "lambdaCut=", lambdaCut)
            print(40 * "-")
            outList = outDict[(lambdaSRG, lambdaCut)]
            print(f"{'theta':>7} | {'Cross Section':>14}")
            for val in outList:
                x, y = val
                print(f"{x:7.2f} , {y:14.6f}")
    # print("\n")
    # print("Gerry experimental Data")
    # print(f"{'theta':>7} | {'val':>10} | {'uncertainty':>12}")
    # print(40 * "-")
    for theta, val, uncer in data:
        print(f"{theta:7.2f} | {val:10.6f} | {uncer:12.6f}")
    yn = input("Save it?  [y/n]" + "\n" + saveString + "\n").lower()
    if yn == "y":
        fig.savefig(saveString, format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
