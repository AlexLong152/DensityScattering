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


alphaDiff = -2.0
betaDiff = 0.0


def main():
    base = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/"

    energy = 86
    # angle = 180
    lambdaSRG = 1.880
    Ntotmax = 14
    lambdaCuts = [450, 500]
    lambdaSRGs = [1.880, 2.236]
    omegaH = 18
    plotting = "theta"

    thetas = np.array([1, 40, 55, 75, 90, 110, 125, 145, 159, 180])
    xs = thetas
    ccs = []
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
    outDict = {}

    for lambdaSRG in lambdaSRGs:
        for lambdaCut in lambdaCuts:
            outDict[(lambdaSRG, lambdaCut)] = []

    for x in xs:
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
                    yVal = cc.crossSection(
                        onebody_dir + onebod,
                        twobody_dir + twobod,
                        delta=[alphaDiff, betaDiff],
                    )["cc"]
                    ccs.append(yVal)
                    outDict[(lambdaSRG, lambdaCut)].append([float(x), yVal])
                    if x == xs[0]:
                        ax_scatter.scatter(
                            x,
                            yVal,
                            marker=markers[i],
                            c=colors[i],
                            label=rf"$\Lambda_{{SRG}}={lambdaSRG}, \Lambda_{{NN}}={lambdaCut}$",
                        )
                    else:
                        ax_scatter.scatter(x, yVal, marker=markers[i], c=colors[i])
                i += 1

        dString = cc.dictString(kwargs)
        MeVtmp = kwargs["energy"]
    # make a figure with 1 row, 2 columns

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
    )
    # factor = 1.0
    data = data2014
    for i, dataset in enumerate(data):
        if i == 0:
            ax_scatter.errorbar(
                dataset[0],
                dataset[1],  # * factor,
                yerr=dataset[2],
                # label=labels[i],
                linestyle=None,
                fmt="o",
                capsize=5,
                c="C0",
                # label="2012 Experiment (Gerry) $\\times" + f"{factor}$",
                label="2014 Experiment (Gerry)",
            )
        else:
            ax_scatter.errorbar(
                dataset[0],
                dataset[1],  # * factor,
                yerr=dataset[2],
                linestyle=None,
                c="C0",
                fmt="o",
                capsize=5,
            )

    ax_scatter.set_ylabel(r"$\mathrm{d} \sigma /\mathrm{d} \Omega$ ", fontsize=14)
    ax_scatter.set_xlabel("$\\theta$", fontsize=14)
    ax_scatter.set_ylim([0, 1.2 * np.max(ccs)])
    ax_scatter.set_title(
        r"$\mathrm{d} \sigma /\mathrm{d} \Omega$ vs "
        + r"$\theta$"
        + r" for Compton Scattering on ${}^6\mathrm{Li}$"
        + f" at ${MeVtmp} \\mathrm{{MeV}}$"
        + rf" with $\omega_H={omegaH}$",
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
    # saveString = "factor=" + str(factor) + "-" + saveString
    saveString = saveString + ".pdf"
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
    print("\n")
    print("Gerry experimental Data")
    print(f"{'theta':>7} | {'val':>10} | {'uncertainty':>12}")
    print(40 * "-")
    for theta, val, uncer in data:
        print(f"{theta:7.2f} | {val:10.6f} | {uncer:12.6f}")

    print("AlphaOffset=", alphaDiff)
    print("BetaOffset=", betaDiff)
    yn = input("Save it?  [y/n]" + "\n" + saveString + "\n").lower()
    if yn == "y":
        fig.savefig(saveString, format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
