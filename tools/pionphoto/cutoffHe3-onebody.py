# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import parseThreshOneBody as pto
from os import listdir
from os.path import isfile, join

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def getPaths(folder):
    if folder[-1] != r"/":
        folder += r"/"

    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    paths = []
    for f in files:
        paths.append(folder + f)
    return paths


def main():
    folder = "/home/alexander/Dropbox/pionphotop/3He/threshold/onebody/"
    paths = getPaths(folder)
    fig, axs = plt.subplots(2)
    Trans = []
    Long = []
    for i, p in enumerate(paths):
        myDict = pto.He3Parse(p)
        # if i == 0:
        #     for k, v in myDict.items():
        #         print(f"{k}: {v}")
        Trans.append(myDict[1])
        Long.append(myDict[3])
        if i == 0:
            axs[0].scatter(
                myDict["lambdaCut"],
                float(myDict[1][0]),
                c="b",
                marker="*",
                label="Transverse",
            )
            axs[0].scatter(
                myDict["lambdaCut"],
                float(myDict[3][0]),
                c="r",
                marker="|",
                label="Longitudinal",
            )

            axs[1].scatter(
                myDict["lambdaCut"],
                float(myDict[1][1]),
                c="b",
                marker="*",
                label="Transverse",
            )
            axs[1].scatter(
                myDict["lambdaCut"],
                float(myDict[3][1]),
                c="r",
                marker="|",
                label="Longitudinal",
            )
        else:
            axs[0].scatter(myDict["lambdaCut"], float(myDict[1][0]), c="b", marker="*")
            axs[0].scatter(myDict["lambdaCut"], float(myDict[3][0]), c="r", marker="|")

            axs[1].scatter(myDict["lambdaCut"], float(myDict[1][1]), c="b", marker="*")
            axs[1].scatter(myDict["lambdaCut"], float(myDict[3][1]), c="r", marker="|")
    axs[0].set_xlabel(r"$\Lambda_{\mathrm{NN}}$")
    axs[1].set_xlabel(r"$\Lambda_{\mathrm{NN}}$")
    fig.supylabel("Matrix Element Value")
    fig.suptitle(r"${}^3\mathrm{He}$ Pion-Photoproduction Marix Elements")
    plt.legend()
    # saveFolder = r"/home/alexander/Dropbox/pionphotop/results/"
    # save_path = saveFolder + "He3_OneBody_cutoff.pdf"
    # fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    Trans = np.array(Trans).T
    Long = np.array(Long).T
    SigmaT = np.array([maxDiff(a) for a in Trans]) / 2
    SigmaL = np.array([maxDiff(a) for a in Long]) / 2
    meanT = np.mean(Trans, axis=1)
    meanL = np.mean(Long, axis=1)
    strAr = ["^S+V=", "^S-V="]
    print("------------ Transverse ------------")
    for i in range(2):
        print("F_T" + strAr[i], np.round(meanT[i], 5), "+/-", np.round(SigmaT[i], 5))
    print("\n")
    print("------------ Longitudinal ------------")
    for i in range(2):
        print("F_L" + strAr[i], np.round(meanL[i], 5), "+/-", np.round(SigmaL[i], 5))
    print("\n\nRaw Values")
    print(Trans)
    print(Long)


def maxDiff(a):
    # print("a=", a)
    vmin = a[0]
    dmax = 0
    for i in range(len(a)):
        if a[i] < vmin:
            vmin = a[i]
        elif a[i] - vmin > dmax:
            dmax = a[i] - vmin
    return dmax


if __name__ == "__main__":
    main()
