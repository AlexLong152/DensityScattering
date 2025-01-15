# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import CrossSection as cc
from os import listdir
from os.path import isfile, join
import sys

sys.path.insert(1, "..")
import readDensity as rd
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda400/twobody/"
    onebody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/newfiles/lambda400/"
    energy = 60
    angle = 159
    lambdaSRG = 1.880
    Ntotmax = 14
    lambdaCut = 400
    omegaH = 14
    omegaH = None

    paramToPlot = "omegaH"

    CompareCross(
        onebody_dir,
        twobody_dir,
        paramToPlot,
        "Odelta3",
        energy=energy,
        angle=angle,
        lambdaCut=lambdaCut,
        lambdaSRG=lambdaSRG,
        Ntotmax=Ntotmax,
        omegaH=omegaH,
    )


def CompareCross(onebody_dir, twobody_dir, paramToPlot, Odeltaonebod, **kwargs):
    """
    Compares cross-sections from one-body and two-body files whose
    parameters match. Plots cross section vs. `paramToPlot`.
    """
    kwargs["theta"] = kwargs["angle"]
    # 1. Gather all one-body and two-body files
    onebody_files = [f for f in listdir(onebody_dir) if isfile(join(onebody_dir, f))]
    twobody_files = [f for f in listdir(twobody_dir) if isfile(join(twobody_dir, f))]

    # 2. Extract parameter dictionaries (WITHOUT reading in the matrix) for each file
    onebody_info = []
    matched_onebody = []
    for f in onebody_files:
        if Odeltaonebod in f:
            # returnMat=False so that we skip reading the actual matrix
            info = rd.getQuantNums(join(onebody_dir, f), returnMat=False)
            if params_match(info, kwargs, paramToPlot):
                onebody_info.append((f, info))
                matched_onebody.append(f)

    twobody_info = []
    matched_twobody = []
    for f in twobody_files:
        info = rd.getQuantNums(join(twobody_dir, f), returnMat=False)
        if params_match(info, kwargs, paramToPlot):
            twobody_info.append((f, info))
            matched_twobody.append(f)

    # 3. Define a function to decide if two parameter dictionaries match
    #    This is just an example. Adjust to match your needs.
    # print("matched_twobody=", matched_twobody)
    # print("matched_onebody=", matched_onebody)
    xs = []
    ys = []

    files = []
    for one in matched_onebody:
        onebod = rd.getQuantNums(onebody_dir + one, returnMat=True)
        x = onebod[paramToPlot]
        for two in matched_twobody:
            twobod = rd.getQuantNums(twobody_dir + two, returnMat=False)
            xTest = twobod[paramToPlot]
            if x == xTest:
                twobod = rd.getQuantNums(twobody_dir + two, returnMat=True)
                xs.append(x)
                files.append(two)
                yVal = cc.crossSection(onebod["file"], twobod["file"])["cc"]
                ys.append(yVal)

    # for i, f in enumerate(files):
    #     if xs[i] == 10:
    #         print(xs[i])
    #         print(f)

    plt.figure()
    plt.scatter(xs, ys)
    plt.xlabel(paramToPlot)
    plt.ylabel(
        r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\mathrm{ sr^{-1}}]$", fontsize=12
    )

    plt.title(
        r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\;\mathrm{ sr^{-1}} ]$ vs "
        + paramToPlot
    )
    plt.show()


def params_match(dictA, dictB, paramToPlot):
    """
    Returns True if the relevant parameters match in both dicts.
    You can compare as many or as few fields as you need.
    """
    # For example, compare lambdaSRG, Ntotmax, omegaH, etc.
    # If you have other constraints (energy, angle, lambdaCut, etc.)
    # you can incorporate them here or pass them in **kwargs.
    keys_to_compare = [
        "lambdaSRG",
        "Ntotmax",
        "omegaH",
        "lambdaCut",
        "theta",
        "energy",
    ]  # example fields
    keys_to_compare.remove(paramToPlot)
    for key in keys_to_compare:
        # If a key doesn't exist or the values differ, return False
        if key not in dictA or key not in dictB:
            return False
        if dictA[key] != dictB[key]:
            # print(dictA[key], dictB[key])
            return False

    # If all checks pass, they match
    return True


if __name__ == "__main__":
    main()
