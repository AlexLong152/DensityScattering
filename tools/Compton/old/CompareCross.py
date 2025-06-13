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
from scipy.optimize import curve_fit

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
printFiles = False  # plots the names of the
fitIt = True


def fitFunc(x, b, c, d):
    """
    Define a function passed to plotMe that will be fitted at the angles
    you specify in thetaPlot
    if you change the number of parameters this function takes
    you will also have to change "bounds" to be the correct length
    """
    return b + c * x**-2.0 + d * x**-3.0


def main():
    varyOmegaHPlot()


def varyNtotmaxPlot():
    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/"
    onebody_dir = twobody_dir + r"../onebody/"

    energy = 60
    angle = 159
    lambdaSRG = 1.880
    Ntotmax = 14
    lambdaCut = 550
    # omegaH = 24
    omegaHs = [14, 18, 22]
    angle = None
    paramToPlot = "angle"

    assert str(lambdaCut) in twobody_dir
    assert str(lambdaCut) in onebody_dir
    # plt.figure()


def varyOmegaHPlot():
    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/"
    onebody_dir = twobody_dir + r"../onebody/"

    energy = 60
    angle = 159
    lambdaSRG = 1.880
    Ntotmax = 14
    lambdaCut = 550
    # omegaH = 24
    omegaHs = [14, 18, 22]
    angle = None
    paramToPlot = "angle"

    assert str(lambdaCut) in twobody_dir
    assert str(lambdaCut) in onebody_dir
    # plt.figure()
    for i, omegaH in enumerate(omegaHs):
        xs, ys, info = compareCross(
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

        # out, _ = curve_fit(fitFunc, xs, ys, maxfev=10000000)  # "out" is an array
        # xsTmp = np.linspace(np.min(xs), np.max(xs), num=100)
        # yDatafit = fitFunc(xsTmp, *out)  # automatically unpacks "out"

        # ScatterPlot(xs, ys, paramToPlot, info, fit=np.array([xsTmp, yDatafit]))
        ScatterPlot(
            xs,
            ys,
            paramToPlot,
            info,
            fit=None,
            iteration=i,
            label="omegaH=" + str(omegaH),
        )
    plt.show()


def ScatterPlot(xs, ys, paramToPlot, info, fit=None, iteration=None, label=None):
    # print("paramToPlot=", paramToPlot)
    # print(paramToPlot == "angle")
    colors = ["b", "g", "r", "c", "m", "y", "k"]
    markers = [".", ",", "o", "v", "^", "<", ">"]
    if iteration is None:
        iteration = 0

    if paramToPlot == "angle":
        # print("in")
        tmp = r"$\theta$"
    else:
        tmp = paramToPlot
    plt.title(
        r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\;\mathrm{ sr^{-1}} ]$ vs " + tmp
    )
    plt.ylabel(r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\;\mathrm{ sr^{-1}} ]$")

    infoStr = ""
    for key, value in info.items():
        if not isinstance(value, type(None)) and key != "theta":
            infoStr += str(key) + "=" + str(value) + ", "

    plt.xlabel(tmp + "\n" + infoStr[:-2])
    if not isinstance(fit, type(None)):
        plt.plot(
            fit[0],
            fit[1],
            label="Fitted Curve",
        )
        plt.scatter(xs, ys, label="True Data")

    else:
        if label is not None:
            plt.scatter(
                xs, ys, c=colors[iteration], marker=markers[iteration], label=label
            )
            plt.legend()
        else:
            plt.scatter(
                xs, ys, c=colors[iteration], marker=markers[iteration], label=label
            )
    plt.tight_layout()
    # plt.show()


def compareCross(onebody_dir, twobody_dir, paramToPlot, Odeltaonebod, **kwargs):
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

    assert len(matched_onebody) > 0
    assert len(matched_twobody) > 0

    # print("matched_twobody=", matched_twobody)
    # print("matched_onebody=", matched_onebody)
    xs = []
    ys = []
    matchedFiles = []
    for one in matched_onebody:
        onebod = rd.getQuantNums(onebody_dir + one, returnMat=True)
        x = onebod[paramToPlot]
        for two in matched_twobody:
            twobod = rd.getQuantNums(twobody_dir + two, returnMat=False)
            xTest = twobod[paramToPlot]
            if x == xTest:
                twobod = rd.getQuantNums(twobody_dir + two, returnMat=True)
                xs.append(x)
                matchedFiles.append((one, two))
                yVal = cc.crossSection(onebod["file"], twobod["file"])["cc"]
                ys.append(yVal)

    xs = np.array(xs)
    ys = np.array(ys)
    matchedFiles = np.array(matchedFiles)

    indx = np.argsort(xs)
    xs = xs[indx]
    ys = ys[indx]

    matchedFiles = matchedFiles[indx]
    fileStr = ""
    for i, f in enumerate(matchedFiles):
        one, two = matchedFiles[i]
        j = i + 1
        fileStr += f"{j}. {one}\n  {two}\n\n"

    if printFiles:
        print(fileStr)

    return xs, ys, kwargs


def params_match(dictA, dictB, paramToPlot):
    """
    Returns True if the relevant parameters match in both dicts.
    You can compare as many or as few fields as you need.
    """
    # For example, compare lambdaSRG, Ntotmax, omegaH, etc.
    # If you have other constraints (energy, angle, lambdaCut, etc.)
    # you can incorporate them here or pass them in **kwargs.
    # print("paramToPlot=", paramToPlot)
    keys_to_compare = [
        "lambdaSRG",
        "Ntotmax",
        "omegaH",
        "lambdaCut",
        "angle",
        "energy",
    ]  # example fields
    eps = 1.0
    keys_to_compare.remove(paramToPlot)
    for key in keys_to_compare:
        # If a key doesn't exist or the values differ, return False
        if key not in dictA or key not in dictB:
            return False
        if key != "theta":
            if dictA[key] != dictB[key]:
                # print(dictA[key], dictB[key])
                return False
        else:
            if abs(dictA["theta"] - dictB["theta"]) > eps:
                return False

    # If all checks pass, they match
    return True


if __name__ == "__main__":
    main()
