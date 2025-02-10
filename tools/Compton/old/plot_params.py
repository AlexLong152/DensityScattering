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


def main():
    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/"
    onebody_dir = twobody_dir + r"../onebody/"

    energy = 60
    angle = 159
    lambdaSRG = 1.880
    lambdaCut = 550
    omegaH = 12
    thetas = np.array([0, 40, 55, 75, 90, 110, 125, 145, 159, 180])
    assert str(lambdaCut) in twobody_dir
    assert str(lambdaCut) in onebody_dir
    # plt.figure()
    Ntotmaxs = np.array([6, 8, 10, 12, 14])
    omegaHs = [14, 18, 22]
    Ntotmax = 14
    plt.figure()
    colors = ["b", "g", "r", "c", "m", "y", "k"]
    markers = [".", ",", "o", "v", "^", "<", ">"]
    for i, Ntotmax in enumerate(Ntotmaxs):
        xs = []
        ys = []
        for j, angle in enumerate(thetas):
            ccVal = ccForDict(
                onebody_dir,
                twobody_dir,
                energy=energy,
                angle=angle,
                lambdaCut=lambdaCut,
                lambdaSRG=lambdaSRG,
                Ntotmax=Ntotmax,
                omegaH=omegaH,
            )
            if ccVal is not None:
                # xs.append(angle + (Ntotmax - 6) / 2)
                xs.append(angle)
                ys.append(ccVal)
        plt.scatter(
            xs,
            ys,
            c=colors[i],
            marker=markers[i],
            label="Ntotmax=" + str(Ntotmax),
        )
    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\;\mathrm{ sr^{-1}} ]$")
    plt.legend()
    plt.show()


def compareCross(onebody_dir, twobody_dir, paramToPlot, Odeltaonebod, **kwargs):
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

    # matchedFiles = matchedFiles[indx]
    # fileStr = ""
    # for i, f in enumerate(matchedFiles):
    #     one, two = matchedFiles[i]
    #     j = i + 1
    #     fileStr += f"{j}. {one}\n  {two}\n\n"

    return xs, ys, kwargs


def ccForDict(onebody_dir, twobody_dir, Odeltaonebod="Odelta3", delta=0, **kwargs):
    kwargs["theta"] = kwargs["angle"]
    # 1. Gather all one-body and two-body files
    onebody_files = [f for f in listdir(onebody_dir) if isfile(join(onebody_dir, f))]
    twobody_files = [f for f in listdir(twobody_dir) if isfile(join(twobody_dir, f))]

    # 2. Extract parameter dictionaries (WITHOUT reading in the matrix) for each file
    onebody_info = []
    matched_onebody = []

    for f in onebody_files:
        # print("f=", f)
        if Odeltaonebod in f:
            # returnMat=False so that we skip reading the actual matrix
            info = rd.getQuantNums(join(onebody_dir, f), returnMat=False)
            if params_match_free(info, kwargs):
                onebody_info.append((f, info))
                matched_onebody.append(f)
    # print("matched_onebody=", matched_onebody)
    twobody_info = []
    matched_twobody = []
    for f in twobody_files:
        info = rd.getQuantNums(join(twobody_dir, f), returnMat=False)
        if params_match_free(info, kwargs):
            twobody_info.append((f, info))
            matched_twobody.append(f)
    # print("matched_twobody=", matched_twobody)
    # assert len(matched_onebody) > 0
    # assert len(matched_twobody) > 0
    if len(matched_twobody) == 0 or len(matched_onebody) == 0:
        return None
    one = matched_onebody[0]
    two = matched_twobody[0]
    onebod = rd.getQuantNums(onebody_dir + one, returnMat=True)
    twobod = rd.getQuantNums(twobody_dir + two, returnMat=False)
    ccVal = cc.crossSection(onebod["file"], twobod["file"], delta=delta)["cc"]
    return ccVal


def params_match_free(dictA, dictB):
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
