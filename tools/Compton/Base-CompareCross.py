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
    base = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/"
    resultSave = base + "results/"
    onebody_dir = base + r"1bod/86MeV/"
    twobody_dir = base + r"2bod/86MeV/"
    energy = 86
    angle = 180
    lambdaSRG = 1.880
    Ntotmax = 10
    lambdaCut = 550
    omegaH = 18
    plotting = "omegaH"

    xs1 = getValuesAvailable(onebody_dir, plotting)
    xs2 = getValuesAvailable(twobody_dir, plotting)
    xs = list(set(xs1) & set(xs2))
    xs = np.sort(np.array(xs))
    ccs = []
    xsPlot = []
    dString = ""
    MeVtmp = 0

    kwargsfind = {
        "energy": energy,
        "angle": angle,
        "lambdaCut": lambdaCut,
        "lambdaSRG": lambdaSRG,
        "Ntotmax": Ntotmax,
        "omegaH": omegaH,
    }
    for x in xs:
        kwargsfind[plotting] = x

        # print(kwargsfind)
        onebod, twobod, kwargs = getMatchingFiles(
            onebody_dir, twobody_dir, **kwargsfind
        )

        if onebod is not None and twobod is not None:
            yVal = cc.crossSection(onebody_dir + onebod, twobody_dir + twobod)["cc"]
            ccs.append(yVal)
            xsPlot.append(x)
            # print("xsPlot=", xsPlot)
            # print("ccs=", ccs)
        # else:
        #     print(f"For x={x}:\n   onebod={onebod}\n   twobod={twobod}\n\n")

        dString = dictString(kwargs)
        MeVtmp = kwargs["energy"]
    # make a figure with 1 row, 2 columns
    fig, (ax_scatter, ax_text) = plt.subplots(
        1, 2, figsize=7 * np.array([1.6, 1]), gridspec_kw={"width_ratios": [2, 1]}
    )

    # left subplot: your scatter
    ax_scatter.scatter(xsPlot, ccs)
    ax_scatter.set_ylabel(r"$\mathrm{d} \sigma /\mathrm{d} \Omega$ ", fontsize=14)
    ax_scatter.set_xlabel(plotting, fontsize=14)
    ax_scatter.set_title(
        r"$\mathrm{d} \sigma /\mathrm{d} \Omega$ vs "
        + plotting
        + r" for Compton Scattering on ${}^6\mathrm{Li}$"
        + f" at ${MeVtmp} \\mathrm{{MeV}}$",
        fontsize=15,
    )

    # right subplot: hide the axes, and draw text
    ax_text.axis("off")
    # place the text; use wrap=True to break long lines
    ax_text.text(
        0.5,
        0.5,
        dString,
        transform=ax_text.transAxes,
        va="center",
        ha="center",
        wrap=True,
        fontsize=16,
    )
    plt.tight_layout()
    saveString = "plot_" + plotting
    for key, value in kwargsfind.items():
        if key != plotting:
            saveString += "-" + key + "=" + str(value)
    saveString = saveString + ".pdf"
    saveString = resultSave + saveString

    plt.show()

    yn = input("Save it?  [y/n]" + "\n" + saveString + "\n").lower()
    if yn == "y":
        fig.savefig(saveString, format="pdf", bbox_inches="tight")


def dictString(d):
    keys_to_compare = ["lambdaSRG", "Ntotmax", "omegaH", "lambdaCut", "theta", "energy"]
    # keys_to_compare.remove(skip)
    out = ""
    for k in keys_to_compare:
        out += k + "=" + str(d[k]) + "\n"
    return out


def getValuesAvailable(dir, param):
    if dir[-1] != r"/":
        dir += r"/"
    files = [f for f in listdir(dir) if isfile(join(dir, f))]
    values = []
    for f in files:
        info = rd.getQuantNums(join(dir, f), returnMat=False)
        if info is not None:
            val = info[param]
            if val not in values:
                values.append(val)

    return values


def getMatchingFiles(onebody_dir, twobody_dir, **kwargs):
    """
    Given a directory onebody_dir and twobody_dir along with the kwargs
    keys_to_compare = [
        "lambdaSRG",
        "Ntotmax",
        "omegaH",
        "lambdaCut",
        "theta",
        "energy",
    ]
    returns the files that match both of these
    """
    Odeltaonebod = "Odelta3"
    if onebody_dir[-1] != r"/":
        onebody_dir += r"/"

    if twobody_dir[-1] != r"/":
        twobody_dir += r"/"

    kwargs["theta"] = kwargs["angle"]
    # 1. Gather all one-body and two-body files
    onebody_files = [f for f in listdir(onebody_dir) if isfile(join(onebody_dir, f))]
    twobody_files = [f for f in listdir(twobody_dir) if isfile(join(twobody_dir, f))]
    oneBodMatch = []
    for f in onebody_files:
        if Odeltaonebod in f:
            info = rd.getQuantNums(join(onebody_dir, f), returnMat=False)
            if params_match(info, kwargs):
                oneBodMatch.append(f)

    twoBodMatch = []
    for f in twobody_files:
        info = rd.getQuantNums(join(twobody_dir, f), returnMat=False)
        if params_match(info, kwargs):
            twoBodMatch.append(f)
    # print("oneBodMatch=", oneBodMatch)
    # print("twoBodMatch=", twoBodMatch)

    if len(oneBodMatch) > 1:
        print("More than 1 1 body file matches, picking first one ")
        for f in oneBodMatch:
            print("file=", f)

    if len(twoBodMatch) > 1:
        print("More than 1 2 body file matches, picking first one ")
        for f in twoBodMatch:
            print("file=", f)
    if len(oneBodMatch) == 0:
        oneBodMatch = None
    else:
        oneBodMatch = oneBodMatch[0]

    if len(twoBodMatch) == 0:
        twoBodMatch = None
    else:
        twoBodMatch = twoBodMatch[0]

    if twoBodMatch is None and oneBodMatch is not None:
        print("twobod is None but onebody isnt")
        print(dictString(kwargs))

    if oneBodMatch is None and twoBodMatch is not None:
        print("onebod is None but twobody isnt")
        print(dictString(kwargs))
    return oneBodMatch, twoBodMatch, kwargs


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
            if params_match_exclude(info, kwargs, paramToPlot):
                onebody_info.append((f, info))
                matched_onebody.append(f)

    twobody_info = []
    matched_twobody = []
    for f in twobody_files:
        info = rd.getQuantNums(join(twobody_dir, f), returnMat=False)
        if params_match_exclude(info, kwargs, paramToPlot):
            twobody_info.append((f, info))
            matched_twobody.append(f)

    # 3. Define a function to decide if two parameter dictionaries match
    #    This is just an example. Adjust to match your needs.
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

    return xs, ys, kwargs


def params_match(dictA, dictB):
    """
    Returns True if the relevant parameters match in both dicts.
    You can compare as many or as few fields as you need.
    """
    # For example, compare lambdaSRG, Ntotmax, omegaH, etc.
    # If you have other constraints (energy, angle, lambdaCut, etc.)
    # you can incorporate them here or pass them in **kwargs.
    keys_to_compare = ["lambdaSRG", "Ntotmax", "omegaH", "lambdaCut", "theta", "energy"]
    for key in keys_to_compare:
        # If a key doesn't exist or the values differ, return False
        if (key in dictA) and (key in dictB):
            if dictA[key] != dictB[key]:
                # print(dictA[key], dictB[key])
                return False
    # If all checks pass, they match
    return True


def params_match_exclude(dictA, dictB, paramToPlot):
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
