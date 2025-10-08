# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import sys

sys.path.insert(1, "..")
import readDensity as rd
from os import listdir
from os.path import isfile, join

# Dictionary for pole values from Lenkewitz 2011
# in units of 10^-3 / m_pi+
poleDict = {}
poleDict["Eprot"] = -1.16
poleDict["Eneut"] = 2.13
poleDict["Lprot"] = -1.35
poleDict["Lneut"] = -2.41


def main():
    # folder = r"/home/alexander/OneDrive/densities-4He/1Ndensities/133MeV/"
    # files = [folder + f for f in listdir(folder) if isfile(join(folder, f))]
    # TwoBod_getQnums(files[0])
    run3He()


def run3He():
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/"
    oBodFiles = [
        "onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=3bd9b5d55f401de84f369724cdeb02eb176fda6928c442e1611c2ba14efdbc8b.v2.0.dat",
        "onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=fac31e8f17bc2b1a6d84411b408e05cca71aae63c04be8f6fbf3ba311acb45a3.v2.0.dat",
    ]
    oBodFiles = [folder + r"1bod/132MeV/" + f for f in oBodFiles]
    print("oBodFiles[0]=", oBodFiles[0])
    tBodFiles = [
        "twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=2-.denshash=e8f1ba5e7fafac8431e05ae0e6b63381f668664f52d39abfcea0830f420f4706.v2.0.dat",
        "twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=2-.denshash=7215fbef6c42aa1ac62c6e1948b052adddcb9b1043996bc9654f23ee00588d72.v2.0.dat",
    ]

    tBodFiles = [folder + r"2bod/132MeV/" + f for f in tBodFiles]
    E0s = []
    scatLens = []
    for oBod, tBod in zip(oBodFiles, tBodFiles):
        E0, scatLen = getScatLen(oBod, tBod)
        E0s.append(E0)
        scatLens.append(scatLen)
    # print("boop")
    # print("scatLens=", scatLens)
    # print("E0s=", E0s)
    E0s = np.array(E0s)
    scatLens = np.array(scatLens)
    E0, E0Spread = meanAndSpread(E0s)
    len, lenSpread = meanAndSpread(scatLens)

    # print("E0=", f"{E0} +/- {E0Spread}   10^-3/m_pi^+ units")
    # print("Scattering length=", f"{len} +/- {lenSpread}   10^-6/m^2_pi^+ units")
    print(f"E0 = {E0:.3f} +/- {E0Spread:.3f}   10^-3/m_pi^+ units")
    print(f"Scattering length = {len:.3f} +/- {lenSpread:.3f}   10^-6/m^2_pi^+ units")


def meanAndSpread(arr):
    x = np.mean(arr)
    spread = np.max(arr) - np.mean(arr)
    return x, spread


def getScatLen(oneBod, twoBod):
    """
    Returns E_{0+} as defined in lenkewitz 2011
    in units of 10^-6/M^2_{pi^+}

    E0plus in slightly more normal units of 10^-3/M_{pi^+}
    """
    OneBodpart = readOneBod(oneBod)
    TwoBodResult = 2 * TwoBod_getQnums(twoBod)[0]
    # OneBodpart = 1.71
    # TwoBodResult = -29.3
    K2N = 0.135
    TwoBodpart = TwoBodResult * K2N
    E0plus = OneBodpart + TwoBodpart
    return E0plus, E0plus**2


def TwoBod_AveAndSpread(folder, verbose=False):
    assert folder[-1] == r"/"
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    valsA = []
    valsB = []
    for f in files:
        a, b = TwoBod_getQnums(folder + f)
        valsA.append(a)
        valsB.append(b)
    vals = np.vstack((valsA, valsB))
    meanvals = np.mean(vals, axis=1)
    stdvals = np.max(vals, axis=1) - np.min(vals, axis=1)

    if verbose:
        print("Qnum A: ", meanvals[0], " +/- ", stdvals[0])
        print("Qnum B: ", meanvals[1], " +/- ", stdvals[1])
        print("Or with factor of 2 from lenkewitz")

        print("Qnum A: ", 2 * meanvals[0], " +/- ", 2 * stdvals[0])
        print("Qnum B: ", 2 * meanvals[1], " +/- ", 2 * stdvals[1])
    return meanvals, stdvals


def TwoBod_getQnums(path):
    """
    Gets the transverse quantum numbers of pion photoproduction
    """
    nuc = ID_nuc(path)
    mat = rd.getQuantNums(path)["MatVals"]
    match nuc:
        case "3He":
            valuesA = np.array([-1 * mat[0, 1, 0].real, mat[1, 0, 1].imag])
            valueA = np.mean(valuesA)
            valueB = -1 * mat[2, 0, 0].real
            return valueA, valueB
        case "4He":
            breakpoint()
            print("impliment me")
            assert False
        case "6Li":
            print("impliment me")
            assert False


def TwoBod_getQnums_Longitud(path):
    """
    Gets the longitudinal quantum numbers of pion photoproduction
    """
    mat = rd.getQuantNums(path)["MatVals"]

    valuesA = np.array([-1 * mat[0, 1, 0].real, mat[1, 0, 1].imag])
    valueA = np.mean(valuesA)
    valueB = -1 * mat[2, 0, 0].real
    return valueA, valueB


def readOneBod(path):
    with open(path, "r", encoding="utf-8") as file:
        content = file.read()
    content = content.split("Average E_0+^1N=")[1]
    content = content.split("\n")
    content = float(content[0])
    return content


def ID_nuc(path):
    path = path.split(r"/")
    name = path[-1]
    subString = name.split("body")[1]
    checkStr = subString.split("MeV")[0]

    match checkStr:
        case s if "3He" in s:
            return "3He"
        case s if "4He" in s:
            return "4He"
        case s if "6Li" in s:
            return "6Li"
        case _:
            raise ValueError(f"checkString={checkStr} did not match pattern")


if __name__ == "__main__":
    main()
