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
import re

# Dictionary for pole values from Lenkewitz 2011
# in units of 10^-3 / m_pi+
poleDict = {}
poleDict["Eprot"] = -1.16
poleDict["Eneut"] = 2.13
poleDict["Lprot"] = -1.35
poleDict["Lneut"] = -2.41


def main():
    run3He()
    run4He()


def run4He():
    print("\n4He results")

    print(50 * "-")
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/1bod/"
    files = [folder + f for f in listdir(folder) if isfile(join(folder, f))]

    tFolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/"
    tFiles = [tFolder + f for f in listdir(tFolder) if isfile(join(tFolder, f))]
    # TwoBod_getQnums(tFiles[0])

    getResult(files, tFiles, "4He")


def run3He():
    print("3He results")
    print(50 * "-")
    oBodFiles = [
        "onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=3bd9b5d55f401de84f369724cdeb02eb176fda6928c442e1611c2ba14efdbc8b.v2.0.dat",
        "onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=fac31e8f17bc2b1a6d84411b408e05cca71aae63c04be8f6fbf3ba311acb45a3.v2.0.dat",
    ]
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/"
    oBodFiles = [folder + r"1bod/132MeV/" + f for f in oBodFiles]
    tBodFiles = [
        "twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=2-.denshash=e8f1ba5e7fafac8431e05ae0e6b63381f668664f52d39abfcea0830f420f4706.v2.0.dat",
        "twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=2-.denshash=7215fbef6c42aa1ac62c6e1948b052adddcb9b1043996bc9654f23ee00588d72.v2.0.dat",
    ]

    tBodFiles = [folder + r"2bod/132MeV/" + f for f in tBodFiles]
    getResult(oBodFiles, tBodFiles, "3He")


def getResult(oFiles, tFiles, nuc):
    """
    oFiles: onebody files
    tFiles: twobody files
    """
    E0s = []
    scatLens = []
    for oBod, tBod in zip(oFiles, tFiles):
        E0, scatLen = getScatLen(oBod, tBod)
        E0s.append(E0)
        scatLens.append(scatLen)
    E0s = np.array(E0s)
    scatLens = np.array(scatLens)
    E0, E0Spread = meanAndSpread(E0s)
    len, lenSpread = meanAndSpread(scatLens)

    # print("E0=", f"{E0} +/- {E0Spread}   10^-3/m_pi^+ units")
    # print("Scattering length=", f"{len} +/- {lenSpread}   10^-6/m^2_pi^+ units")
    print(f"E0 = {E0:.8f} +/- {E0Spread:.8f}   10^-3/m_pi^+ units")
    print(f"Scattering length = {len:.8f} +/- {lenSpread:.8f}   10^-6/m^2_pi^+ units\n")

    valsA = []
    valsB = []
    for f in tFiles:
        a, b = TwoBod_getQnums(f)
        valsA.append(a)
        valsB.append(b)
    vals = np.vstack((valsA, valsB))
    meanvals = np.mean(vals, axis=1)
    stdvals = np.max(vals, axis=1) - np.min(vals, axis=1)

    # print("Qnum A: ", meanvals[0], " +/- ", stdvals[0])
    # print("Qnum B: ", meanvals[1], " +/- ", stdvals[1])
    print("Onebody Results")
    FormFacts = []
    for f in oFiles:
        result = rd.getQuantNums(f)
        file = result["file"]
        tmp = np.array(parse_onebody(file))
        FTMinus, FTPlus, FLMinus, FLPlus = tmp
        FormFacts.append(tmp)
    FormFacts = np.vstack(FormFacts)
    if nuc == "6Li":
        print(FormFacts.real)
        print("CHECK: Do you need to get real part or imag part?")
        assert False

    means = np.mean(FormFacts, axis=0)
    maxs = np.max(FormFacts, axis=0)
    mins = np.min(FormFacts, axis=0)
    assert np.allclose(means.imag, 0, atol=1e-5)

    means = means.real
    spreads = (maxs - mins).real
    outStrs = ["F_T^{S-V}", "F_T^{S+V}", "F_L^{S-V}", "F_L^{S+V}"]

    for s, val, uncer in zip(outStrs, means, spreads):
        print(f"{s} | {val:8.4f} +/-{uncer:8.4f}")
    print("\nTwobody Results")
    if "3He" == nuc:
        factor = 2
    elif "4He" == nuc or "6Li" == nuc:
        factor = 1  # 4He is spin 0, but ignore that
    if "4He" == nuc:
        print("4He spin is zero, printing raw numbers")

    print(f"Including factor or {factor} from spin")

    print("Qnum A: ", 2 * meanvals[0], " +/- ", 2 * stdvals[0])
    print("Qnum B: ", 2 * meanvals[1], " +/- ", 2 * stdvals[1])


def parse_onebody(filepath):
    if "6Li" in filepath.split(r"/")[-1]:
        print("Check this works for 6Li")

    with open(filepath) as f:
        text = f.read()

    def get_block(eps_tag):
        """Extract the block for eps=..."""
        if f"eps={eps_tag}" not in text:
            raise ValueError(f"Section eps={eps_tag} not found")
        block = text.split(f"eps={eps_tag}")[1]
        block = block.split("############################################")[0]
        return block

    def extract_complex(block, key):
        """Extract a complex value like FMinus[...] = a + bi"""
        for line in block.splitlines():
            if key in line:
                _, rhs = line.split("=")
                rhs = rhs.strip().replace("i", "")
                parts = rhs.split()
                real = float(parts[0])
                imag = float(parts[1])
                return complex(real, imag)
        raise ValueError(f"{key} not found in block")

    # Get blocks
    eps100 = get_block("1,0,0")
    eps010 = get_block("0,1,0")
    eps001 = get_block("0,0,1")

    # Extract values
    FMinus_100 = extract_complex(eps100, "FMinus")
    FPlus_100 = extract_complex(eps100, "FPlus")
    FMinus_010 = extract_complex(eps010, "FMinus")
    FPlus_010 = extract_complex(eps010, "FPlus")
    FMinus_Long = extract_complex(eps001, "FMinus")
    FPlus_Long = extract_complex(eps001, "FPlus")

    # Average over x and y (1,0,0 and 0,1,0)
    FMinus_avg = 0.5 * (FMinus_100 + FMinus_010)
    FPlus_avg = 0.5 * (FPlus_100 + FPlus_010)

    # Return the averages plus z-polarization (0,0,1)
    return FMinus_avg, FPlus_avg, FMinus_Long, FPlus_Long


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
    OneBodpart = getE0_onebod(oneBod)
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
        if "3He" in folder:
            factor = 2
        elif "4He" or "6Li" in folder:
            factor = 1  # 4He is spin 0, but ignore that
        print(f"Or with factor or {factor} 2 from spin")

        print("Qnum A: ", 2 * meanvals[0], " +/- ", 2 * stdvals[0])
        print("Qnum B: ", 2 * meanvals[1], " +/- ", 2 * stdvals[1])
    return meanvals, stdvals


def TwoBod_getQnums(path):
    """
    Gets the quantum numbers of pion photoproduction twobody
    """
    nuc = ID_nuc(path)
    mat = rd.getQuantNums(path)["MatVals"]
    match nuc:
        case "3He":
            valuesA = np.array([-1 * mat[0, 1, 0].real, mat[1, 0, 1].imag])
            valueA = np.mean(valuesA)
            valueB = -1 * mat[2, 0, 0].real
        case "4He":
            # breakpoint()
            a, b = mat[0, 0, 0], mat[1, 0, 0]
            valueA = np.mean([a, b])
            valueB = mat[2, 0, 0]
        case "6Li":
            print("impliment me")
            assert False
    return valueA, valueB


def getE0_onebod(path):
    """
    Don't need to actually read the numerical values, since the
    fortran code does the work already
    """
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
