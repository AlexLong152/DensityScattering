# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np

# from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join
from matplotlib import rcParams
from copy import copy
# import re

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"

mypath = "."
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-3:] == "dat"]


def main():
    breakStr = 20 * "#" + "\n"
    breakStr = breakStr + breakStr

    onlyfiles = [
        r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda400/twobody/twobody-6Li.060MeV-040deg.dens-chiralsmsN4LO+3nfN2LO-lambda400-lambdaSRG1.880setNtotmax06omegaH14.Odelta2-j12max=2-.v2.0.dat"
    ]
    f = onlyfiles[0]
    f = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/onebody/onebody-6Li.060MeV-040deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax06omegaH10.Odelta2-.denshash=0d7417d10a4d4a7ac0a3bf3f1ba864813227ec751d7c2b6300fcd29a0515dd1f.v2.0.dat"
    # f = r"/home/alexander/OneDrive/DensityScattering/varsub-PionPhotoProdThresh.twobody/results/4He-trivial-test.dat"
    # for f in onlyfiles:
    #     print(breakStr)
    #     getTheta(f)
    #     out = output(f)
    #     print(out)
    out = getQuantNums(f)
    print("out=", out)

    # varyFile = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/onebody/onebody-6Li.060MeV-040deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax06omegaH10.VaryA1p-.denshash=0d7417d10a4d4a7ac0a3bf3f1ba864813227ec751d7c2b6300fcd29a0515dd1f.v2.0.dat"


# def getVaryA(filename):
#     # A regex pattern to capture lines like:
#     # A1p =    137.03599000000000
#     pattern = re.compile(r"(A\d[pn])\s*=\s*([0-9.+\-E]+)")
#
#     with open(filename, "r") as f:
#         data = f.read()
#
#     for match in pattern.finditer(data):
#         amp_name = match.group(1)  # e.g. "A1p"
#         amp_value_str = match.group(2)  # e.g. "137.03599000000000"
#         try:
#             amp_value = float(amp_value_str)
#         except ValueError:
#             # If for some reason conversion fails, skip
#             continue
#
#         # Store only non-zero values
#         if amp_value != 0.0:
#             return amp_value
#
#     raise ValueError("Value not found")


def getQuantNums(filename, returnMat=True):
    """
    Reads in quantum numbers as real numbers from filename

    output values are in the order:
    out[extNum, mzp, mz]
    """
    out = {}
    if returnMat:
        with open(filename, "r") as f:
            contents = f.read()
            lines = np.array(contents.splitlines())
            indx = -1
            for i in range(len(lines)):
                if "(" in lines[i]:
                    indx = i
                    break
            if indx == -1:
                raise BadDataError(f"Bad file {filename}")

            # print("filename=", filename)
            lines = lines[indx:]
            lines = "".join(lines)
            lines = lines.split("=")
            # print("lines[0]=", lines[0])
            lines = lines[0]
            lines = lines.replace("(", "")
            lines = lines.split(")")
            lines = np.array([x.strip() for x in lines if x.strip() != ""])
            vals = np.zeros(len(lines), dtype=np.complex128)
            for i in range(len(vals)):
                tmp = lines[i].split(",")
                try:
                    vals[i] = float(tmp[0]) + 1j * float(tmp[1])
                except ValueError:
                    return None
        out["MatVals"] = vals2matrix(filename, vals)

    out["name"] = filename
    out["file"] = filename

    out["hash"] = getHash(filename)

    out["omega"] = getOmega(filename)
    out["energy"] = getOmega(filename)

    out["theta"] = getTheta(filename)
    out["angle"] = getTheta(filename)

    out["Ntotmax"] = getNtotmax(filename)
    out["omegaH"] = getOmegaH(filename)
    out["lambda"] = getLambda(filename)
    out["lambdaCut"] = getLambda(filename)
    # print("filename=", filename)
    out["lambdaSRG"] = getLambdaSRG(filename)
    out["numBodies"] = getNumBodies(filename)
    out["Odelta"] = getOdelta(filename)

    # print("len(vals)=", len(vals))
    # print(np.shape(out["MatVals"]))
    return out


class BadDataError(Exception):
    def __init__(self, message):
        super().__init__(message)


def vals2matrix(filename, vals):
    name = filename.split(r"/")[-1]

    match name:
        case _ if "6Li" in name:
            twoSpin = 2
        case _ if "4He" in name:
            twoSpin = 0
        case _ if "3He" in name:
            twoSpin = 1
        case _:
            raise ValueError("Nucleus not found")

    # if twoSpin==1:
    #     numStates=2
    # if twoSpin==2:
    #     numStates=3
    numStates = twoSpin + 1  # 2l+1 quantum states

    # twoMz = np.arange(-1 * twoSpin, twoSpin + 1, 1, dtype=int)
    # twoMzp = np.arange(-1 * twoSpin, twoSpin + 1, 1, dtype=int)
    # mzStates = np.arange(0, numStates + 1, dtype=int)
    # mzpStates = np.arange(0, numStates + 1, dtype=int)
    numextQnums = len(vals) // (numStates**2)
    out = np.zeros((numextQnums, numStates, numStates), dtype=np.complex128)
    # print(np.shape(vals))
    # print("np.shape(out)=", np.shape(out))
    # print("twoMz=", twoMz)
    allZeros = True
    i = 0
    for extNum in range(numextQnums):
        for mzp in range(numStates):
            for mz in range(numStates):
                # print("vals[i]=", vals[i])
                out[extNum, mzp, mz] = vals[i]
                if vals[i] != 0:
                    allZeros = False
                i += 1
    if allZeros:
        raise ValueError(f"All zero entries for file: {filename}")

    return out


def getStringbtwn(myString, prefix, suffix):
    myString = copy(myString)
    myString = " " + myString
    myString = myString.split(prefix)[1]
    myString = myString.split(suffix)[0]
    return myString


def getOdelta(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    btwn = getStringbtwn(filename, "omegaH", "-.")
    return btwn[3:]


def getOmega(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    omega = filename.split("MeV")[0][-3:]
    omega = int(omega)
    return omega


def getNumBodies(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    if "twobody" in filename:
        return 1
    elif "onebody" in filename:
        return 2
    else:
        raise ValueError(f"Num bodies not found for {filename}")


def getTheta(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    theta = filename.split("deg")[0][-3:]
    theta = int(theta)
    return theta


def getNtotmax(filename):
    """
    'onebody-6Li.060MeV-180deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax014omegaH24.denshash=b7a1bdaca9e2d2ec9bedc928a210baa4fc442b27027a63fbb2512686170f40f9.v2.0.h5'
    """
    filename = copy(filename)
    filename = filename.split(r"/")[-1]

    Ntotmax = getStringbtwn(filename, "setNtotmax", "omegaH")
    Ntotmax = int(Ntotmax)
    return Ntotmax


def getOmegaH(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    omegaH = filename.split("omegaH")[-1][:2]
    omegaH = int(omegaH)
    return omegaH


def getLambdaSRG(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    lambdaSRG = filename.split("lambdaSRG")[1]
    lambdaSRG = lambdaSRG[:5]
    return float(lambdaSRG)


def getLambda(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    lambdaVal = filename.split("lambda")[1]
    lambdaVal = lambdaVal[:-1]
    return int(lambdaVal)


def getHash(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    if "denshash" in filename:
        hash = filename.split("denshash=")[-1]
        hash = hash.split(".v2")[0]
        return hash
    else:
        return ""


if __name__ == "__main__":
    main()
