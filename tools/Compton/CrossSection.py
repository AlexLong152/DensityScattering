# -*- coding: utf-8 -*-
"""
@author: alexl
"""

from copy import copy
import numpy as np
import sys
from os.path import isfile, join
from os import listdir

sys.path.insert(1, "..")
import readDensity as rd

########################
# Constants
########################
MeVtofm = 1.0 / 197.327
M6Li = 6.0151228874 * 931.49432
# Based on snippet:
# BKM values at Odelta3
BKMProtonalphaE1 = 12.52
BKMNeutronalphaE1 = 12.52
BKMProtonbetaM1 = 1.25
BKMNeutronbetaM1 = 1.25

# Best-fit values given:
centralprotonalphaE1 = 10.65
centralneutronalphaE1 = 11.55
centralprotonbetaM1 = 3.15
centralneutronbetaM1 = 3.65

# Compute deltas
deltaAlphaE1p = centralprotonalphaE1 - BKMProtonalphaE1  # = 10.65 - 12.52 = -1.87
deltaAlphaE1n = centralneutronalphaE1 - BKMNeutronalphaE1  # = 11.55 - 12.52 = -0.97
deltaBetaM1p = centralprotonbetaM1 - BKMProtonbetaM1  # = 3.15 - 1.25 = 1.90
deltaBetaM1n = centralneutronbetaM1 - BKMNeutronbetaM1  # = 3.65 - 1.25 = 2.40

alphap = -1.87
alphan = -0.97
betap = 1.9
betan = 2.4


def main():
    twobody_file = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/twobody-6Li.060MeV-075deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax14omegaH24.Odelta2-j12max=2-.denshash=19f5bc14eff49c008d12dbb0bfd3d63a7c03ee4999ccc3407af828479f37a2cb.v2.0.dat"

    # Odelta3 onebody file
    # VaryAXX files inferred from onebody file string
    onebody_file = "/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/newfiles/lambda550/onebody-6Li.060MeV-075deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax14omegaH24.Odelta3-.denshash=53cdbc7a4128c4b6b4c00f980eaafb7abca39033fdf70f90d4054a2534774285.v2.0.dat"

    dSigmadOmega = crossSection(onebody_file, twobody_file)

    print("dσ/dΩ=", dSigmadOmega["cc"], "μBarn")


def crossSection(onebody_file, twobody_file, delta=0):
    """
    Calculates the differential cross section given two output files
    onebody_file is expected to be the Odelta version. The varyA file names are then automatically

    Parameters
    -----------
    onebody_file: str
        Full path to onebody output file
    twobody_file: str
        Full path to twobody output file

    Returns
    ------------
    dSigmadOmega: float
        The differential cross section in nanobarns
    """

    Z, N, spin, name = getNuc(onebody_file)

    # Construct corresponding onebody and varyA filenames
    varyA_files = getVaryAFiles(onebody_file, Z, N)
    # Load amplitudes

    onebody_data, twobody_data, varyA_data = loadAmplitudes(
        onebody_file, twobody_file, varyA_files
    )

    energy = onebody_data["omega"]  # in MeV
    energy_twobod = twobody_data["omega"]  # in MeV

    theta = onebody_data["theta"]  # in MeV
    theta_twobod = twobody_data["theta"]  # in MeV
    assert theta == theta_twobod
    assert energy == energy_twobod

    matrixValues = computeMatrix(onebody_data, twobody_data, varyA_data, delta=delta)
    dSigmadOmega = computeCrossSection(matrixValues, energy, spin, M6Li)

    returnObject = {}
    returnObject["cc"] = dSigmadOmega
    returnObject["dSigmadOmega"] = dSigmadOmega
    returnObject["onebody_file"] = onebody_file
    returnObject["twobody_file"] = twobody_file
    returnObject["varyA_files"] = varyA_files
    returnObject["onebody"] = onebody_data["MatVals"]
    returnObject["twobody"] = twobody_data["MatVals"]
    returnObject["varyA_data"] = varyA_data
    returnObject["nuc"] = name

    return returnObject


def getVaryAFiles(onebody_file, Z, N):
    file = copy(onebody_file)
    num = Z + N
    # file = onebody_file.split(r"/")[-1]
    out = []
    for n in ["n", "p"]:
        for i in range(1, num + 1):
            varyStr = "VaryA" + str(i) + n
            varyFile = file.replace("Odelta3", varyStr)
            out.append(varyFile)

    return out


def printMathematica(x):
    """
    prints an array x in the form that mathematica can read in
    """
    y = x.flatten()
    formatted_elements = [
        f"{z.real:.10f} + {z.imag:.10f} I"
        if z.imag >= 0
        else f"{z.real:.10f} - {-z.imag:.10f} I"
        for z in y
    ]
    mathematica_output = "{" + ", ".join(formatted_elements) + "}"
    print(mathematica_output)
    print("")


def loadAmplitudes(onebody_file, twobody_file, varyA_files):
    """
    Load amplitudes using getQuantNums for onebody, twobody, and varyA files.
    """
    onebody_data = rd.getQuantNums(onebody_file)
    twobody_data = rd.getQuantNums(twobody_file)
    theta_one = onebody_data["theta"]
    theta_two = twobody_data["theta"]

    omega_one = onebody_data["omega"]
    omega_two = twobody_data["omega"]

    assert theta_one == theta_two, f"theta_one={theta_one}, theta_two={theta_two}"
    assert omega_one == omega_two, f"theta_one={omega_one}, theta_two={omega_two}"

    varyA_data = {}
    for fname in varyA_files:
        key = getVaryStrFromName(fname)
        tmp = rd.getQuantNums(fname)
        varyA_data[key] = tmp
        theta = tmp["theta"]
        omega = tmp["omega"]
        assert tmp["theta"] == theta_one, (
            f"tmp[theta]={theta}, theta_onebody={theta_one}"
        )
        assert tmp["omega"] == omega_one, (
            f"tmp[theta]={omega}, theta_onebody={theta_one}"
        )

    return onebody_data, twobody_data, varyA_data


def getVaryStrFromName(file):
    file = copy(file)
    file = file.split(r"/")[-1]
    tmp = file.split("omegaH")[-1]
    assert "Vary" in tmp
    return tmp[3:10]


def computeMatrix(onebody_data, twobody_data, varyA_data, delta=0):
    """
    Compute total as per the given formula (for scalar polarisabilities only):

    total[Mf,Mi,λf,λi] = onebody + twobody
      + 10^-4 * [ ω^2 * MeVtofm^3 * (
          (δalphaE1p + cos(θ)*δbetaM1p)*varyA[proton,1] + (-δbetaM1p)*varyA[proton,2]
          + (δalphaE1n + cos(θ)*δbetaM1n)*varyA[neutron,1] + (-δbetaM1n)*varyA[neutron,2]
        ) ]

    Note:
    - proton corresponds to 'p' in our keys, neutron corresponds to 'n'.
    - varyA_data[(A,pn)]["MatVals"] gives the amplitude array.
    - We assume θ and ω come from onebody_data.

    We do not use A=3,... since spin polarisabilities = 0 now.
    """
    onebody = onebody_data["MatVals"]
    twobody = twobody_data["MatVals"]

    omega = onebody_data["omega"]  # in MeV
    theta_deg = onebody_data["theta"]  # in degrees
    theta_rad = theta_deg * np.pi / 180

    cos_theta = np.cos(theta_rad)

    varyA_1p = varyA_data["VaryA1p"]["MatVals"]
    varyA_2p = varyA_data["VaryA2p"]["MatVals"]
    varyA_1n = varyA_data["VaryA1n"]["MatVals"]
    varyA_2n = varyA_data["VaryA2n"]["MatVals"]

    # varyStrs = ["VaryA1p", "VaryA2p", "VaryA1n", "VaryA2n"]
    # for i, strV in enumerate(varyStrs):
    #     print(varyA_data[strV]["name"])
    #     print(vals[i].flatten())
    #     print("\n")
    if isinstance(delta, (int, float)):
        alphap_tmp = alphap + delta
        alphan_tmp = alphan + delta
        betap_tmp = betap + delta
        betan_tmp = betan + delta
    else:
        deltap, deltan = delta
        alphap_tmp = alphap + deltap
        alphan_tmp = alphan + deltan
        betap_tmp = betap + deltap
        betan_tmp = betan + deltan

    tmp1 = (alphap_tmp + cos_theta * betap_tmp) * varyA_1p - betap_tmp * varyA_2p
    tmp2 = (alphan_tmp + cos_theta * betan_tmp) * varyA_1n - betan_tmp * varyA_2n
    polarizability = (omega**2) * (MeVtofm**3) * (10**-4) * (tmp1 + tmp2)
    total = onebody + twobody + polarizability
    total = total / MeVtofm
    return total


def computeCrossSection(matrix, omega, Snucl, Mnucl):
    """
    Compute the cross section based on the provided amplitude 'total' in units of microbarn

    Parameters:
    -----------
    matrix : np.ndarray
        The total amplitude array. Assumed shape: (9, 3, 3) for λf,λi ∈ {-1,1,2} and Mf,Mi ∈ {-1,0,1}.
        Here:
        - The first dimension: all (λf, λi) combinations. (3 polarizations each → 3*3=9)
        - The second dimension: Mf states (3 states for S=1)
        - The third dimension: Mi states (3 states for S=1)

    omega : float
        The photon energy in MeV.

    Snucl : int
        Spin of the nucleus.

    Mnucl : float
        Mass of the nucleus in MeV.

    fineStructure : float
        Fine structure constant (default = 1/137.03599)

    Returns:
    --------
    cross_section : float
        Computed cross section value.
    """
    fineStructure = 1 / 137.03599
    # Number of nuclear spin states: 2Snucl+1
    # Wnucl(ω) = ω + sqrt(Mnucl² + ω²)
    Wnucl_omega = omega + np.sqrt(Mnucl**2 + omega**2)

    mat_flat = matrix.flatten()
    # total = np.sum(mat_flat)
    # sum_amp_sq = np.float64(total * total.conj())
    sum_amp_sq = np.sum(mat_flat * mat_flat.conj())
    sum_amp_sq = sum_amp_sq.real

    factor = (
        (fineStructure**2)
        * (10**7)
        * (1.0 / (2 * (2 * Snucl + 1)))
        * ((Mnucl / Wnucl_omega) ** 2)
    )
    cross_section = factor * sum_amp_sq
    return float(cross_section)


def getNuc(filename):
    file = filename.split(r"/")[-1]
    prefix = file[:30]
    if "6Li" in prefix:
        name = "6Li"
        Z = 3
        N = 3
        S = 1
    elif "3He" in prefix:
        name = "3He"
        Z = 2
        N = 1
        S = 0.5
    else:
        raise ValueError("Something went wrong identifying nucleus")
    return (Z, N, S, name)


def ccForDict(onebody_dir, twobody_dir, Odeltaonebod="Odelta3", delta=0, **kwargs):
    """
    Given a set of parameters in kwargs and the directories the output files are in
    returns the cross section for the given parameters, for example:

    ccVal = ccForDict(
        onebody_dir,
        twobody_dir,
        delta=delta,
        energy=energy,
        angle=theta,
        lambdaCut=lambdaCut,
        lambdaSRG=lambdaSRG,
        Ntotmax=Ntotmax,
        omegaH=omegaH,
    )

    Parameters
    ----------
    onebody_dir: str
        the onebody directory
    twobody_dir: str
        the twobody directory
    Odeltaonebod: str,optional
        "Odelta3" or "Odelta2"
    delta: float or int
        the value to shift the polarizabilities
    Returns
    -------
    ccVal: float
        the cross section

    """
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
    # if kwargs["theta"] == 55 and kwargs["Ntotmax"] == 12:
    #     print("In CrossSection.py ccForDict debug")
    #     print("onebody_dir=", onebody_dir)
    #     print("twobody_dir=", twobody_dir)
    #     print("matched_onebody=", matched_onebody)
    #     print("matched_twobody=", matched_twobody)
    if len(matched_twobody) == 0 or len(matched_onebody) == 0:
        return None
    one = matched_onebody[0]
    two = matched_twobody[0]
    onebod = rd.getQuantNums(onebody_dir + one, returnMat=True)
    twobod = rd.getQuantNums(twobody_dir + two, returnMat=False)
    ccVal = crossSection(onebod["file"], twobod["file"], delta=delta)["cc"]
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
