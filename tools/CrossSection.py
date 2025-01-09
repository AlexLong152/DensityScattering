from copy import copy
import numpy as np
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
    Z, N, spin = getNuc(onebody_file)

    # Construct corresponding onebody and varyA filenames
    varyA_files = getVaryAFiles(onebody_file, Z, N)
    # Load amplitudes

    onebody_data, twobody_data, varyA_data = loadAmplitudes(
        onebody_file, twobody_file, varyA_files
    )

    # print("varyA_data=", varyA_data)
    onebod = onebody_data["MatVals"]
    twobod = twobody_data["MatVals"]
    """
    print("Path to onebody file:\n", onebody_file, sep="")
    print("onebod=\n", onebod.flatten())
    print("\n")

    print("Path to twobody file:\n", twobody_file, sep="")
    print("twobod=\n", twobod.flatten())
    print("\n")

    # for key, _ in varyA_data.items():
    #     print(f"Path to {key} file:\n", varyA_data[key]["name"], sep="")
    #     print(f"{key} data:\n", varyA_data[key]["MatVals"].flatten())
    #     print("\n")

    key = "VaryA1n"
    print(f"Path to {key} file:\n", varyA_data[key]["name"], sep="")
    print(f"{key} data:\n", varyA_data[key]["MatVals"].flatten())
    print("\n")
    """
    energy = onebody_data["omega"]  # in MeV
    energy_twobod = twobody_data["omega"]  # in MeV

    theta = onebody_data["theta"]  # in MeV
    theta_twobod = twobody_data["theta"]  # in MeV
    assert theta == theta_twobod
    assert energy == energy_twobod

    matrixValues = computeMatrix(onebody_data, twobody_data, varyA_data)
    crossSection = computeCrossSection(matrixValues, energy, spin, M6Li)

    # print("Mathematica friendly output of matrix values")
    # printMathematica(matrixValues)

    print(f"At {theta} degrees, {energy} MeV")
    print("dσ/dΩ=", crossSection, "μBarn")


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

    varyA_data = {}
    for fname in varyA_files:
        key = getVaryStrFromName(fname)
        varyA_data[key] = rd.getQuantNums(fname)

    return onebody_data, twobody_data, varyA_data


def getVaryStrFromName(file):
    file = copy(file)
    file = file.split(r"/")[-1]
    tmp = file.split("omegaH")[-1]
    assert "Vary" in tmp
    return tmp[3:10]


def computeMatrix(onebody_data, twobody_data, varyA_data):
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
    # manual data insert from Greisshammer
    # 550MeV at 75 degrees, 60MeV

    omega = 60
    theta_deg = 75
    theta_rad = theta_deg * np.pi / 180
    cos_theta = np.cos(theta_rad)

    # varyStrs = ["VaryA1p", "VaryA2p", "VaryA1n", "VaryA2n"]
    # for i, strV in enumerate(varyStrs):
    #     print(varyA_data[strV]["name"])
    #     print(vals[i].flatten())
    #     print("\n")

    tmp1 = (alphap + cos_theta * betap) * varyA_1p - betap * varyA_2p
    tmp2 = (alphan + cos_theta * betan) * varyA_1n - betan * varyA_2n
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
    return cross_section


def getNuc(filename):
    file = filename.split(r"/")[-1]
    prefix = file[:30]
    if "6Li" in prefix:
        Z = 3
        N = 3
        S = 1
    elif "3He" in prefix:
        Z = 2
        N = 1
        S = 0.5
    else:
        raise ValueError("Something went wrong identifying nucleus")
    return Z, N, S


if __name__ == "__main__":
    main()
