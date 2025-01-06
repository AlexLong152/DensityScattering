import os
import re
import sys
import glob
import numpy as np
import readDensity as rd

folder = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/"
folder1N = folder + r"onebody/"
folder2N = folder + r"twobody/"
energy = 60
angle = 55

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


def main():
    twobody_file = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/twobody-6Li.060MeV-075deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax06omegaH24.Odelta2-j12max=2-.denshash=08cbddd96d0928d56b11000bf420a33e22d9200bf4f68972a75d8f7b57d60b4c.v2.0.dat"

    onebody_file = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/newfiles/onebody-6Li.060MeV-075deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG2.236setNtotmax14omegaH24.Odelta3-.denshash=35dc0cfe7cc3560fe2da9cd9603584f51bd0644e6f9914c4146d1d8b99deafc7.v2.0.dat"
    base_dir = os.path.dirname(os.path.dirname(twobody_file))

    Snucl = 1

    # Parse twobody filename
    params = parse_twobody_filename(twobody_file)
    nucleus = params["nucleus"]
    energy_str = params["energy"]
    angle_str = params["angle"]
    interaction = params["interaction"]

    # Construct corresponding onebody and varyA filenames
    onebody_file, varyA_files = construct_onebody_filenames(
        base_dir, nucleus, energy_str, angle_str, interaction
    )

    # Load amplitudes
    onebody_data, twobody_data, varyA_data = load_amplitudes(
        onebody_file, twobody_file, varyA_files
    )

    # Compute total amplitude
    total = compute_total(onebody_data, twobody_data, varyA_data)

    # Print shape and a sample element
    print("total shape:", total.shape)
    print("Sample total amplitude element [0,0,0]:", total[0, 0, 0])
    omega = float(energy_str[:3])
    crossSection = compute_cross_section(total, omega, Snucl, M6Li)
    print("crossSection=", crossSection, "micro Barn")
    print("Result should be on the order of 150 nanobarn=0.15 microbarn")


def parse_twobody_filename(twobody_file):
    """
    Adjusted regex to handle the denshash part.
    Example filename (from user):
    twobody-6Li.060MeV-055deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax14omegaH20.Odelta2-j12max=2-.denshash=...v2.0.dat
    """
    name = os.path.basename(twobody_file)
    pattern = (
        r"^twobody-(?P<nucleus>\w+)\.(?P<energy>\d+MeV)-(?P<angle>\d+deg)\.dens-"
        r"(?P<interaction>.+?)\.Odelta2-j12max=2-\.denshash=.*\.v2\.0\.dat$"
    )

    match = re.match(pattern, name)
    if not match:
        raise ValueError("twobody filename doesn't match expected pattern")
    return match.groupdict()


def construct_onebody_filenames(onebodyname, Odelta="Odelta3"):
    """
    Construct onebody filenames (Odelta2 and VaryA files).
    """


def load_amplitudes(onebody_file, twobody_file, varyA_files):
    """
    Load amplitudes using getQuantNums for onebody, twobody, and varyA files.
    """
    onebody_data = rd.getQuantNums(onebody_file)
    twobody_data = rd.getQuantNums(twobody_file)

    varyA_data = {}
    for key, fname in varyA_files.items():
        varyA_data[key] = rd.getQuantNums(fname)

    return onebody_data, twobody_data, varyA_data


def compute_total(onebody_data, twobody_data, varyA_data):
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
    theta_rad = np.deg2rad(theta_deg)

    cos_theta = np.cos(theta_rad)

    # Extract varyA arrays we need:
    # varyA[proton,1] = varyA_data[(1,'p')]["MatVals"]
    # varyA[proton,2] = varyA_data[(2,'p')]["MatVals"]
    # varyA[neutron,1] = varyA_data[(1,'n')]["MatVals"]
    # varyA[neutron,2] = varyA_data[(2,'n')]["MatVals"]

    varyA_p1 = varyA_data[(1, "p")]["MatVals"]
    varyA_p2 = varyA_data[(2, "p")]["MatVals"]
    varyA_n1 = varyA_data[(1, "n")]["MatVals"]
    varyA_n2 = varyA_data[(2, "n")]["MatVals"]

    # Construct the additional term:
    factor = (
        (10 ** (-4))
        * (omega**2)
        * (MeVtofm**3)
        * (
            (deltaAlphaE1p + cos_theta * deltaBetaM1p) * varyA_p1
            + (-deltaBetaM1p) * varyA_p2
            + (deltaAlphaE1n + cos_theta * deltaBetaM1n) * varyA_n1
            + (-deltaBetaM1n) * varyA_n2
        )
    )

    total = onebody + twobody + factor

    return total


def compute_cross_section(total, omega, Snucl, Mnucl):
    """
    Compute the cross section based on the provided amplitude 'total'.

    Parameters:
    -----------
    total : np.ndarray
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

    # sum over all Mf, Mi, λf, λi
    # total has dimensions (9, 3, 3)
    # sum |total|^2 over all indices
    sum_amp_sq = np.sum(np.abs(total) ** 2)

    # Formula from the snippet:
    # crosssection = fineStructure^2 * 10^7 * (1 / [2*(2Snucl+1)]) * (Mnucl/Wnucl(ω))^2 * sum(|total|^2)
    factor = (
        (fineStructure**2)
        * (10**7)
        * (1.0 / (2 * (2 * Snucl + 1)))
        * ((Mnucl / Wnucl_omega) ** 2)
    )
    cross_section = factor * sum_amp_sq

    return cross_section


if __name__ == "__main__":
    main()
