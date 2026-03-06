import numpy as np
import glob
import readDensity as rd

RESULTS_DIR = "/home/alexander/Dropbox/PionPion"


def get_scattering_lengths(filepath):
    """
    Extract the three pion scattering lengths from a single output file.

    Returns (pi_minus, pi_zero, pi_plus) as real floats, taken from the
    diagonal element MatVals[extQnum-1, 0, 0] for extQnum = 1, 2, 3.
    """
    data = rd.getQuantNums(filepath, kind="Result")
    mat = data["MatVals"]
    pi_minus = mat[0, 0, 0].real  # extQnum=1
    pi_zero = mat[1, 0, 0].real  # extQnum=2
    pi_plus = mat[2, 0, 0].real  # extQnum=3
    return pi_minus, pi_zero, pi_plus


def get_scattering_length_array(nucleus):
    """
    Build a 3 x n array of pion scattering lengths for a given nucleus.

    pi^- and pi^+ come from the thresh folder (charged-pion threshold),
    pi^0 comes from the threshpi0 folder (neutral-pion threshold).

    Parameters
    ----------
    nucleus : str
        Nucleus label, e.g. "3He", "4He", "6Li", "12C".

    Returns
    -------
    numpy.ndarray
        Shape (3, n) where n is the number of density files.
        Row 0: pi^-   Row 1: pi^0   Row 2: pi^+
    """
    thresh_pattern = f"{RESULTS_DIR}/results-{nucleus}/1bod/thresh/*.dat"
    pi0_pattern = f"{RESULTS_DIR}/results-{nucleus}/1bod/threshpi0/*.dat"

    thresh_files = sorted(glob.glob(thresh_pattern))
    pi0_files = sorted(glob.glob(pi0_pattern))

    if not thresh_files:
        raise FileNotFoundError(f"No files found matching {thresh_pattern}")
    if not pi0_files:
        raise FileNotFoundError(f"No files found matching {pi0_pattern}")

    pi_minus = [get_scattering_lengths(f)[0] for f in thresh_files]
    pi_zero = [get_scattering_lengths(f)[1] for f in pi0_files]
    pi_plus = [get_scattering_lengths(f)[2] for f in thresh_files]

    return np.array([pi_minus, pi_zero, pi_plus])


def getMeanUncer(arr):
    """
    Compute the uncertainty (standard deviation) of an array of scattering lengths.

    Parameters
    ----------
    arr : numpy.ndarray
        Array of scattering lengths.

    Returns
    -------
    float
        Standard deviation of the input array.
    """
    error = 3 / 100
    out = np.zeros((3, 2))
    for i, xs in enumerate(arr):
        stdDev = np.std(xs)
        mean = np.mean(xs)
        totalStd = np.sqrt((error * mean) ** 2 + stdDev**2)
        out[i][0] = mean
        out[i][1] = totalStd
    return out


if __name__ == "__main__":
    for nuc in ["3He", "4He", "6Li"]:
        arr = get_scattering_length_array(nuc)
        arr = getMeanUncer(arr)
        print(f"{nuc}:")
        print(f"  pi-:  {arr[0]}")
        print(f"  pi0:  {arr[1]}")
        print(f"  pi+:  {arr[2]}")
