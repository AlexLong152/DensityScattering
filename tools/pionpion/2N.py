import numpy as np
import glob
import readDensity as rd

RESULTS_DIR = "/home/alexander/Dropbox/PionPion"


def get_scattering_lengths(filepath):
    """
    Extract the three pion scattering lengths from a single 2N output file.

    Returns (pi_minus, pi_zero, pi_plus) as real floats, taken from the
    diagonal element MatVals[extQnum-1, 0, 0] for extQnum = 1, 2, 3.
    Imaginary parts are discarded.
    """
    data = rd.getQuantNums(filepath, kind="Scattering Length")
    mat = data["MatVals"]
    pi_minus = mat[0, 0, 0].real  # extQnum=1
    pi_zero = mat[1, 0, 0].real  # extQnum=2
    pi_plus = mat[2, 0, 0].real  # extQnum=3
    return pi_minus, pi_zero, pi_plus


def filter_files(files):
    """Remove files with lambda400 or lambda550 from the list."""
    return [f for f in files if "lambda400" not in f and "lambda550" not in f]


def get_scattering_length_array(nucleus):
    """
    Build a 3 x n array of 2N pion scattering lengths for a given nucleus.

    All three pion charges (pi-, pi0, pi+) come from the same thresh folder.

    Parameters
    ----------
    nucleus : str
        Nucleus label, e.g. "3He", "4He", "6Li".

    Returns
    -------
    numpy.ndarray
        Shape (3, n) where n is the number of density files.
        Row 0: pi^-   Row 1: pi^0   Row 2: pi^+
    """
    thresh_pattern = f"{RESULTS_DIR}/results-{nucleus}/2bod/thresh/*.dat"
    thresh_files = sorted(glob.glob(thresh_pattern))
    thresh_files = filter_files(thresh_files)

    if not thresh_files:
        raise FileNotFoundError(f"No files found matching {thresh_pattern}")

    pi_minus = [get_scattering_lengths(f)[0] for f in thresh_files]
    pi_zero = [get_scattering_lengths(f)[1] for f in thresh_files]
    pi_plus = [get_scattering_lengths(f)[2] for f in thresh_files]

    return np.array([pi_minus, pi_zero, pi_plus])


def getMeanUncer(arr):
    """
    Compute the mean and uncertainty of an array of scattering lengths.

    Parameters
    ----------
    arr : numpy.ndarray
        Shape (3, n) array of scattering lengths.

    Returns
    -------
    numpy.ndarray
        Shape (3, 2) with columns [mean, total_uncertainty].
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
        print()
