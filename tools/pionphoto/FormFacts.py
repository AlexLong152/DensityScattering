# -*- coding: utf-8 -*-
"""Form factor calculation utilities for pion photoproduction.

This module provides functions to extract and calculate form factors from
nuclear matrix elements and density files. It supports processing of output
files from two-body calculations with different chiral orders (Odelta2 and
Odelta4).

Author: alexl
"""

import numpy as np
from math import isinf, isnan
import readDensity as rd
from os import listdir
from os.path import isfile, join, dirname, basename
# from CrossSection import getStringBetween as gsb


def He4():
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/"

    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files1 = [
        join(folder, f) for f in files if (("Odelta2" in f) and ("j12max=1" in f))
    ]
    print("\nHelium 4 Results")
    print("With j12max=1")
    ffForFiles(files1, precision=6)


def Li6():
    folder = "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod"

    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files1 = [
        join(folder, f) for f in files if (("Odelta2" in f) and ("j12max=1" in f))
    ]
    print("\nLithium 6 Results")
    print("With j12max=1")
    ffForFiles(files1)


def He3():
    # filepath = r"twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=1-.denshash=ebcf00005072b8042425e0f57af9b3bbb3995273ed88d3533d4e1e058562dac3.v2.0.dat"
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV"
    files = [f for f in listdir(folder) if isfile(join(folder, f))]

    files1 = [
        join(folder, f) for f in files if (("Odelta2" in f) and ("j12max=1" in f))
    ]
    files2 = [
        join(folder, f) for f in files if (("Odelta2" in f) and ("j12max=2" in f))
    ]
    print("Helium 3 Results")
    print("With j12max=1")
    ffForFiles(files1)

    print("\nWith j12max=2")
    ffForFiles(files2)


def ffForFiles(files, precision=3):
    """
    Given a list of Odelta2 paths, compute form factors for
    Odelta2, Odelta4, and their differences.
    Prints a summary table with means and uncertainties.
    """
    o2_vals = []
    o4_vals = []
    diff_vals = []

    for path in files:
        diff, diffs = getDiffs(path, returnAll=True)
        # diffs has shape (2, N): [Odelta2, Odelta4]
        o2_vals.append(diffs[0])
        o4_vals.append(diffs[1])
        diff_vals.append(diff)

    o2_vals = np.array(o2_vals)
    o4_vals = np.array(o4_vals)
    diff_vals = np.array(diff_vals)

    # Means and uncertainties
    o2_mean = np.mean(o2_vals, axis=0)
    o2_unc = np.std(o2_vals, axis=0)
    o4_mean = np.mean(o4_vals, axis=0)
    o4_unc = np.std(o4_vals, axis=0)
    d_mean = np.mean(diff_vals, axis=0)
    d_unc = np.std(diff_vals, axis=0)

    print_ff_table(o2_mean, o2_unc, o4_mean, o4_unc, d_mean, d_unc, precision=precision)


def getDiffs(Odelta2Path, returnAll=False):
    """
    given the path Odelta2, returns the static contribution
    """
    diffs = diffFFS(Odelta2Path)
    # Odelta2+staticContribution=Odelta4 total
    # staticContribution=diff
    diff = diffs[1] - diffs[0]
    if not returnAll:
        return diff
    else:
        return diff, diffs


def diffFFS(path):
    """Gets form factors between Odelta2 and Odelta4 orders.

    Given a path to either an Odelta2 or Odelta4 result file, this function
    finds the corresponding file with the other order and calculates form
    factors for both.

    Parameters
    ----------
    path : str
        Path to a twobody density file containing either 'Odelta2' or
        'Odelta4' in the filename.

    Returns
    -------
    ffs : numpy.ndarray
        Array of shape (2, N) containing form factors for both Odelta2
        (first row) and Odelta4 (second row) calculations.

    Raises
    ------
    ValueError
        If the filename does not contain 'Odelta2' or 'Odelta4'.
    AssertionError
        If no unique matching file is found in the directory.
    """
    file = basename(path)
    folder = dirname(path)

    # matchString = gsb(file, "twobody", "denshash")
    matchString = file
    if "Odelta2" in matchString:
        assert "Odelta4" not in matchString
        swapString = "Odelta4"
        oldStr = "Odelta2"
        pathIsO2 = True
    elif "Odelta4" in matchString:
        # This way you can pass either Odelta2 or Odelta4 variant
        swapString = "Odelta2"
        oldStr = "Odelta4"
        pathIsO2 = False
    else:
        raise ValueError(
            f"file with name:\n{file}\n did not contain Odelta4 or Odelta2 substrings"
        )

    matchFile = matchString.replace(oldStr, swapString)
    # print("folder=", folder)
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    matchFile = [f for f in files if matchFile in f]
    assert len(matchFile) == 1, "matchfile=\n" + "\n".join(matchFile)
    matchFile = matchFile[0]

    files = [file, matchFile]
    if not pathIsO2:  # Odelta2 file always comes first
        files = [files[1], files[0]]
    # print("FormFacts.py:107 files[0]=\n", join(folder, files[0]))
    # print("")
    # print("FormFacts.py:107 files[1]=\n", join(folder, files[1]))
    # print("\n")

    # length 2 array with Odelta2 and Odelta4 results in that order
    files = [join(folder, f) for f in files]
    ffs = []

    for f in files:
        ffs.append(ffFromPath(f))

    ffs = np.array(ffs)
    return ffs


def ffFromPath(path):
    """Extract form factors from a density file path.

    This is a convenience wrapper that reads quantum numbers and matrix
    elements from a file and calculates the corresponding form factors.

    Parameters
    ----------
    path : str
        Path to the density file containing matrix elements.

    Returns
    -------
    ffs : numpy.ndarray
        Array of form factors extracted from the matrix elements.
    """
    matOut = rd.getQuantNums(path)
    mat = matOut["MatVals"]
    nuc = matOut["nuc"]
    oper = rd.spin4Nuc(nuc)
    ffs = ffFromMat(mat, oper)
    return ffs


def ffFromMat(mat, oper):
    """Calculate form factors from matrix elements and spin operators.

    Given the output matrix and the spin operator, this function computes
    the form factors by averaging over all non-zero matrix elements weighted
    by the corresponding spin operator elements.

    Parameters
    ----------
    mat : numpy.ndarray
        Matrix of nuclear matrix elements with shape (3, N, N) where N is
        the dimension of the spin space.
    oper : numpy.ndarray
        Spin operator matrices. If 2D or higher, contains operators for
        x, y, z components. Shape is (3, N, N) for multi-dimensional case.

    Returns
    -------
    formfacts : numpy.ndarray
        Array of length 2 containing the calculated form factors.
        First element is the average of x and y components, second is z.

    Notes
    -----
    The function takes the real part of the computed form factors and
    combines the x and y components by averaging them. The 4He case is
    not yet implemented.
    """

    # This just ignores divides by zero warnings, like try catch for warnings
    # ignoring invalid also supresses 0/0 errors
    with np.errstate(divide="ignore", invalid="ignore"):
        tmpMat = (mat / oper).real
    formfacts = np.zeros(3, dtype=complex)
    # if len(np.shape(oper)) > 1:
    for i in range(3):
        tmp = tmpMat[i]
        aves = [f for f in tmp.flatten() if not (isinf(f) or isnan(f))]
        formfacts[i] = np.mean(np.array(aves))

    # else:
    #     print("impliment 4He case")
    formfacts = formfacts.real
    formfacts = np.array([(formfacts[0] + formfacts[1]) / 2, formfacts[2]])
    # multiply by -1 in order to match conventions from literature
    # factor of -1 on the matrix element doesn't affect physical observables
    return formfacts * -1


def print_ff_table(o2_mean, o2_unc, o4_mean, o4_unc, d_mean, d_unc, precision=3):
    """
    Print a nice terminal table of form-factor statistics.
    """

    COLW = 24  # width of each numeric column
    PREC = precision  # number of decimals

    def fmt_pm(mu, sigma):
        s = f"{mu: .{PREC}f} ± {sigma:.{PREC}f}"
        return f"{s:^{COLW}}"  # CENTERED

    print("-" * 96)
    print(
        f"{'Component':<11} | "
        f"{'Odelta2 (mean ± σ)':^{COLW}} | "  # format string character ^ is center align, < is left align, > is right align
        f"{'Odelta4 (mean ± σ)':^{COLW}} | "
        f"{'Diff (mean ± σ)':^{COLW}}"
    )
    print("-" * 96)

    indexList = ["T", "L"]
    for i in range(len(o2_mean)):
        lab = f"F_{indexList[i]} [fm^-1]"
        print(
            f"{lab:<10} | "
            f"{fmt_pm(o2_mean[i], o2_unc[i])} | "
            f"{fmt_pm(o4_mean[i], o4_unc[i])} | "
            f"{fmt_pm(d_mean[i], d_unc[i])}"
        )

    print("-" * 96)
    print("")


if __name__ == "__main__":
    He3()
    print(96 * "#" + "\n" + 96 * "#")
    He4()
    print(96 * "#" + "\n" + 96 * "#")
    Li6()
