# -*- coding: utf-8 -*-
"""Form factor calculation utilities for pion photoproduction.

This module provides functions to extract and calculate form factors from
nuclear matrix elements and density files. It supports processing of output
files from two-body calculations with different chiral orders (Odelta2 and
Odelta4).

Author: alexl
"""

import numpy as np
import readDensity as rd
from CrossSection import getStringBetween as gsb
from os import listdir
from os.path import isfile, join, dirname, basename


def main():
    filepath = r"twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=1-.denshash=ebcf00005072b8042425e0f57af9b3bbb3995273ed88d3533d4e1e058562dac3.v2.0.dat"
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV"

    # files = [f for f in listdir(folder) if isfile(join(folder, f))]
    # for f in files:
    #     path = join(folder, f)
    #     ffs = ffFromPath(path)
    #     # print(path)
    #     print(ffs)

    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files = [f for f in files if (("Odelta2" in f) and ("j12max=1" in f))]
    diffs = []
    for f in files:
        # print("f=", f, "\n")
        path = join(folder, f)
        diff = getDiffs(path)
        diffs.append(diff)
        print("diff=", diff)
    print("")
    diffs = np.array(diffs)
    means = np.mean(diffs, axis=0)
    uncer = np.std(diffs, axis=0)
    print("means=", means)
    print("uncer=", uncer)


def getDiffs(Odelta2Path):
    """
    given the path Odelta2, returns the static contribution
    """
    diffs = diffFFS(Odelta2Path)
    # Odelta2+staticContribution=Odelta4 total
    # staticContribution=diff
    diff = diffs[1] - diffs[0]
    return diff


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
    formfacts = np.zeros(3, dtype=complex)
    if len(np.shape(oper)) > 1:
        for k, xyz in enumerate(oper):
            aves = []

            # print("###################################################################")
            # print("###################################################################")
            # print(xyz)
            # print(mat[k])
            # print("-------------------------------------------------------------------")
            for i in range(len(xyz)):
                for j in range(len(xyz)):
                    denom = xyz[i, j]
                    numer = mat[k, i, j]
                    if denom != 0:
                        tmp = numer / denom
                        # print("numer,denom=", numer, denom)
                        # print("tmp=", tmp, "\n")
                        aves.append(tmp)

            # print("")
            formfacts[k] = np.mean(np.array(aves))
    else:
        print("impliment 4He case")
    formfacts = formfacts.real
    formfacts = np.array([(formfacts[0] + formfacts[1]) / 2, formfacts[2]])
    # Magic number here, tests for spin 1/2, and if so multiplies the form factors by -1.
    # in order to match conventions from literature
    if np.shape(oper) == (3, 2, 2):
        formfacts = formfacts * -1
    return formfacts


if __name__ == "__main__":
    main()
