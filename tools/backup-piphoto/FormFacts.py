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
import array_to_latex as a2l
# from CrossSection import getStringBetween as gsb


def OneBodFFs(folder, name):
    """
    Gets and prints the onebody form factors
    """
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files1 = [join(folder, f) for f in files]
    print(f"\n{name} One Body Results")
    data = OneBodFFsForFiles(files1, divide=False)
    print("data=", data)
    createOneBodTable(data)


def TwoBodFFs(folder, name):
    """
    Gets and prints the twobody form factors
    """
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files1 = [
        join(folder, f) for f in files if (("Odelta2" in f) and ("j12max=1" in f))
    ]
    print(f"\n{name} Two Body Results")
    ffForFiles(
        files1,
        precision=6,
        func=combinedFFs,
        labels=["F_T", "F_L", "E_0+", "L_0+"],
        name=name,
    )

    # print("\nWith j12max=2")
    # ffForFiles(files2)


def ELFromPath(path):
    """
    Extract the E_{0+} and L_{0+} values from the output
    """
    matOut = rd.getQuantNums(path, kind="Matrix")
    mat = matOut["MatVals"]
    nuc = matOut["nuc"]
    oper = rd.spin4Nuc(nuc)
    ELs = ffFromMat(mat, oper)
    return ELs


def combinedFFs(path):
    """
    Extract both form factors and E_{0+}, L_{0+} values from the output
    """
    ffs = ffFromPath(path)
    els = ELFromPath(path)
    return np.concatenate([ffs, els])


def scatMatFromPath(path):
    """Extract E_0+^1N and L_0+^1N directly from output file text.

    These values are pre-computed and printed in the output file as:
    E_0+^1N=  1.751478 (averaged)
    L_0+^1N= -1.932005

    Parameters
    ----------
    path : str
        Path to the output file.

    Returns
    -------
    ffs : numpy.ndarray
        Array of length 2 containing [E_0+^1N, L_0+^1N].
    """
    with open(path, "r") as f:
        text = f.read()

    # Find the lines with E_0+^1N and L_0+^1N
    import re

    e0p_match = re.search(r"E_0\+\^1N=\s*([+-]?\d+\.?\d*)", text)
    l0p_match = re.search(r"L_0\+\^1N=\s*([+-]?\d+\.?\d*)", text)

    if not e0p_match or not l0p_match:
        raise ValueError(f"Could not find E_0+^1N or L_0+^1N in file {path}")

    E0p = float(e0p_match.group(1))
    L0p = float(l0p_match.group(1))

    return np.array([E0p, L0p])


def ffFromPath(path, kind="FormFactors", divide=True):
    """Extract form factors from an output file path

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
    # Special case: ScatMat values are pre-computed in the file
    if kind == "ScatMat" and "onebody in path":
        return scatMatFromPath(path)

    matOut = rd.getQuantNums(path, kind=kind)
    mat = matOut["MatVals"]
    nuc = matOut["nuc"]
    oper = rd.spin4Nuc(nuc)
    ffs = ffFromMat(mat, oper, divide=divide)
    return ffs


def OneBodFFsForFiles(files, divide=False):
    """
    Given a list of Odelta2 paths, compute form factors for
    Odelta2, Odelta4, and their differences.
    Prints a summary table with means and uncertainties.
    """
    Forms = ["F^{S+V}", "F^{S-V}", "ScatMat"]
    Results = []
    for form in Forms:
        ffValue = []
        for path in files:
            tmp = ffFromPath(path, kind=form, divide=divide)
            if form != "ScatMat":
                tmp = tmp * -1  # cancel out ffFrom -1 factor
            ffValue.append(tmp)
        Results.append(ffValue)
    Results = np.array(Results)
    # print("Results=\n", Results)
    # print("Results[2]=\n", Results[2])
    out = []
    for i, mat in enumerate(Results):
        # print("\n\n")
        # print("mat=\n", mat)
        mean = np.mean(Results[i], axis=0)
        mean1, mean2 = mean
        stdDev = spread(Results[i])
        std1, std2 = stdDev
        out.append([mean1, std1])
        out.append([mean2, std2])
        # print("FormFacts.py:158 mean=\n", mean)
        # print("stdDev=", stdDev)
    return np.array(out)


def ffForFiles(files, precision=3, func=ffFromPath, labels=None, name=None):
    """
    Given a list of Odelta2 paths, compute form factors for
    Odelta2, Odelta4, and their differences.
    Prints a summary table with means and uncertainties.
    """
    o2_vals = []
    o4_vals = []
    diff_vals = []

    for path in files:
        # diff, diffs = getDiffs(path, returnAll=True, func=func)
        diff, o2, o4 = getDiffs(path, returnAll=True, func=func)
        # diffs has shape (2, N): [Odelta2, Odelta4]
        o2_vals.append(o2)
        o4_vals.append(o4)
        diff_vals.append(diff)

    o2_vals = np.array(o2_vals)
    o4_vals = np.array(o4_vals)
    diff_vals = np.array(diff_vals)

    # Means and uncertainties
    o2_mean = np.mean(o2_vals, axis=0)
    o4_mean = np.mean(o4_vals, axis=0)
    d_mean = np.mean(diff_vals, axis=0)
    # o2_unc = np.std(o2_vals, axis=0)
    # o4_unc = np.std(o4_vals, axis=0)
    # d_unc = np.std(diff_vals, axis=0)

    o2_unc = spread(o2_vals)
    o4_unc = spread(o4_vals)
    d_unc = spread(diff_vals)
    print_ff_table(
        o2_mean,
        o2_unc,
        o4_mean,
        o4_unc,
        d_mean,
        d_unc,
        precision=precision,
        labels=labels,
    )
    headerMap = {Li6: r"\LiS", He4: r"\HeF", He3: r"\HeT", H3: r"\HThree"}
    height = 4
    width = 3
    out = np.empty(
        (height, width), dtype="U50"
    )  # U50 allows up to 50 Unicode characters
    plusMinus = r"\pm"
    meanUncers = [(o2_mean, o2_unc), (o4_mean, o4_unc), (d_mean, d_unc)]
    for i, (mean, unc) in enumerate(meanUncers):
        colStrings = [
            rf"${tmpMean:.3f} {plusMinus} {tmpUnc:.3f} $"
            for tmpMean, tmpUnc in zip(mean, unc)
        ]
        # print("colStrings=", colStrings)
        out[:, i] = np.array(colStrings)
        # print("o2Str=", o2Str)
        # print("out=", out)

    row_labels = [
        r"$F_T\ [\mathrm{fm}^{-1}]$",
        r"$F_L\ [\mathrm{fm}^{-1}]$",
        r"$E^{1N}_{0+}\ [10^{-3}/m_\pi]$",
        r"$L^{1N}_{0+}\ [10^{-3}/m_\pi]$",
    ]
    row_labels = [r"\hline" + row for row in row_labels]
    column_labels = [
        r"$\calO(\delta^2)$",
        r"$\calO(\delta^4)$ ",
        r"Static contribution to $\calO(\delta^4)$",
    ]
    out2 = np.empty((height + 2, width + 1), dtype="U50")
    out2[:, :] = ""
    out2[0, 0] = "$" + headerMap[name] + "$"
    out2[2:, 0] = row_labels
    out2[1, 1:] = column_labels
    out2[2:, 1:] = out
    strOut = a2l.to_ltx(out2, arraytype="tabular", print_out=False)
    strOut = strOut.replace("begin{tabular}", "begin{tabular}{l|l|l|l}")
    # strOut = strOut.replace(r"\\", r"\hline\\")
    print(strOut)


def spread(arr):
    # print("np.shape(arr)=", np.shape(arr))
    # print("np.shape(arr[0,:])=", np.shape(arr[0, :]))
    length = len(arr[0, :])
    out = np.zeros(length)
    for i in range(length):
        subArr = arr[:, i]
        out[i] = (np.max(subArr) - np.min(subArr)) / 2

    return out


def getDiffs(Odelta2Path, returnAll=False, func=ffFromPath):
    """
    given the path Odelta2, returns the static contribution
    """
    diffs = diffFFS(Odelta2Path, func=func)
    # Odelta2+staticContribution=Odelta4 total
    # staticContribution=diff
    diff = diffs[1] - diffs[0]
    if not returnAll:
        return diff
    # else:
    #     return diff, diffs

    else:
        o2 = diffs[0]
        o4 = diffs[1]
        return diff, o2, o4


def diffFFS(path, func=ffFromPath):
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
    if "twobody" in matchString:
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
    else:
        raise ValueError(f"Wrong file type silly \n file={file}")
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
        ffs.append(func(f))

    ffs = np.array(ffs)
    return ffs


def ffFromMat(mat, oper, divide=True):
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
    if divide:
        with np.errstate(divide="ignore", invalid="ignore"):
            tmpMat = (mat / oper).real
    else:
        tmpMat = mat.real

    formfacts = np.zeros(3, dtype=complex)
    for i in range(3):
        tmp = tmpMat[i]
        aves = [
            f for f in tmp.flatten() if (not (isinf(f) or isnan(f))) and abs(f) > 0.0001
        ]

        # Handle empty aves case
        if len(aves) == 0:
            formfacts[i] = 0.0
        else:
            formfacts[i] = np.mean(np.array(aves))

    formfacts = formfacts.real
    formfacts = np.array([(formfacts[0] + formfacts[1]) / 2, formfacts[2]])
    # multiply by -1 in order to match conventions from literature
    # factor of -1 on the matrix element doesn't affect physical observables
    return formfacts * -1


def createOneBodTable(data):
    """
    Print a nice terminal table of form-factor statistics.
    """

    labels = [
        "F_T^{S+V}",
        "F_L^{S+V}",
        "F_T^{S-V}",
        "F_L^{S-V}",
        "E_{0+}^{1N}",
        "L_{0+}^{1N}",
    ]
    COLW = 25  # width of each numeric column
    COl_DIF = 35
    PREC = 3  # number of decimals
    lWidth = 54

    def fmt_pm(mu, sigma, width):
        s = f"{mu: .{PREC}f} ± {sigma:.{PREC}f}"
        return f"{s:^{width}}"  # CENTERED

    # labels = [element + "  [fm^-1]" for element in labels]
    labels[4] += "  [10^-3/m_π]"
    labels[5] += "  [10^-3/m_π]"
    labelWidth = np.max([len(element) for element in labels])

    print("-" * lWidth)
    print(f"{'One Body Results':<{labelWidth}} | {' (mean ± σ)':^{COLW}} | ")
    print("-" * lWidth)

    for i in range(len(data)):
        if i == 2:
            print("-" * lWidth)
        lab = f"{labels[i]}"
        print(f"{lab:<{labelWidth}} | {fmt_pm(data[i, 0], data[i, 1], COLW)} | ")
    print("-" * lWidth)
    print("")


def print_ff_table(
    o2_mean, o2_unc, o4_mean, o4_unc, d_mean, d_unc, precision=3, labels=None
):
    """
    Print a nice terminal table of form-factor statistics.
    """

    COLW = 25  # width of each numeric column
    COl_DIF = 35
    PREC = precision  # number of decimals
    lWidth = 115

    def fmt_pm(mu, sigma, width):
        s = f"{mu: .{PREC}f} ± {sigma:.{PREC}f}"
        return f"{s:^{width}}"  # CENTERED

    # labels = [element + "  [fm^-1]" for element in labels]
    labels[0] += "   [fm^-1]"
    labels[1] += "   [fm^-1]"
    labels[2] += "  [10^-3/m_π]"
    labels[3] += "  [10^-3/m_π]"
    labelWidth = np.max([len(element) for element in labels])
    print("-" * lWidth)
    print(
        f"{'Two Body Results':<{labelWidth}} | "
        f"{'Oδ2 (mean ± σ)':^{COLW}} | "  # format string character ^ is center align, < is left align, > is right align
        f"{'Oδ4 (mean ± σ)':^{COLW}} | "
        f"{('Diff (mean ± σ),  Oδ2 + Diff = Oδ4'):^{COl_DIF}}"
    )
    print("-" * lWidth)

    if labels is None:
        labels = ["T", "L"]
    for i in range(len(o2_mean)):
        if i == 2:
            print("-" * lWidth)
        lab = f"{labels[i]}"
        print(
            f"{lab:<{labelWidth}} | "
            f"{fmt_pm(o2_mean[i], o2_unc[i], COLW)} | "
            f"{fmt_pm(o4_mean[i], o4_unc[i], COLW)} | "
            f"{fmt_pm(d_mean[i], d_unc[i], COl_DIF)}"
        )

    print("-" * lWidth)
    print("")


He3 = "Helium 3"
H3 = "Hydrogen 3"
He4 = "Helium 4"
Li6 = "Lithium 6"


def TwoBods():
    He3Folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/"

    H3Folder = "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod/"

    He4Folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/"

    Li6Folder = "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/"

    params = [(He3Folder, He3), (H3Folder, H3), (He4Folder, He4), (Li6Folder, Li6)]
    print("2 Body Threshold Pion Photoproduction Results")
    for folder, name in params:
        TwoBodFFs(folder, name)


def OneBods():
    He3Folder = (
        r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/thresh/132MeV/"
    )
    # params = [(He3Folder, He3), (H3Folder, H3), (He4Folder, He4), (Li6Folder, Li6)]
    params = [(He3Folder, He3)]

    for folder, name in params:
        OneBodFFs(folder, name)


if __name__ == "__main__":
    # TwoBods()
    OneBods()
