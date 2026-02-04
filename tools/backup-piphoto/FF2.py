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
import re


# ============================================================================
# Core extraction functions
# ============================================================================


def getSpinOperFromFile(path):
    """Extract spin operator matrices from SpinVector entries in output file.

    The output files now contain SpinVector entries that represent the spin
    operator matrices for the nucleus. This function parses these entries
    and converts them to the matrix format used for form factor calculations.

    Parameters
    ----------
    path : str
        Path to the output file containing SpinVector entries.

    Returns
    -------
    numpy.ndarray
        Spin operator matrices with shape (3, N, N) where N is the number
        of spin states for the nucleus.
    """
    with open(path, 'r') as f:
        contents = f.read()

    # Extract nucleus name to determine dimensions
    nucName = rd.getNucName(path)

    # Extract and parse SpinVector block
    block = rd.getBlock(contents, "SpinVector")
    parsed_data = rd.parseBlock(block)

    # Convert to matrix format
    oper = rd.vals2matrix(nucName, parsed_data)

    return oper


def extractMatrixFromFile(path, kind="ScatMat"):
    """Extract matrix elements of a given kind from output file.

    Parameters
    ----------
    path : str
        Path to the output file.
    kind : str
        Type of matrix to extract. Options: "ScatMat", "F^{S+V}", "F^{S-V}",
        "FormFactors", etc.

    Returns
    -------
    numpy.ndarray
        Matrix elements with shape (3, N, N).
    """
    matOut = rd.getQuantNums(path, kind=kind)
    return matOut["MatVals"]


def extractScalarFromFile(path, pattern, label):
    """Extract a scalar value from output file using regex pattern.

    Parameters
    ----------
    path : str
        Path to the output file.
    pattern : str
        Regex pattern to match the value.
    label : str
        Label for error messages.

    Returns
    -------
    float
        Extracted value.
    """
    with open(path, "r") as f:
        text = f.read()

    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"Could not find {label} in file {path}")

    return float(match.group(1))


# ============================================================================
# Form factor calculation
# ============================================================================


def calculateFormFactors(mat, oper, divide=True):
    """Calculate form factors from matrix elements and spin operators.

    Given the output matrix and the spin operator, this function computes
    the form factors by averaging over all non-zero matrix elements weighted
    by the corresponding spin operator elements.

    Parameters
    ----------
    mat : numpy.ndarray
        Matrix of nuclear matrix elements with shape (3, N, N).
    oper : numpy.ndarray
        Spin operator matrices with shape (3, N, N).
    divide : bool
        If True, divide matrix elements by operator elements before averaging.

    Returns
    -------
    numpy.ndarray
        Array of length 2 containing [F_T, F_L] where F_T is the average
        of x and y components and F_L is the z component.

    Notes
    -----
    The function returns the real part and multiplies by -1 to match
    literature conventions. This factor doesn't affect physical observables.
    """
    # Divide matrix elements by operator elements if requested
    if divide:
        with np.errstate(divide='ignore', invalid='ignore'):
            tmpMat = (mat / oper).real
    else:
        tmpMat = mat.real

    # Average over non-zero, finite elements for each component
    formfacts = np.zeros(3)
    for i in range(3):
        tmp = tmpMat[i].flatten()
        # Filter out inf, nan, and near-zero values
        valid = [f for f in tmp if not (isinf(f) or isnan(f)) and abs(f) > 0.0001]
        formfacts[i] = np.mean(valid) if valid else 0.0

    # Combine x and y components into transverse component
    f_transverse = (formfacts[0] + formfacts[1]) / 2
    f_longitudinal = formfacts[2]

    # Multiply by -1 to match literature conventions
    return np.array([f_transverse, f_longitudinal]) * -1


def getFormFactorsFromFile(path, kind="FormFactors", divide=True):
    """Extract form factors from an output file.

    Parameters
    ----------
    path : str
        Path to the output file.
    kind : str
        Type of matrix to extract for form factor calculation.
    divide : bool
        Whether to divide by spin operator elements.

    Returns
    -------
    numpy.ndarray
        Array containing form factors.
    """
    mat = extractMatrixFromFile(path, kind=kind)
    oper = getSpinOperFromFile(path)
    return calculateFormFactors(mat, oper, divide=divide)


def getScatMatFromFile(path):
    """Extract E_0+^1N and L_0+^1N directly from output file.

    These values are pre-computed and printed in the output file as:
    E_0+^1N=  1.751478 (averaged)
    L_0+^1N= -1.932005

    Parameters
    ----------
    path : str
        Path to the output file.

    Returns
    -------
    numpy.ndarray
        Array [E_0+^1N, L_0+^1N].
    """
    E0p = extractScalarFromFile(path, r"E_0\+\^1N=\s*([+-]?\d+\.?\d*)", "E_0+^1N")
    L0p = extractScalarFromFile(path, r"L_0\+\^1N=\s*([+-]?\d+\.?\d*)", "L_0+^1N")
    return np.array([E0p, L0p])


def getCombinedResults(path):
    """Extract both form factors and scattering amplitudes from output file.

    Parameters
    ----------
    path : str
        Path to the output file.

    Returns
    -------
    numpy.ndarray
        Array containing [F_T, F_L, E_0+, L_0+].
    """
    ffs = getFormFactorsFromFile(path)
    scatMat = getScatMatFromFile(path)
    return np.concatenate([ffs, scatMat])


# ============================================================================
# Two-body calculations (Odelta2 vs Odelta4)
# ============================================================================


def findMatchingFile(path, swap_order):
    """Find the matching Odelta2 or Odelta4 file for a given path.

    Parameters
    ----------
    path : str
        Path to a file containing either 'Odelta2' or 'Odelta4'.
    swap_order : bool
        If True, find Odelta4 for given Odelta2 (and vice versa).

    Returns
    -------
    str
        Path to the matching file.

    Raises
    ------
    ValueError
        If the file doesn't contain order information.
    AssertionError
        If no unique match is found.
    """
    file = basename(path)
    folder = dirname(path)

    if "twobody" not in file:
        raise ValueError(f"Wrong file type - expected twobody file, got: {file}")

    # Determine which order to swap to
    if "Odelta2" in file:
        old_str = "Odelta2"
        new_str = "Odelta4"
        path_is_o2 = True
    elif "Odelta4" in file:
        old_str = "Odelta4"
        new_str = "Odelta2"
        path_is_o2 = False
    else:
        raise ValueError(f"File doesn't contain Odelta2 or Odelta4: {file}")

    # Find matching file
    match_file = file.replace(old_str, new_str)
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    matches = [f for f in files if match_file in f]

    assert len(matches) == 1, f"Expected 1 match, found {len(matches)}:\n" + "\n".join(matches)

    return join(folder, matches[0])


def getOrderedPaths(path):
    """Get paths for both Odelta2 and Odelta4 files in that order.

    Parameters
    ----------
    path : str
        Path to either an Odelta2 or Odelta4 file.

    Returns
    -------
    tuple
        (path_o2, path_o4) - paths to Odelta2 and Odelta4 files.
    """
    file = basename(path)

    if "Odelta2" in file:
        path_o2 = path
        path_o4 = findMatchingFile(path, swap_order=True)
    else:
        path_o2 = findMatchingFile(path, swap_order=True)
        path_o4 = path

    return path_o2, path_o4


def calculateStaticContribution(path, extraction_func=getFormFactorsFromFile):
    """Calculate static contribution to form factors.

    The static contribution is defined as: Odelta4 - Odelta2

    Parameters
    ----------
    path : str
        Path to either Odelta2 or Odelta4 file.
    extraction_func : callable
        Function to extract values from files.

    Returns
    -------
    tuple
        (static_contrib, o2_result, o4_result)
    """
    path_o2, path_o4 = getOrderedPaths(path)

    o2_result = extraction_func(path_o2)
    o4_result = extraction_func(path_o4)
    static_contrib = o4_result - o2_result

    return static_contrib, o2_result, o4_result


# ============================================================================
# Statistical analysis
# ============================================================================


def calculateSpread(values):
    """Calculate spread (half of range) for each column.

    Parameters
    ----------
    values : numpy.ndarray
        2D array where each row is a sample and each column is a quantity.

    Returns
    -------
    numpy.ndarray
        1D array of spreads (max - min) / 2 for each column.
    """
    return (np.max(values, axis=0) - np.min(values, axis=0)) / 2


def analyzeFiles(files, extraction_func, return_individuals=False):
    """Analyze multiple files and compute statistics.

    Parameters
    ----------
    files : list
        List of file paths to analyze.
    extraction_func : callable
        Function to extract values from each file.
    return_individuals : bool
        If True, also return individual values.

    Returns
    -------
    tuple or dict
        If return_individuals=False: (mean, spread)
        If return_individuals=True: dict with 'mean', 'spread', 'values'
    """
    values = np.array([extraction_func(f) for f in files])
    mean = np.mean(values, axis=0)
    spread = calculateSpread(values)

    if return_individuals:
        return {'mean': mean, 'spread': spread, 'values': values}
    return mean, spread


def analyzeTwoBodyFiles(files, extraction_func=getCombinedResults):
    """Analyze two-body files for Odelta2, Odelta4, and static contributions.

    Parameters
    ----------
    files : list
        List of Odelta2 file paths.
    extraction_func : callable
        Function to extract values from files.

    Returns
    -------
    dict
        Dictionary with keys 'o2', 'o4', 'static', each containing 'mean' and 'spread'.
    """
    o2_vals = []
    o4_vals = []
    static_vals = []

    for path in files:
        static, o2, o4 = calculateStaticContribution(path, extraction_func)
        o2_vals.append(o2)
        o4_vals.append(o4)
        static_vals.append(static)

    o2_vals = np.array(o2_vals)
    o4_vals = np.array(o4_vals)
    static_vals = np.array(static_vals)

    return {
        'o2': {'mean': np.mean(o2_vals, axis=0), 'spread': calculateSpread(o2_vals)},
        'o4': {'mean': np.mean(o4_vals, axis=0), 'spread': calculateSpread(o4_vals)},
        'static': {'mean': np.mean(static_vals, axis=0), 'spread': calculateSpread(static_vals)}
    }


def analyzeOneBodyFiles(files):
    """Analyze one-body files for multiple form factor types.

    Parameters
    ----------
    files : list
        List of output file paths.

    Returns
    -------
    numpy.ndarray
        Array of shape (6, 2) containing [mean, spread] for each of:
        F_T^{S+V}, F_L^{S+V}, F_T^{S-V}, F_L^{S-V}, E_0+^1N, L_0+^1N
    """
    kinds = ["F^{S+V}", "F^{S-V}", "ScatMat"]
    results = []

    for kind in kinds:
        values = []
        for path in files:
            if kind == "ScatMat":
                val = getScatMatFromFile(path)
            else:
                val = getFormFactorsFromFile(path, kind=kind, divide=False) * -1
            values.append(val)

        values = np.array(values)
        mean = np.mean(values, axis=0)
        spread = calculateSpread(values)

        # Flatten into [mean1, spread1, mean2, spread2] format
        for m, s in zip(mean, spread):
            results.append([m, s])

    return np.array(results)


# ============================================================================
# Output formatting
# ============================================================================


def printOneBodyTable(data):
    """Print formatted table of one-body form factor results.

    Parameters
    ----------
    data : numpy.ndarray
        Array of shape (6, 2) with [mean, spread] for each quantity.
    """
    labels = [
        "F_T^{S+V}",
        "F_L^{S+V}",
        "F_T^{S-V}",
        "F_L^{S-V}",
        "E_{0+}^{1N}  [10^-3/m_π]",
        "L_{0+}^{1N}  [10^-3/m_π]",
    ]

    col_width = 25
    prec = 3
    label_width = max(len(label) for label in labels)
    total_width = label_width + col_width + 10

    def format_val(mean, spread):
        s = f"{mean: .{prec}f} ± {spread:.{prec}f}"
        return f"{s:^{col_width}}"

    print("-" * total_width)
    print(f"{'One Body Results':<{label_width}} | {' (mean ± σ)':^{col_width}} | ")
    print("-" * total_width)

    for i, label in enumerate(labels):
        if i == 2 or i == 4:
            print("-" * total_width)
        print(f"{label:<{label_width}} | {format_val(data[i, 0], data[i, 1])} | ")

    print("-" * total_width)
    print()


def printTwoBodyTable(results, labels, precision=3):
    """Print formatted table of two-body results.

    Parameters
    ----------
    results : dict
        Dictionary with 'o2', 'o4', 'static' keys.
    labels : list
        Labels for each row.
    precision : int
        Number of decimal places.
    """
    col_width = 25
    col_diff = 35
    label_width = max(len(label + "   [fm^-1]") for label in labels[:2])
    total_width = label_width + 2 * col_width + col_diff + 15

    def format_val(mean, spread, width):
        s = f"{mean: .{precision}f} ± {spread:.{precision}f}"
        return f"{s:^{width}}"

    # Add units to labels
    full_labels = [
        labels[0] + "   [fm^-1]",
        labels[1] + "   [fm^-1]",
        labels[2] + "  [10^-3/m_π]",
        labels[3] + "  [10^-3/m_π]",
    ]

    print("-" * total_width)
    print(
        f"{'Two Body Results':<{label_width}} | "
        f"{'Oδ2 (mean ± σ)':^{col_width}} | "
        f"{'Oδ4 (mean ± σ)':^{col_width}} | "
        f"{'Diff (mean ± σ),  Oδ2 + Diff = Oδ4':^{col_diff}}"
    )
    print("-" * total_width)

    o2 = results['o2']
    o4 = results['o4']
    static = results['static']

    for i, label in enumerate(full_labels):
        if i == 2:
            print("-" * total_width)
        print(
            f"{label:<{label_width}} | "
            f"{format_val(o2['mean'][i], o2['spread'][i], col_width)} | "
            f"{format_val(o4['mean'][i], o4['spread'][i], col_width)} | "
            f"{format_val(static['mean'][i], static['spread'][i], col_diff)}"
        )

    print("-" * total_width)
    print()


def createLatexTable(results, nuc_name):
    """Create LaTeX table for two-body results.

    Parameters
    ----------
    results : dict
        Dictionary with 'o2', 'o4', 'static' keys.
    nuc_name : str
        Nucleus name for header.

    Returns
    -------
    str
        LaTeX table code.
    """
    header_map = {
        'Lithium 6': r"\LiS",
        'Helium 4': r"\HeF",
        'Helium 3': r"\HeT",
        'Hydrogen 3': r"\HThree"
    }

    height = 4
    width = 3
    out = np.empty((height, width), dtype='U50')
    plus_minus = r"\pm"

    o2 = results['o2']
    o4 = results['o4']
    static = results['static']

    mean_uncs = [
        (o2['mean'], o2['spread']),
        (o4['mean'], o4['spread']),
        (static['mean'], static['spread'])
    ]

    for i, (mean, unc) in enumerate(mean_uncs):
        col_strings = [
            rf"${m:.3f} {plus_minus} {u:.3f} $"
            for m, u in zip(mean, unc)
        ]
        out[:, i] = np.array(col_strings)

    row_labels = [
        r"\hline$F_T\ [\mathrm{fm}^{-1}]$",
        r"\hline$F_L\ [\mathrm{fm}^{-1}]$",
        r"\hline$E^{1N}_{0+}\ [10^{-3}/m_\pi]$",
        r"\hline$L^{1N}_{0+}\ [10^{-3}/m_\pi]$",
    ]

    column_labels = [
        r"$\calO(\delta^2)$",
        r"$\calO(\delta^4)$ ",
        r"Static contribution to $\calO(\delta^4)$",
    ]

    out2 = np.empty((height + 2, width + 1), dtype='U50')
    out2[:, :] = ""
    out2[0, 0] = "$" + header_map.get(nuc_name, nuc_name) + "$"
    out2[2:, 0] = row_labels
    out2[1, 1:] = column_labels
    out2[2:, 1:] = out

    latex_str = a2l.to_ltx(out2, arraytype="tabular", print_out=False)
    latex_str = latex_str.replace("begin{tabular}", "begin{tabular}{l|l|l|l}")

    return latex_str


# ============================================================================
# High-level interface functions
# ============================================================================


def processOneBodyFolder(folder, name):
    """Process all one-body files in a folder and print results.

    Parameters
    ----------
    folder : str
        Path to folder containing output files.
    name : str
        Name of the nucleus for display.
    """
    files = [join(folder, f) for f in listdir(folder) if isfile(join(folder, f))]

    print(f"\n{name} One Body Results")
    data = analyzeOneBodyFiles(files)
    printOneBodyTable(data)


def processTwoBodyFolder(folder, name):
    """Process all two-body files in a folder and print results.

    Parameters
    ----------
    folder : str
        Path to folder containing output files.
    name : str
        Name of the nucleus for display.
    """
    all_files = [join(folder, f) for f in listdir(folder) if isfile(join(folder, f))]
    files = [f for f in all_files if "Odelta2" in f and "j12max=1" in f]

    print(f"\n{name} Two Body Results")
    labels = ["F_T", "F_L", "E_0+", "L_0+"]

    results = analyzeTwoBodyFiles(files, extraction_func=getCombinedResults)
    printTwoBodyTable(results, labels, precision=6)

    # Print LaTeX table
    latex_table = createLatexTable(results, name)
    print(latex_table)


# ============================================================================
# Main execution
# ============================================================================

# Nucleus names
He3 = "Helium 3"
H3 = "Hydrogen 3"
He4 = "Helium 4"
Li6 = "Lithium 6"


def runTwoBodyAnalysis():
    """Run analysis for all two-body calculations."""
    folders = {
        He3: r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/",
        H3: "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod/",
        He4: r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/",
        Li6: "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/",
    }

    print("2 Body Threshold Pion Photoproduction Results")
    for name, folder in folders.items():
        processTwoBodyFolder(folder, name)


def runOneBodyAnalysis():
    """Run analysis for all one-body calculations."""
    folders = {
        He3: r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/thresh/132MeV/",
    }

    for name, folder in folders.items():
        processOneBodyFolder(folder, name)


if __name__ == "__main__":
    runTwoBodyAnalysis()
    runOneBodyAnalysis()
