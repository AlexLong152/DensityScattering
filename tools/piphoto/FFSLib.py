# -*- coding: utf-8 -*-
"""Shared library for form factor and scattering amplitude extraction.

Core extraction functions used by both one-body and two-body analyses.

Author: alexl
"""

import numpy as np
from math import isinf, isnan
import readDensity as rd
from os import listdir
from os.path import isfile, join
import re


# ============================================================================
# Constants
# ============================================================================

percent = 2  # assumed percent error from densities

# Nucleus names
He3 = "Helium 3"
H3 = "Hydrogen 3"
He4 = "Helium 4"
Li6 = "Lithium 6"


# ============================================================================
# Core extraction functions
# ============================================================================


def getSpinOperFromFile(path):
    """Extract spin operator matrices from output file.

    For one-body files: reads SpinVec from file.
    For two-body files: uses rd.spin4Nuc() as fallback.
    """
    with open(path, "r") as f:
        contents = f.read()

    nucName = rd.getNucName(path)

    expected = rd.spin4Nuc(nucName)

    # Try to read SpinVec (one-body files)
    if "onebody" in path:
        try:
            block = rd.getBlock(contents, "SpinVec")
            if block.strip():
                parsed_data = rd.parseBlock(block)
                if parsed_data:
                    result = rd.vals2matrix(nucName, parsed_data)
                    passed = np.allclose(result, expected)
                    if not passed:
                        print(f"SpinVec from file doesn't match spin4Nuc for {nucName}")
                        # print("result=\n", result)
                        # print("expected=\n", expected, "\n")
                        for i in range(3):
                            res = result[i]
                            exp = expected[i]
                            if not np.allclose(res, exp):
                                print("i=", i)
                                print(res, "\n")
                                print(exp, "\n")
                                print(res - exp)

                        print("\n\n")
                        # print("result/expected=", result / expected)

                    return result
        except AssertionError:
            raise  # checks for the spin vectors
        except Exception:
            pass

    # Fallback to rd.spin4Nuc() for two-body files
    return expected


def extractMatrixFromFile(path, kind="ScatMat"):
    """Extract matrix elements from output file."""
    matOut = rd.getQuantNums(path, kind=kind)
    return matOut["MatVals"]


def extractScalarFromFile(path, pattern, label):
    """Extract scalar value using regex pattern."""
    with open(path, "r") as f:
        match = re.search(pattern, f.read())

    if not match:
        raise ValueError(f"Could not find {label} in file {path}")

    return float(match.group(1))


def calculateFormFactors(mat, oper, divide=True):
    """Calculate form factors from matrix elements and spin operators.

    Returns [F_T, F_L] where F_T is average of x,y components and F_L is z.
    Multiplies by -1 to match literature conventions.
    """
    if divide:
        with np.errstate(divide="ignore", invalid="ignore"):
            fullMat = mat / oper
            fullMat = np.nan_to_num(fullMat, nan=0.0, posinf=0.0, neginf=0.0)
            tmpMat = fullMat.real
            # assert np.allclose(fullMat, tmpMat, atol=1e-5, rtol=1e-3), (
            #     f"fullMat\n={fullMat},\ntmpMat=\ntmpMat"
            # )
    else:
        tmpMat = mat.real

    formfacts = np.zeros(3)
    for i in range(3):
        flattened = tmpMat[i].flatten()
        valid = [f for f in flattened if not (isinf(f) or isnan(f)) and abs(f) > 0.0001]
        formfacts[i] = np.mean(valid) if valid else 0.0

    # print("formfacts\n=", formfacts)
    # Combine x,y into transverse and return [F_T, F_L]
    return np.array([(formfacts[0] + formfacts[1]) / 2, formfacts[2]]) * -1


def getFormFactorsFromFile(path, kind="FormFactors", divide=True):
    """Extract form factors [F_T, F_L] from output file."""
    mat = extractMatrixFromFile(path, kind=kind)
    oper = getSpinOperFromFile(path)
    return calculateFormFactors(mat, oper, divide=divide)


def getScatMatFromFile(path):
    """Extract E_0+^1N and L_0+^1N from output file.

    One-body: reads directly from text.
    Two-body: reads Matrix section and divides by spin operator.
    """
    if "onebody" in path:
        with open(path, "r") as f:
            contents = f.read()

        e_match = re.search(r"E_0\+\^1N=\s*([-\d.]+)", contents)
        l_match = re.search(r"L_0\+\^1N=\s*([-\d.]+)", contents)

        if not e_match or not l_match:
            raise ValueError(f"Could not find E_0+^1N or L_0+^1N in file {path}")

        return np.array([float(e_match.group(1)), float(l_match.group(1))])
    else:
        assert "twobody" in path
        mat = extractMatrixFromFile(path, kind="Matrix")
        oper = getSpinOperFromFile(path)

        return calculateFormFactors(mat, oper, divide=True)


# ============================================================================
# Statistical analysis (shared)
# ============================================================================


def calculateSpread(values):
    """Calculate spread (max - min) / 2 for each column."""
    return (np.max(values, axis=0) - np.min(values, axis=0)) / 2


# ============================================================================
# Output formatting (shared)
# ============================================================================


def _format_parens(mean, spread, precision=3):
    """Format mean(uncertainty) string for LaTeX, e.g. -0.047(5) or -0.047(13)."""
    scale = 10**precision
    unc = int(round(abs(spread) * scale))
    return rf"${mean:.{precision}f}({unc})$"


def _format_value(mean, spread, width, precision):
    """Format mean ± spread string centered in given width."""
    return f"{mean: .{precision}f} ± {spread:.{precision}f}".center(width)


# ============================================================================
# High-level interface (shared)
# ============================================================================


def _get_files_from_folder(folder):
    """Get all files from a folder."""
    return [join(folder, f) for f in listdir(folder) if isfile(join(folder, f))]
