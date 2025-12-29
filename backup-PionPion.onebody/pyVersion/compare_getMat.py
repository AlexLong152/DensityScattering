#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare getMmat (Python) and getMat (Fortran) functions

This script tests whether the Python implementation of getMmat matches
the Fortran implementation of getMat.

@author: test script
"""

import numpy as np
import sys
import os
import subprocess

# Import the Python implementation
from PionScatLib import getMmat, getCS as getCS_pl


def ensure_fortran_executable():
    """
    Ensure the Fortran testGetMat executable exists. Compile it if necessary.
    """
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fortran_exe = os.path.join(parent_dir, "testGetMat")

    if not os.path.exists(fortran_exe):
        print("testGetMat executable not found. Compiling...")

        # Change to parent directory and compile
        compile_commands = [
            f"cd {parent_dir}",
            "rm -f testGetMat testGetMat.o PionScat.o parseFile.o utility_suite.o",
            "gfortran -c -ffixed-form -ffixed-line-length-none parseFile.f",
            "gfortran -c -ffixed-form -ffixed-line-length-none utility_suite.f",
            "gfortran -c -ffixed-form -ffixed-line-length-none PionScat.f",
            "gfortran -c -ffixed-form -ffixed-line-length-none testGetMat.f",
            "gfortran -o testGetMat testGetMat.o PionScat.o utility_suite.o parseFile.o",
        ]

        compile_cmd = " && ".join(compile_commands)
        result = os.system(compile_cmd)

        if result != 0:
            raise RuntimeError(f"Failed to compile testGetMat. Exit code: {result}")

        print("Compilation successful!")

    return fortran_exe


def ensure_fortran_cs_executable():
    """
    Ensure the Fortran testGetCS executable exists. Compile it if necessary.
    """
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fortran_exe = os.path.join(parent_dir, "testGetCS")

    if not os.path.exists(fortran_exe):
        print("testGetCS executable not found. Compiling...")

        # Change to parent directory and compile
        compile_commands = [
            f"cd {parent_dir}",
            "rm -f testGetCS testGetCS.o",
            "gfortran -c -ffixed-form -ffixed-line-length-none parseFile.f",
            "gfortran -c -ffixed-form -ffixed-line-length-none utility_suite.f",
            "gfortran -c -ffixed-form -ffixed-line-length-none PionScat.f",
            "gfortran -c -ffixed-form -ffixed-line-length-none testGetCS.f",
            "gfortran -o testGetCS testGetCS.o PionScat.o utility_suite.o parseFile.o",
        ]

        compile_cmd = " && ".join(compile_commands)
        result = os.system(compile_cmd)

        if result != 0:
            raise RuntimeError(f"Failed to compile testGetCS. Exit code: {result}")

        print("Compilation successful!")

    return fortran_exe


def call_fortran_getMat(sqrtS, theta_deg, isospin, piCharge, sqrtSReal=None):
    """
    Call the Fortran testGetMat executable to get the matrix result.

    Parameters
    ----------
    sqrtS : float
        Center-of-mass energy in MeV
    theta_deg : float
        Scattering angle in degrees
    isospin : int
        Nucleon isospin (1=proton, -1=neutron)
    piCharge : int
        Pion charge (-1, 0, +1)
    sqrtSReal : float, optional
        Real sqrt(S) value (defaults to sqrtS if not provided)

    Returns
    -------
    numpy.ndarray
        2x2 complex matrix in Python indexing
    """
    if sqrtSReal is None:
        sqrtSReal = sqrtS

    # Path to the Fortran executable
    fortran_exe = os.path.join(os.path.dirname(__file__), "..", "testGetMat")

    # Prepare input string
    input_str = f"{sqrtS} {theta_deg} {isospin} {piCharge}\n"

    # Call the Fortran program
    try:
        result = subprocess.run(
            [fortran_exe], input=input_str, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error calling Fortran program: {e}")
        print(f"Stderr: {e.stderr}")
        raise

    # Parse the output
    # Format: mat(-1,-1)= real imag
    lines = result.stdout.strip().split("\n")

    mat_fortran = {}
    for line in lines:
        if "mat" in line:
            # Parse line like: mat(-1,-1)=  0.512457124244193E+00   0.331526143042571E+00
            parts = line.split("=")
            key = parts[0].strip()  # e.g., "mat(-1,-1)"
            values = parts[1].strip().split()
            real_part = float(values[0])
            imag_part = float(values[1])

            # Extract indices from key like "mat(-1,-1)"
            indices = key.replace("mat(", "").replace(")", "").split(",")
            i = int(indices[0].strip())
            j = int(indices[1].strip())

            mat_fortran[(i, j)] = complex(real_part, imag_part)

    # Convert to Python 2x2 matrix (0-based indexing)
    mat = np.zeros((2, 2), dtype=complex)
    mat[0, 0] = mat_fortran[(-1, -1)]
    mat[0, 1] = mat_fortran[(-1, 1)]
    mat[1, 0] = mat_fortran[(1, -1)]
    mat[1, 1] = mat_fortran[(1, 1)]

    return mat


def call_fortran_getCS(sqrtS, theta_deg, isospin, piCharge):
    """
    Call the Fortran testGetCS executable to get the cross section.

    Parameters
    ----------
    sqrtS : float
        Center-of-mass energy in MeV
    theta_deg : float
        Scattering angle in degrees
    isospin : int
        Nucleon isospin (1=proton, -1=neutron)
    piCharge : int
        Pion charge (-1, 0, +1)

    Returns
    -------
    float
        Cross section in millibarns
    """
    # Path to the Fortran executable
    fortran_exe = os.path.join(os.path.dirname(__file__), "..", "testGetCS")

    # Prepare input string
    input_str = f"{sqrtS} {theta_deg} {isospin} {piCharge}\n"

    # Call the Fortran program
    try:
        result = subprocess.run(
            [fortran_exe], input=input_str, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error calling Fortran program: {e}")
        print(f"Stderr: {e.stderr}")
        raise

    # Parse the output
    # Format: CrossSec= value
    lines = result.stdout.strip().split("\n")

    for line in lines:
        if "CrossSec=" in line:
            parts = line.split("=")
            value_str = parts[1].strip()
            # Fortran scientific notation format doesn't always include 'E'
            # e.g., "-0.396659013727452-153" should be "-0.396659013727452E-153"
            # We need to handle cases like "1.23456789012345+100" or "-1.23-100"
            # Look for + or - signs after the first character (and after any initial sign)
            import re

            # Pattern: optional sign, digits and dot, then +/- followed by more digits
            value_str = re.sub(r"([0-9])([+-][0-9])", r"\1E\2", value_str)
            cs_value = float(value_str)
            return cs_value

    raise RuntimeError("Could not parse CrossSec from Fortran output")


def compare_cross_sections(cs1, cs2, label1="CS1", label2="CS2", tol=1e-3):
    """
    Compare two cross section values and report differences.

    Parameters
    ----------
    cs1, cs2 : float
        Cross sections to compare
    label1, label2 : str
        Labels for the cross sections
    tol : float
        Tolerance for comparison

    Returns
    -------
    bool
        True if cross sections match within tolerance
    """
    print(f"\n{label1}: {cs1:.10e} mb")
    print(f"{label2}: {cs2:.10e} mb")

    # Calculate difference
    diff = abs(cs1 - cs2)
    print(f"\nAbsolute difference: {diff:.10e} mb")

    if abs(cs2) > 1e-10:
        rel_diff = diff / abs(cs2)
        print(f"Relative difference: {rel_diff:.10e}")

    if diff < tol:
        print("✓ Cross sections MATCH within tolerance")
        return True
    else:
        print("✗ Cross sections DO NOT MATCH")
        return False


def compute_cs_from_matrix(mat):
    """
    Compute cross section from a 2x2 scattering matrix.

    Calculates: CrossSec = 10.0 * Real(Tr(mat * mat^dagger)) / 4.0 * 2.0

    Parameters
    ----------
    mat : numpy.ndarray
        2x2 complex scattering matrix

    Returns
    -------
    float
        Cross section in millibarns
    """
    # Calculate Hermitian conjugate
    matDag = np.conjugate(mat.T)

    # Calculate trace of mat * matDag
    trace = 0.0 + 0.0j
    for i in [0, 1]:  # Spin indices
        for j in [0, 1]:
            trace += mat[i, j] * matDag[j, i]

    # Apply formula
    CrossSec = 10.0 * np.real(trace) / 4.0
    fudgeFactor = 2.0
    CrossSec = CrossSec * fudgeFactor

    return CrossSec


def compare_matrices(mat1, mat2, label1="Mat1", label2="Mat2", tol=1e-3):
    """
    Compare two 2x2 matrices and report differences.

    Parameters
    ----------
    mat1, mat2 : numpy.ndarray
        2x2 matrices to compare
    label1, label2 : str
        Labels for the matrices
    tol : float
        Tolerance for comparison

    Returns
    -------
    bool
        True if matrices match within tolerance
    """
    print(f"\n{label1}:")
    for i in range(2):
        for j in range(2):
            # Map from Python indexing to spin indexing
            spin_i = 1 if i == 1 else -1
            spin_j = 1 if j == 1 else -1
            print(f"  [{spin_i:2d},{spin_j:2d}] = {mat1[i, j]:.10e}")

    print(f"\n{label2}:")
    for i in range(2):
        for j in range(2):
            spin_i = 1 if i == 1 else -1
            spin_j = 1 if j == 1 else -1
            print(f"  [{spin_i:2d},{spin_j:2d}] = {mat2[i, j]:.10e}")

    # Calculate differences
    diff = np.abs(mat1 - mat2)
    max_diff = np.max(diff)

    print(f"\nMaximum absolute difference: {max_diff:.10e}")

    # Compute and display cross sections from matrices
    cs1 = compute_cs_from_matrix(mat1)
    cs2 = compute_cs_from_matrix(mat2)
    print(f"\nCross sections from matrices:")
    print(f"  {label1}: {cs1:.10e} mb")
    print(f"  {label2}: {cs2:.10e} mb")
    print(f"  Difference:   {abs(cs1 - cs2):.10e} mb")

    if max_diff < tol:
        print("\n✓ Matrices MATCH within tolerance")
        return True
    else:
        print("\n✗ Matrices DO NOT MATCH")
        print("\nElement-wise absolute differences:")
        for i in range(2):
            for j in range(2):
                spin_i = 1 if i == 1 else -1
                spin_j = 1 if j == 1 else -1
                print(f"  [{spin_i:2d},{spin_j:2d}] diff = {diff[i, j]:.10e}")
        print("\nRelative differences (where |fortran| > 1e-10):")
        for i in range(2):
            for j in range(2):
                if abs(mat2[i, j]) > 1e-10:
                    rel_diff = diff[i, j] / abs(mat2[i, j])
                    spin_i = 1 if i == 1 else -1
                    spin_j = 1 if j == 1 else -1
                    print(f"  [{spin_i:2d},{spin_j:2d}] rel diff = {rel_diff:.10e}")
        return False


def test_matrix_single_case(sqrtS, theta_deg, isospin, piCharge):
    """
    Test matrix elements for a single case.

    Parameters
    ----------
    sqrtS : float
        Center-of-mass energy in MeV
    theta_deg : float
        Scattering angle in degrees
    isospin : int
        Nucleon isospin (1=proton, -1=neutron)
    piCharge : int
        Pion charge (-1, 0, +1)

    Returns
    -------
    tuple
        (py_mat, fort_mat, match)
    """
    x = np.cos(theta_deg * np.pi / 180)

    print("=" * 70)
    print(
        f"Test case: sqrtS={sqrtS} MeV, theta={theta_deg}°, "
        f"isospin={isospin}, piCharge={piCharge}"
    )
    print("=" * 70)

    py_mat = getMmat(sqrtS, x, isospin, piCharge)
    fort_mat = call_fortran_getMat(sqrtS, theta_deg, isospin, piCharge, sqrtS)
    mat_match = compare_matrices(
        py_mat, fort_mat, label1="Python getMmat", label2="Fortran getMat"
    )

    return py_mat, fort_mat, mat_match


def test_cs_single_case(sqrtS, theta_deg, isospin, piCharge):
    """
    Test cross sections for a single case.

    Parameters
    ----------
    sqrtS : float
        Center-of-mass energy in MeV
    theta_deg : float
        Scattering angle in degrees
    isospin : int
        Nucleon isospin (1=proton, -1=neutron)
    piCharge : int
        Pion charge (-1, 0, +1)

    Returns
    -------
    bool
        True if cross sections match within tolerance
    """
    x = np.cos(theta_deg * np.pi / 180)

    print("=" * 70)
    print(
        f"Test case: sqrtS={sqrtS} MeV, theta={theta_deg}°, "
        f"isospin={isospin}, piCharge={piCharge}"
    )
    print("=" * 70)

    py_cs = getCS_pl(sqrtS, x, isospin, piCharge)
    fort_cs = call_fortran_getCS(sqrtS, theta_deg, isospin, piCharge)
    cs_match = compare_cross_sections(
        py_cs, fort_cs, label1="Python getCS", label2="Fortran getCS", tol=0.01
    )

    return cs_match


def main():
    """
    Run comparison tests for several test cases.
    """
    # Ensure the Fortran executables exist
    ensure_fortran_executable()
    ensure_fortran_cs_executable()

    print("\n" + "=" * 70)
    print("Comparing Python vs Fortran Implementations")
    print("  - Matrix elements: getMmat() vs getMat()")
    print("  - Cross sections: getCS() vs getCS()")
    print("=" * 70)

    # Test cases: (sqrtS, theta_deg, isospin, piCharge)
    test_cases = [
        (1162, 30, 1, 1),  # π+ on proton at 30°
        (1162, 90, 1, -1),  # π- on proton at 90°
        (1298, 45, 1, 0),  # π0 on proton at 45°
        (1200, 60, -1, 1),  # π+ on neutron at 60°
    ]

    # Phase 1: Test all matrix elements
    print("\n" + "=" * 70)
    print("PHASE 1: MATRIX ELEMENT COMPARISONS")
    print("=" * 70 + "\n")

    mat_results = []
    for sqrtS, theta, iso, charge in test_cases:
        py_mat, fort_mat, mat_match = test_matrix_single_case(sqrtS, theta, iso, charge)
        mat_results.append(mat_match)
        print("\n")

    # Phase 2: Test all cross sections
    print("\n" + "=" * 70)
    print("PHASE 2: CROSS SECTION COMPARISONS")
    print("=" * 70 + "\n")

    cs_results = []
    for sqrtS, theta, iso, charge in test_cases:
        cs_match = test_cs_single_case(sqrtS, theta, iso, charge)
        cs_results.append(cs_match)
        print("\n")

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    n_mat_pass = sum(1 for r in mat_results if r)
    n_cs_pass = sum(1 for r in cs_results if r)
    n_total = len(test_cases)

    print(f"Matrix element tests passed: {n_mat_pass}/{n_total}")
    print(f"Cross section tests passed:  {n_cs_pass}/{n_total}")

    if n_mat_pass == n_total and n_cs_pass == n_total:
        print("\n✓ All tests PASSED!")
    else:
        n_failed = (n_total - n_mat_pass) + (n_total - n_cs_pass)
        print(f"\n✗ {n_failed} test(s) FAILED")

    modPath = r"../pionscatlib.mod"
    if os.path.exists(modPath):
        os.remove(modPath)


if __name__ == "__main__":
    sys.exit(main())
