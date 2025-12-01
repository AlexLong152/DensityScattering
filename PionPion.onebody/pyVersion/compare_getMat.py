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
from PionScatLib import getMmat


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
            "rm -f testGetMat testGetMat.o PionScat.o",
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

    if max_diff < tol:
        print("✓ Matrices MATCH within tolerance")
        return True
    else:
        print("✗ Matrices DO NOT MATCH")
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


def test_single_case(sqrtS, theta_deg, isospin, piCharge):
    """
    Test a single case comparing Python and Fortran implementations.

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
    """
    x = np.cos(theta_deg * np.pi / 180)

    print("=" * 70)
    print(
        f"Test case: sqrtS={sqrtS} MeV, theta={theta_deg}°, "
        f"isospin={isospin}, piCharge={piCharge}"
    )
    print("=" * 70)

    # Get Python result
    py_mat = getMmat(sqrtS, x, isospin, piCharge)

    # Get Fortran result
    fort_mat = call_fortran_getMat(sqrtS, theta_deg, isospin, piCharge, sqrtS)

    # Compare
    match = compare_matrices(
        py_mat, fort_mat, label1="Python getMmat", label2="Fortran getMat"
    )
    return match


def main():
    """
    Run comparison tests for several test cases.
    """
    # Ensure the Fortran executable exists
    ensure_fortran_executable()

    print("\n" + "=" * 70)
    print("Comparing Python getMmat() vs Fortran getMat()")
    print("=" * 70)

    # Test cases: (sqrtS, theta_deg, isospin, piCharge)
    test_cases = [
        (1162, 30, 1, 1),  # π+ on proton at 30°
        (1162, 90, 1, -1),  # π- on proton at 90°
        (1298, 45, 1, 0),  # π0 on proton at 45°
        (1200, 60, -1, 1),  # π+ on neutron at 60°
    ]

    results = []
    for sqrtS, theta, iso, charge in test_cases:
        result = test_single_case(sqrtS, theta, iso, charge)
        results.append(result)
        print("\n")

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    n_pass = sum(1 for r in results if r)
    n_total = len(results)
    print(f"Tests passed: {n_pass}/{n_total}")

    if n_pass == n_total:
        print("\n✓ All tests PASSED!")
    else:
        print(f"\n✗ {n_total - n_pass} test(s) FAILED")


if __name__ == "__main__":
    sys.exit(main())
