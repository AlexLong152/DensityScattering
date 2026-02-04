# -*- coding: utf-8 -*-
"""
Compare getRawM: Python vs Fortran implementation.

This script tests whether the Python implementation of getRawM matches
the Fortran implementation by calling the testRawM executable.

@author: alexl
"""

import numpy as np
import os
import subprocess
import PionPhotoLib as ppl

# Constants (matching Fortran constants.def)
mpi0 = 134.97
mpi = 139.5675
Mproton = 938.27231
Mneutron = 939.56563


def ensure_fortran_executable():
    """
    Compile the Fortran test executable using make.
    """
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    print("Compiling testRawM executable...")
    result = os.system(f"cd {parent_dir} && make testRawM")

    if result != 0:
        raise RuntimeError(f"Failed to compile testRawM. Exit code: {result}")

    print("Compilation successful!")
    return os.path.join(parent_dir, "testRawM")


def call_fortran_getRawM(sqrtS, x, nucs, epsVec):
    """
    Call the Fortran testRawM executable.

    Parameters
    ----------
    sqrtS : float
        Center-of-mass energy in MeV
    x : float
        cos(theta) scattering angle
    nucs : str
        Reaction type ('pp0', 'nn0', 'pn+', 'np-')
    epsVec : array-like
        Polarization vector (3 real components)

    Returns
    -------
    numpy.ndarray
        2x2 complex matrix
    """
    fortran_exe = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "testRawM"
    )

    # Prepare input: sqrtS x nucs eps1 eps2 eps3
    input_str = (
        f"{sqrtS} {x} {nucs} {epsVec[0].real} {epsVec[1].real} {epsVec[2].real}\n"
    )

    try:
        result = subprocess.run(
            [fortran_exe], input=input_str, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error calling Fortran program: {e}")
        print(f"Stderr: {e.stderr}")
        raise

    # Parse output
    mat_fortran = {}
    for line in result.stdout.strip().split("\n"):
        if "Mout(" in line:
            parts = line.split("=")
            key = parts[0].strip()
            values = parts[1].strip().split()
            real_part = float(values[0])
            imag_part = float(values[1])

            # Extract indices from "Mout(-1,-1)" or "Mout( 1, 1)"
            indices = key.replace("Mout(", "").replace(")", "").split(",")
            i = int(indices[0].strip())
            j = int(indices[1].strip())
            mat_fortran[(i, j)] = complex(real_part, imag_part)

    # Convert to Python 2x2 matrix
    mat = np.zeros((2, 2), dtype=complex)
    mat[0, 0] = mat_fortran[(-1, -1)]
    mat[0, 1] = mat_fortran[(-1, 1)]
    mat[1, 0] = mat_fortran[(1, -1)]
    mat[1, 1] = mat_fortran[(1, 1)]

    return mat


def get_threshold_sqrtS(nucs):
    """Calculate threshold sqrtS for a given reaction."""
    if nucs == "pp0":
        return mpi0 + Mproton
    elif nucs == "nn0":
        return mpi0 + Mneutron
    elif nucs == "pn+":
        return mpi + Mproton
    elif nucs == "np-":
        return mpi + Mneutron
    else:
        raise ValueError(f"Unknown reaction: {nucs}")


def compare_matrices(py_mat, fort_mat, sqrtS, x, nucs, tol=2e-3):
    """Compare Python and Fortran matrices. Only prints output if test fails."""
    # Calculate differences
    diff = np.abs(py_mat - fort_mat)
    max_diff = np.max(diff)

    # Relative differences for non-zero elements
    rel_diffs = []
    for i in range(2):
        for j in range(2):
            if np.abs(fort_mat[i, j]) > 1e-15:
                rel_diff = np.abs(py_mat[i, j] - fort_mat[i, j]) / np.abs(
                    fort_mat[i, j]
                )
                rel_diffs.append(rel_diff)

    max_rel_diff = max(rel_diffs) if rel_diffs else 0.0
    passed = max_rel_diff < tol

    # Only print details for failing tests
    if not passed:
        print(f"\nsqrtS={sqrtS:.2f} MeV, x={x:.2f}, nucs='{nucs}'")
        print("-" * 50)

        print("Python:")
        for i in range(2):
            spin_i = -1 if i == 0 else 1
            for j in range(2):
                spin_j = -1 if j == 0 else 1
                print(f"  [{spin_i:2d},{spin_j:2d}] = {py_mat[i, j]:.10e}")

        print("Fortran:")
        for i in range(2):
            spin_i = -1 if i == 0 else 1
            for j in range(2):
                spin_j = -1 if j == 0 else 1
                print(f"  [{spin_i:2d},{spin_j:2d}] = {fort_mat[i, j]:.10e}")

        print(f"Max abs diff: {max_diff:.6e}, Max rel diff: {max_rel_diff:.6e}")
        print("FAILED")

    return passed


def main():
    """Run comparison tests."""
    print("=" * 60)
    print("Testing getRawM: Python vs Fortran")
    print("=" * 60)

    # Compile Fortran executable
    ensure_fortran_executable()

    # Load pole data
    print("\nLoading SAID pole data...")
    poleData = ppl.SaidPoles()

    # Polarization vector
    eps_x = np.array([1, 0, 0], dtype=complex)

    # Test cases: (nucs, sqrtS_values, x_values)
    reactions = ["pp0", "nn0"]
    x_values = np.linspace(-1, 1, num=11, endpoint=True)

    results = []

    for nucs in reactions:
        sqrtS_thresh = get_threshold_sqrtS(nucs)
        sqrtS_values = [
            sqrtS_thresh + 0.1,
            sqrtS_thresh + 2,
            sqrtS_thresh + 5,
            1100.0,
            1150.0,
            1200.0,
            1250.0,
        ]

        for sqrtS in sqrtS_values:
            for x in x_values:
                # Python result (roundThresh=True to match Fortran behavior)
                py_mat = ppl.getRawM(
                    sqrtS, x, nucs, poleData, epsVec=eps_x, roundThresh=True
                )

                # Fortran result
                fort_mat = call_fortran_getRawM(sqrtS, x, nucs, eps_x)

                # Compare
                passed = compare_matrices(py_mat, fort_mat, sqrtS, x, nucs)
                results.append(passed)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    n_passed = sum(results)
    n_total = len(results)
    print(f"Tests passed: {n_passed}/{n_total}")

    if n_passed == n_total:
        print("\nALL TESTS PASSED!")
    else:
        print(f"\n{n_total - n_passed} test(s) FAILED")


if __name__ == "__main__":
    main()
