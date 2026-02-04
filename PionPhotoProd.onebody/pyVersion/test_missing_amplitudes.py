# -*- coding: utf-8 -*-
"""
Test that Python code returns zero/negligible amplitudes for combinations
where the Fortran code reports missing amplitude data.

This verifies the behavior seen in Fortran warnings like:
    WARNING: No amplitude data found for:
      target = 32q
      ell = 3
      sqrtS = 1073.8889700000000

The Fortran warning occurs when sqrtS is below the minimum energy in the
SAID data tables. The Python code extrapolates, giving very small values.

@author: alexl
"""

import sys
import numpy as np
import PionPhotoLib as ppl


def get_sqrtS_range(poleData, sign, target, ell, pole_type="E"):
    """Get the sqrtS range for a given amplitude data array."""
    if pole_type == "E":
        arr = poleData.Epole[sign][target].get(ell)
    else:
        arr = poleData.Mpole[sign][target].get(ell)

    if arr is None or arr.shape[1] == 0:
        return None, None

    sqrtS_values = np.real(arr[0])
    return sqrtS_values[0], sqrtS_values[-1]


def test_missing_amplitude_32q_ell3():
    """
    Test that target=32q, ell=3 returns effectively zero amplitudes
    at sqrtS=1073.88897 MeV.

    This matches the Fortran warning:
        WARNING: No amplitude data found for:
          target = 32q
          ell = 3
          sqrtS = 1073.8889700000000

    The Fortran code returns exactly zero when sqrtS is outside the data range.
    The Python code extrapolates using the first two points, giving very small
    values that are effectively zero for physical purposes.
    """
    # Load SAID pole data
    poleData = ppl.SaidPoles()

    # Test parameters from Fortran warning
    target = "32q"
    ell = 3
    sqrtS = 1073.8889700000000

    # Check the data range first
    print(f"Testing amplitude data for:")
    print(f"  target = {target}")
    print(f"  ell    = {ell}")
    print(f"  sqrtS  = {sqrtS}")
    print()

    # Show data ranges
    print("Data ranges (sqrtS in MeV):")
    for sign in ["plus", "minus"]:
        for pole_type in ["E", "M"]:
            min_s, max_s = get_sqrtS_range(poleData, sign, target, ell, pole_type)
            if min_s is not None:
                in_range = "IN RANGE" if min_s <= sqrtS <= max_s else "OUT OF RANGE"
                print(f"  {pole_type}{sign}: {min_s:.2f} - {max_s:.2f}  [{in_range}]")
            else:
                print(f"  {pole_type}{sign}: no data")
    print()

    # Get poles using the Python function
    Eplus, Mplus, Eminus, Mminus = ppl.getPoles(poleData, target, ell, sqrtS)

    print(f"Results from Python getPoles:")
    print(f"  Eplus  = {Eplus}")
    print(f"  Mplus  = {Mplus}")
    print(f"  Eminus = {Eminus}")
    print(f"  Mminus = {Mminus}")
    print()

    # Check that all amplitudes are effectively zero
    # Use a tolerance that considers these as "negligible" for physics purposes
    # Values < 1e-6 are essentially zero in the context of these amplitudes
    all_negligible = True
    tol = 1e-6  # More realistic tolerance for "effectively zero"

    if abs(Eplus) > tol:
        print(f"UNEXPECTED: Eplus = {Eplus} (expected negligible)")
        all_negligible = False
    if abs(Mplus) > tol:
        print(f"UNEXPECTED: Mplus = {Mplus} (expected negligible)")
        all_negligible = False
    if abs(Eminus) > tol:
        print(f"UNEXPECTED: Eminus = {Eminus} (expected negligible)")
        all_negligible = False
    if abs(Mminus) > tol:
        print(f"UNEXPECTED: Mminus = {Mminus} (expected negligible)")
        all_negligible = False

    if all_negligible:
        print("PASS: All amplitudes are negligible (< 1e-6)")
        print("      This matches the Fortran behavior of returning zero")
        print("      when sqrtS is below the data range.")

    return all_negligible


def test_check_available_data_32q():
    """
    Check what ell values have data for target=32q.
    This helps understand why ell=3 is missing.
    """
    poleData = ppl.SaidPoles()
    target = "32q"

    print("\nAvailable data for target='32q':")
    print("-" * 50)

    for sign in ["plus", "minus"]:
        print(f"\n{sign.upper()}:")
        for ell in range(6):
            Epole_arr = poleData.Epole[sign][target].get(ell, None)
            Mpole_arr = poleData.Mpole[sign][target].get(ell, None)

            E_has_data = (
                Epole_arr is not None
                and isinstance(Epole_arr, np.ndarray)
                and Epole_arr.shape[1] > 0
            )
            M_has_data = (
                Mpole_arr is not None
                and isinstance(Mpole_arr, np.ndarray)
                and Mpole_arr.shape[1] > 0
            )

            E_status = f"E: {Epole_arr.shape[1]} points" if E_has_data else "E: no data"
            M_status = f"M: {Mpole_arr.shape[1]} points" if M_has_data else "M: no data"

            print(f"  ell={ell}: {E_status}, {M_status}")


def test_multiple_ell_values():
    """
    Show amplitude values for all ell values at sqrtS below the data range.

    This is informational - shows that higher ell values extrapolate
    to smaller values when sqrtS is just below the data range.
    """
    poleData = ppl.SaidPoles()
    target = "32q"
    sqrtS = 1073.8889700000000

    print("\nAmplitude values for target='32q' at sqrtS below data range:")
    print("-" * 70)
    print(f"  Test sqrtS = {sqrtS:.2f} MeV")
    print()

    for ell in range(6):
        Eplus, Mplus, Eminus, Mminus = ppl.getPoles(poleData, target, ell, sqrtS)

        # Get data range
        min_s, max_s = get_sqrtS_range(poleData, "plus", target, ell, "E")
        if min_s is not None:
            range_str = f"data: {min_s:.1f}-{max_s:.1f}"
        else:
            range_str = "no data"

        print(f"  ell={ell} ({range_str}):")
        print(
            f"    Eplus={Eplus:.2e}, Mplus={Mplus:.2e}, "
            f"Eminus={Eminus:.2e}, Mminus={Mminus:.2e}"
        )

    # Always return True - this is informational
    return True


def test_viewData():
    """
    Demonstrate the viewData function for the Fortran warning case.
    """
    poleData = ppl.SaidPoles()
    print("\nUsing viewData to show detailed lookup information:")
    print()
    # s1 = 1079.9894492877338
    s1 = 1073.8889700000000
    ell = 0
    for ell in range(2):
        ppl.viewData(poleData, "p12", ell, s1)
        # ppl.viewData(poleData, "p12", ell, s2)


def main():
    print("=" * 60)
    print("Testing Missing Amplitude Data (Python vs Fortran behavior)")
    print("=" * 60)

    test_viewData()

    # Run the main test for the specific Fortran warning case
    # test1_passed = test_missing_amplitude_32q_ell3()

    # Show detailed lookup info using viewData

    # Show available data structure
    # test_check_available_data_32q()
    #
    # # Show all ell values (informational)
    # test_multiple_ell_values()
    #
    # # Summary
    # print("\n" + "=" * 60)
    # print("SUMMARY")
    # print("=" * 60)
    #
    # if test1_passed:
    #     print("Test PASSED - Python returns negligible amplitudes for")
    #     print("  target=32q, ell=3, sqrtS=1073.88897")
    #     print()
    #     print("This matches the Fortran warning about missing amplitude data.")
    #     print("The Fortran code returns exactly 0, Python extrapolates to ~1e-9.")
    #     print("Both are effectively zero for physics purposes.")
    # else:
    #     print("Test FAILED")
    #     return False
    #
    # return True


if __name__ == "__main__":
    success = main()
