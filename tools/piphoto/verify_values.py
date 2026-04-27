#!/usr/bin/env python3
"""Print every stored value from FFSlib.get_literature_results() for verification."""

import math
import numpy as np
from FFSlib import (
    get_literature_results,
    NoResult,
    nuclei,
    quantities,
    _analyze_above_files,
    _calc_multipole_unc,
    _analyze_onebody_files,
    _get_files_from_folder,
    above_folders,
    onebody_folders,
    K1N,
    E0p_pi0p,
    E0p_pi0n,
    L0p_pi0p,
    L0p_pi0n,
    amp_frac,
    err,
)

_spin_values = {"3H": 0.5, "3He": 0.5, "6Li": 1.0}
_twobody_orders = {"TDA": " O(q^4)", "lenke": " O(q^4)", "braun": " O(q^3)"}


def fmt(arr):
    return f"mean={arr[0]:.4f}, unc1={arr[1]:.4f}, unc2={arr[2]:.4f}"


def main():
    d = get_literature_results()

    # ── Section 1: All stored values ──
    print("=" * 70)
    print("ALL STORED VALUES")
    print("=" * 70)
    for author in ["TDA", "lenke", "braun"]:
        for nuc in nuclei:
            for q in quantities:
                val = d[author][nuc][q]
                if isinstance(val, NoResult):
                    print(f"  {author} / {nuc} / {q}: [NoResult]")
                else:
                    print(f"  {author} / {nuc} / {q}: {fmt(val)}")

    # ── Section 2: TDA E1N/L1N breakdown ──
    print()
    print("=" * 70)
    print("TDA E1N/L1N UNCERTAINTY BREAKDOWN")
    print("=" * 70)
    for nuc in nuclei:
        print(f"\n  --- {nuc} ---")

        # Raw above-threshold results (theory_spread and numeric from files)
        above_files = _get_files_from_folder(above_folders[nuc])
        above_data = _analyze_above_files(above_files)

        # Onebody data needed for multipole uncertainty
        onebody_files = _get_files_from_folder(onebody_folders[nuc])
        onebody_data = _analyze_onebody_files(onebody_files)

        # Multipole uncertainty from elementary amplitude variation
        multipole_unc = _calc_multipole_unc(nuc, onebody_data)

        k = K1N(nuc)
        print(f"  K1N({nuc}) = {k:.6f}")

        for amp in ["E1N", "L1N"]:
            raw_mean = above_data[amp][0]
            raw_theory_spread = above_data[amp][1]
            raw_numeric = above_data[amp][2]
            m_unc = multipole_unc[amp]
            combined_theory = math.sqrt(raw_theory_spread**2 + m_unc**2)

            print(f"  {amp}:")
            print(f"    mean (from above files)     = {raw_mean:.6f}")
            print(f"    above_theory_spread          = {raw_theory_spread:.6f}")
            print(f"    above_numeric (1.5%)         = {raw_numeric:.6f}")
            print(f"    multipole_unc (5% amp var)   = {m_unc:.6f}")
            print(f"    combined_theory (quadrature) = {combined_theory:.6f}")
            print(
                f"    => stored: [{raw_mean:.4f}, {combined_theory:.4f}, {raw_numeric:.4f}]"
            )
            stored = d["TDA"][nuc][amp]
            print(f"    => verify: [{stored[0]:.4f}, {stored[1]:.4f}, {stored[2]:.4f}]")

    # ── Section 3: Summary / derived quantities ──
    print()
    print("=" * 70)
    print("SUMMARY TABLE DERIVED QUANTITIES")
    print("=" * 70)
    for nuc in nuclei:
        print(f"\n  === {nuc} (j={_spin_values[nuc]}) ===")
        j = _spin_values[nuc]

        for author in ["TDA", "lenke", "braun"]:
            order = _twobody_orders.get(author, " O(q^3)")
            e1n = d[author][nuc]["E1N"]
            e2n_key = f"E2N{order}"
            e2n = d[author][nuc][e2n_key]
            l1n = d[author][nuc]["L1N"]
            l2n_key = f"L2N{order}"
            l2n = d[author][nuc][l2n_key]

            # For 6Li braun, scale 2-body by sqrt(2)
            scale = math.sqrt(2) if (nuc == "6Li" and author == "braun") else 1.0

            has_e = not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult)
            has_l = not isinstance(l1n, NoResult) and not isinstance(l2n, NoResult)

            if not has_e and not has_l:
                continue

            print(f"\n  {author} (2N order: {order.strip()}, scale={scale:.4f}):")

            if has_e:
                e2n_scaled = e2n * scale
                e_mean = e1n[0] + e2n_scaled[0]
                e_unc = math.sqrt(
                    e1n[1] ** 2 + e1n[2] ** 2 + e2n_scaled[1] ** 2 + e2n_scaled[2] ** 2
                )
                print(
                    f"    E1N          = {e1n[0]:.4f}  (unc1={e1n[1]:.4f}, unc2={e1n[2]:.4f})"
                )
                print(
                    f"    E2N{order} = {e2n_scaled[0]:.4f}  (unc1={e2n_scaled[1]:.4f}, unc2={e2n_scaled[2]:.4f})"
                )
                print(f"    E0+ = E1N+E2N = {e_mean:.4f}  (combined_unc={e_unc:.4f})")

                a0 = (4 / 3) * j * (j + 1) * e_mean**2
                a0_unc = (4 / 3) * j * (j + 1) * 2 * abs(e_mean) * e_unc
                print(f"    a0 = (4/3)*j*(j+1)*E0+^2 = {a0:.4f}  (unc={a0_unc:.4f})")

            if has_l:
                l2n_scaled = l2n * scale if not isinstance(l2n, NoResult) else l2n
                l_mean = l1n[0] + l2n_scaled[0]
                l_unc = math.sqrt(
                    l1n[1] ** 2 + l1n[2] ** 2 + l2n_scaled[1] ** 2 + l2n_scaled[2] ** 2
                )
                print(
                    f"    L1N          = {l1n[0]:.4f}  (unc1={l1n[1]:.4f}, unc2={l1n[2]:.4f})"
                )
                print(
                    f"    L2N{order} = {l2n_scaled[0]:.4f}  (unc1={l2n_scaled[1]:.4f}, unc2={l2n_scaled[2]:.4f})"
                )
                print(f"    L0+ = L1N+L2N = {l_mean:.4f}  (combined_unc={l_unc:.4f})")


if __name__ == "__main__":
    main()
