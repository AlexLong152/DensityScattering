# -*- coding: utf-8 -*-
"""Generate LaTeX tables for pion photoproduction results with split uncertainties.

v2: Standardized uncertainty slots.
  unc1 = nuclear structure / form factor uncertainty
  unc2 = elementary amplitude uncertainty
  TDA results always show both parentheses; literature shows only nonzero ones.
"""

import math
import sys
import numpy as np
from FFSlib import get_literature_results, NoResult, nuclei


_nuc_labels = {
    "3H": r"${}^3\mathrm{H}$",
    "3He": r"${}^3\mathrm{He}$",
    "6Li": r"${}^6\mathrm{Li}$",
}

_spin_values = {"3H": 0.5, "3He": 0.5, "6Li": 1.0}


def _determine_precision(unc):
    """Determine number of decimal places from an uncertainty value (2 sig figs, trim trailing zero).

    Only trims trailing zeros when n_dec > 1 to avoid losing the decimal point
    for larger uncertainties (e.g. 1.98 → keep (20) at 1 decimal, don't reduce to (2) at 0).
    """
    mag = math.floor(math.log10(abs(unc)))
    n_dec = -(mag - 1)  # 2 sig figs
    unc_int = int(round(unc * 10**n_dec))
    if unc_int % 10 == 0 and unc_int >= 10 and n_dec > 1:
        n_dec -= 1
    return max(0, n_dec)


def _format_split(val, precision=None, force_two_parens=False):
    """Format [mean, unc1, unc2] as $mean(unc1)(unc2)$.

    When force_two_parens=True: always show both parentheses, even if one rounds to 0.
    When force_two_parens=False: show parentheses only for nonzero uncertainties.
    Returns '--' for NoResult.
    """
    if isinstance(val, NoResult):
        return "--"

    mean, unc1, unc2 = val[0], val[1], val[2]
    has_1 = unc1 > 0
    has_2 = unc2 > 0

    if not has_1 and not has_2:
        p = precision if precision is not None else 3
        if force_two_parens:
            return f"${mean:.{p}f}(0)(0)$"
        return f"${mean:.{p}f}$"

    if precision is not None:
        n_dec = precision
    else:
        main_unc = max(unc1, unc2) if (has_1 and has_2) else (unc1 if has_1 else unc2)
        n_dec = _determine_precision(main_unc)

    scale = 10**n_dec
    mean_r = round(mean, n_dec)
    u1_int = int(round(unc1 * scale))
    u2_int = int(round(unc2 * scale))

    if force_two_parens:
        return f"${mean_r:.{n_dec}f}({u1_int})({u2_int})$"

    # Non-forced: show parentheses only for nonzero uncertainties
    if has_1 and has_2:
        if u1_int == 0 and u2_int == 0:
            return f"${mean_r:.{n_dec}f}$"
        if u1_int == 0:
            return f"${mean_r:.{n_dec}f}({u2_int})$"
        if u2_int == 0:
            return f"${mean_r:.{n_dec}f}({u1_int})$"
        return f"${mean_r:.{n_dec}f}({u1_int})({u2_int})$"
    elif has_1:
        if u1_int == 0:
            return f"${mean_r:.{n_dec}f}$"
        return f"${mean_r:.{n_dec}f}({u1_int})$"
    else:
        if u2_int == 0:
            return f"${mean_r:.{n_dec}f}$"
        return f"${mean_r:.{n_dec}f}({u2_int})$"


def _format_val_unc(mean, unc, precision=None):
    """Format a single value with one combined uncertainty as $mean(unc)$."""
    if unc <= 0:
        p = precision if precision is not None else 3
        return f"${mean:.{p}f}$"

    if precision is not None:
        n_dec = precision
    else:
        n_dec = _determine_precision(unc)

    scale = 10**n_dec
    mean_r = round(mean, n_dec)
    u_int = int(round(unc * scale))
    return f"${mean_r:.{n_dec}f}({u_int})$"


def _format_central(mean, precision=2):
    """Format just a central value with no uncertainty."""
    return f"${mean:.{precision}f}$"


def a0Value(E0_mean, E0_unc, j):
    """Compute a_0 = (4/3)*j*(j+1)*E0^2 and propagated uncertainty.

    da_0/dE_0 = (4/3)*j*(j+1)*2*|E_0|, so da_0 = |da_0/dE_0| * dE_0.
    """
    a0 = (4 / 3) * j * (j + 1) * E0_mean**2
    a0_unc = (4 / 3) * j * (j + 1) * 2 * abs(E0_mean) * E0_unc
    return a0, a0_unc


_ff_quantities_1bod = {"F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V"}
_ff_quantities_2bod = {"F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L"}


def _tda_precision(q, val, base_prec=2):
    """Return precision for a TDA value: use 3 for small form factors, else base_prec."""
    ff_keys = _ff_quantities_1bod | _ff_quantities_2bod
    if q in ff_keys and not isinstance(val, NoResult) and abs(val[0]) < 0.1:
        return 3
    return base_prec


# ============================================================
# One-body tables
# ============================================================


def _onebody_table(d, nuc):
    """Generate LaTeX one-body table for a given nucleus."""
    nuc_label = _nuc_labels[nuc]

    if nuc == "6Li":
        authors = [("TDA Method", "TDA"), ("Braun", "braun")]
    else:
        authors = [("TDA Method", "TDA"), ("Lenkewitz", "lenke"), ("Braun", "braun")]

    quantities = ["F_T^S+V", "F_T^S-V", "E1N", "F_L^S+V", "F_L^S-V", "L1N"]
    col_labels = [
        r"$F_T^{S+V}$",
        r"$F_T^{S-V}$",
        r"$E_{0+}^{1\mathrm{N}}$",
        r"$F_L^{S+V}$",
        r"$F_L^{S-V}$",
        r"$L_{0+}^{1\mathrm{N}}$",
    ]

    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"{\setlength{\extrarowheight}{3pt}")
    lines.append(r"\begin{tabular}{|l||l|l|l||l|l|l|}\hline")

    header_cells = [rf"{nuc_label} \rule[-1.1ex]{{0pt}}{{2.8ex}}"] + col_labels
    lines.append(" & ".join(header_cells) + r" \\")
    lines.append(r"\thickline{0.15em}")

    for author_label, author_key in authors:
        cells = [author_label]
        ftp = author_key == "TDA"
        for q in quantities:
            val = d[author_key][nuc][q]
            prec = _tda_precision(q, val) if author_key == "TDA" else None
            cells.append(_format_split(val, precision=prec, force_two_parens=ftp))
        lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append("}")

    caption = (
        rf"One-body results for {nuc_label}. Form"
        rf" factors $F_{{T/L}}$ are unitless. $E_{{0+}}^{{1\mathrm{{N}}}}$"
        rf" and $L_{{0+}}^{{1\mathrm{{N}}}}$ in units of"
        rf" $[10^{{-3}}/m_{{\pi^+}}]$."
        r" For TDA form factors, the first (second) uncertainty is from"
        r" the potential variation (numerical precision)."
        r" For $E_{0+}^{1\mathrm{N}}$ and $L_{0+}^{1\mathrm{N}}$,"
        r" the first (second) uncertainty is from the nuclear structure"
        r" (elementary single-nucleon amplitudes)."
    )
    label_map = {"3H": "tab:1N_3H", "3He": "tab:1N_3He", "6Li": "tab:1N_6Li"}
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label_map[nuc]}}}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ============================================================
# Two-body tables
# ============================================================


def _twobody_table(d, nuc):
    """Generate LaTeX two-body table for a given nucleus."""
    nuc_label = _nuc_labels[nuc]

    if nuc == "6Li":
        authors_q3 = [("TDA Method", "TDA"), (r"Braun$\times\sqrt{2}$", "braun")]
        authors_q4 = [("TDA Method", "TDA")]
        scale = {"braun": math.sqrt(2)}
    else:
        authors_q3 = [("TDA Method", "TDA"), ("Lenkewitz", "lenke"), ("Braun", "braun")]
        authors_q4 = [("TDA Method", "TDA"), ("Lenkewitz", "lenke")]
        scale = None

    base_quantities = ["F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L", "E2N", "L2N"]
    col_labels = [
        r"$F_T^{(a)}-F_T^{(b)}$",
        r"$F_L^{(a)}-F_L^{(b)}$",
        r"$E_{0+}^{2\mathrm{N}}$",
        r"$L_{0+}^{2\mathrm{N}}$",
    ]
    ncols = len(base_quantities)
    cline_end = 2 + ncols

    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"{\setlength{\extrarowheight}{3pt}")
    lines.append(r"\begin{tabular}{|l|l||l|l|l|l|}\hline")

    merged = rf"\multicolumn{{2}}{{|c||}}{{{nuc_label} \rule[-1.1ex]{{0pt}}{{2.8ex}}}}"
    header = " & ".join([merged] + col_labels) + r" \\"
    lines.append(header)
    lines.append(r"\thickline{0.15em}")

    # O(q^3) group
    active_q3 = [
        (label, key)
        for label, key in authors_q3
        if any(
            not isinstance(d[key][nuc][q + " O(q^3)"], NoResult)
            for q in base_quantities
        )
    ]
    for j, (author_label, author_key) in enumerate(active_q3):
        if j == 0:
            group_cell = rf"\multirow{{{len(active_q3)}}}{{*}}{{$\mathcal{{O}}(q^3)$}}"
        else:
            group_cell = ""
        cells = [group_cell, author_label]
        ftp = author_key == "TDA"
        for q in base_quantities:
            val = d[author_key][nuc][q + " O(q^3)"]
            if scale and author_key in scale and not isinstance(val, NoResult):
                val = val * scale[author_key]
            prec = _tda_precision(q, val) if author_key == "TDA" else None
            cells.append(_format_split(val, precision=prec, force_two_parens=ftp))
        lines.append(" & ".join(cells) + r" \\")
        if j < len(active_q3) - 1:
            lines.append(rf"\cline{{2-{cline_end}}}")

    lines.append(r"\hline\hline")

    # O(q^4) group
    active_q4 = [
        (label, key)
        for label, key in authors_q4
        if any(
            not isinstance(d[key][nuc][q + " O(q^4)"], NoResult)
            for q in base_quantities
        )
    ]
    for j, (author_label, author_key) in enumerate(active_q4):
        if j == 0:
            group_cell = rf"\multirow{{{len(active_q4)}}}{{*}}{{$\mathcal{{O}}(q^4)$}}"
        else:
            group_cell = ""
        cells = [group_cell, author_label]
        ftp = author_key == "TDA"
        for q in base_quantities:
            val = d[author_key][nuc][q + " O(q^4)"]
            if scale and author_key in scale and not isinstance(val, NoResult):
                val = val * scale[author_key]
            prec = _tda_precision(q, val) if author_key == "TDA" else None
            cells.append(_format_split(val, precision=prec, force_two_parens=ftp))
        lines.append(" & ".join(cells) + r" \\")
        if j < len(active_q4) - 1:
            lines.append(rf"\cline{{2-{cline_end}}}")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append("}")

    caption = (
        rf"Two-body results for {nuc_label}."
        rf" $F_{{T/L}}^{{(a)}}-F_{{T/L}}^{{(b)}}$ in units of"
        rf" $[\mathrm{{fm}}^{{-1}}]$. $E_{{0+}}^{{2\mathrm{{N}}}}$ and"
        rf" $L_{{0+}}^{{2\mathrm{{N}}}}$ in units of"
        rf" $[10^{{-3}}/m_{{\pi^+}}]$. $\mathcal{{O}}(q^4)$ results omit"
        r" recoil contributions."
        r" For TDA form factors, the first (second) uncertainty is from"
        r" the potential variation (numerical precision)."
    )
    if nuc == "6Li":
        caption += (
            r" Values from Braun multiplied by $\sqrt{2}$ in accordance with"
            r" discussion in sec.~\ref{ss:BraunSqrt2}."
        )

    label_map = {"3H": "tab:2N_3H", "3He": "tab:2N_3He", "6Li": "tab:2N_6Li"}
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label_map[nuc]}}}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ============================================================
# Summary table
# ============================================================


def _summary_table(d):
    """Generate LaTeX summary table combining 1N and 2N results.

    All columns except a_0 show central values only (no uncertainties).
    a_0 shows one combined uncertainty propagated from E_{0+}.
    Detailed uncertainty breakdowns are in the individual tables.
    """
    twobody_orders = {"TDA": " O(q^4)", "lenke": " O(q^4)", "braun": " O(q^3)"}

    col_labels = [
        r"$E_{0+}^{1\mathrm{N}}$",
        r"$E_{0+}^{2\mathrm{N}}$",
        r"$E_{0+}$",
        r"$L_{0+}^{1\mathrm{N}}$",
        r"$L_{0+}^{2\mathrm{N}}$",
        r"$L_{0+}$",
        r"$a_0\;(\gamma \pi^0)$",
    ]
    ncols = len(col_labels)
    cline_end = 2 + ncols

    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"{\setlength{\extrarowheight}{3pt}")
    lines.append(rf"\begin{{tabular}}{{|l|l||{'l|' * ncols}}}\hline")

    merged = r"\multicolumn{2}{|c||}{\rule[-1.1ex]{0pt}{2.8ex}}"
    header = " & ".join([merged] + col_labels) + r" \\"
    lines.append(header)
    lines.append(r"\thickline{0.15em}")

    for ni, nuc in enumerate(nuclei):
        nuc_label = _nuc_labels[nuc]

        if nuc == "6Li":
            authors = [("TDA", "TDA"), ("Braun*", "braun")]
            scale_2b = {"braun": math.sqrt(2)}
        else:
            authors = [("TDA", "TDA"), ("Lenkewitz", "lenke"), ("Braun", "braun")]
            scale_2b = None

        # Filter to authors that have at least some data
        active = []
        for label, key in authors:
            order = twobody_orders.get(key, " O(q^3)")
            has_data = any(
                not isinstance(d[key][nuc][q], NoResult)
                for q in ["E1N", f"E2N{order}", "L1N", f"L2N{order}"]
            )
            if has_data:
                active.append((label, key))

        for j, (author_label, author_key) in enumerate(active):
            if j == 0:
                group_cell = rf"\multirow{{{len(active)}}}{{*}}{{{nuc_label}}}"
            else:
                group_cell = ""

            cells = [group_cell, author_label]
            order = twobody_orders.get(author_key, " O(q^3)")
            is_tda = author_key == "TDA"

            # E1N (central value only)
            e1n = d[author_key][nuc]["E1N"]
            cells.append("--" if isinstance(e1n, NoResult) else _format_central(e1n[0]))

            # E2N (central value only, optional scaling)
            e2n = d[author_key][nuc][f"E2N{order}"]
            if scale_2b and author_key in scale_2b and not isinstance(e2n, NoResult):
                e2n = e2n * scale_2b[author_key]
            cells.append("--" if isinstance(e2n, NoResult) else _format_central(e2n[0]))

            # E_0+ = E1N + E2N (central value only)
            if not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult):
                e_mean = e1n[0] + e2n[0]
                cells.append(_format_central(e_mean))
            else:
                cells.append("--")

            # L1N (central value only)
            l1n = d[author_key][nuc]["L1N"]
            cells.append("--" if isinstance(l1n, NoResult) else _format_central(l1n[0]))

            # L2N (central value only)
            l2n = d[author_key][nuc][f"L2N{order}"]
            cells.append("--" if isinstance(l2n, NoResult) else _format_central(l2n[0]))

            # L_0+ = L1N + L2N (central value only)
            if not isinstance(l1n, NoResult) and not isinstance(l2n, NoResult):
                l_mean = l1n[0] + l2n[0]
                cells.append(_format_central(l_mean))
            else:
                cells.append("--")

            # a_0 (from E_0+ only, single combined uncertainty)
            if not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult):
                e_mean = e1n[0] + e2n[0]
                e_total = math.sqrt(
                    e1n[1] ** 2 + e1n[2] ** 2 + e2n[1] ** 2 + e2n[2] ** 2
                )
                a0, a0_unc = a0Value(e_mean, e_total, _spin_values[nuc])
                cells.append(_format_val_unc(a0, a0_unc, precision=2))
            else:
                cells.append("--")

            lines.append(" & ".join(cells) + r" \\")
            if j < len(active) - 1:
                lines.append(rf"\cline{{2-{cline_end}}}")

        if ni < len(nuclei) - 1:
            lines.append(r"\hline\hline")
        else:
            lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append("}")

    caption = (
        r"Summary of photoproduction results. $E_{0+}$ and"
        r" $L_{0+}$ in units of $[10^{-3}/m_{\pi^+}]$. $a_0$"
        r" in units of $[10^{-6}/m_{\pi^+}^2]$."
        r" $E_{0+} = E_{0+}^{1\mathrm{N}} + E_{0+}^{2\mathrm{N}}$,"
        r" $L_{0+} = L_{0+}^{1\mathrm{N}} + L_{0+}^{2\mathrm{N}}$."
        r" All $1\mathrm{N}$ values at $\mathcal{O}(q^4)$."
        r" TDA and Lenkewitz $2\mathrm{N}$ results at"
        r" $\mathcal{O}(q^4)$; Braun at $\mathcal{O}(q^3)$."
        r" ${}^*$The ${}^6\mathrm{Li}$ results of"
        r" Braun have been multiplied by $\sqrt{2}$ to correct"
        r" for the spin-1 matrix normalization error identified"
        r" in sec.~\ref{ss:BraunSqrt2}. See"
        r" Tables~\ref{tab:1N_3H}--\ref{tab:2N_6Li} for"
        r" uncertainty breakdowns."
    )
    lines.append(rf"\caption{{{caption}}}")
    lines.append(r"\label{piphoto_summary}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ============================================================
# Test mode
# ============================================================


def _test_formatting():
    """Test formatting with mock data to verify _format_split behavior."""
    passed = 0
    failed = 0

    def check(label, result, expected):
        nonlocal passed, failed
        ok = result == expected
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {label}: {result}")
        if not ok:
            print(f"         expected: {expected}")
            failed += 1
        else:
            passed += 1

    print("Testing _format_split formatting:\n")

    # Form factor: theory=0.06, numeric=0.02 → (6)(2) at precision=2
    val = np.array([1.23, 0.06, 0.02])
    check(
        "FF theory=0.06, numeric=0.02",
        _format_split(val, precision=2, force_two_parens=True),
        "$1.23(6)(2)$",
    )

    # Form factor: theory=0.001, numeric=0.02 → (0)(2) at precision=2
    val = np.array([1.23, 0.001, 0.02])
    check(
        "FF theory=0.001, numeric=0.02",
        _format_split(val, precision=2, force_two_parens=True),
        "$1.23(0)(2)$",
    )

    # E1N: unc1=0.10, unc2=0.09 → (10)(9)
    val = np.array([1.23, 0.10, 0.09])
    check(
        "E1N unc1=0.10, unc2=0.09",
        _format_split(val, precision=2, force_two_parens=True),
        "$1.23(10)(9)$",
    )

    # Literature: unc1=0.05, unc2=0.0 → (5) single paren
    val = np.array([1.23, 0.05, 0.0])
    check(
        "Lit unc1=0.05, unc2=0.0",
        _format_split(val, precision=None, force_two_parens=False),
        "$1.23(5)$",
    )

    # E2N TDA: unc1=0.43, unc2=0 → (43)(0) with force_two_parens
    val = np.array([-3.55, 0.43, 0.0])
    check(
        "E2N TDA unc1=0.43, unc2=0",
        _format_split(val, precision=2, force_two_parens=True),
        "$-3.55(43)(0)$",
    )

    # Literature with both nonzero
    val = np.array([1.71, 0.04, 0.09])
    check(
        "Lit unc1=0.04, unc2=0.09",
        _format_split(val, precision=None, force_two_parens=False),
        "$1.71(4)(9)$",
    )

    # NoResult
    check(
        "NoResult", _format_split(NoResult(), precision=2, force_two_parens=True), "--"
    )

    # Both zero, not forced
    val = np.array([1.23, 0.0, 0.0])
    check(
        "Both zero, not forced",
        _format_split(val, precision=2, force_two_parens=False),
        "$1.23$",
    )

    # Both zero, forced
    check(
        "Both zero, forced",
        _format_split(val, precision=2, force_two_parens=True),
        "$1.23(0)(0)$",
    )

    print(f"\n{passed} passed, {failed} failed")
    return failed == 0


def main():
    d = get_literature_results()
    for nuc in nuclei:
        print(_onebody_table(d, nuc))
        print()
        print(_twobody_table(d, nuc))
        print()
    print(_summary_table(d))


if __name__ == "__main__":
    if "--test" in sys.argv:
        success = _test_formatting()
        sys.exit(0 if success else 1)
    else:
        main()
