# -*- coding: utf-8 -*-
"""Print pion photoproduction results to terminal in readable format."""

import math
import numpy as np
from FFSlib import get_literature_results, NoResult, nuclei
from latex_tables import a0Value


_spin_values = {"3H": 0.5, "3He": 0.5, "6Li": 1.0}


def _fmt(val, precision=None, force_two_parens=False):
    """Format [mean, unc1, unc2] as parenthetical string for terminal output.

    force_two_parens=True: always show both parens (for TDA).
    force_two_parens=False: show parens only for nonzero uncertainties (for literature).
    """
    if isinstance(val, NoResult):
        return "---"

    mean, u1, u2 = val[0], val[1], val[2]
    has_1 = u1 > 0
    has_2 = u2 > 0

    if not has_1 and not has_2:
        p = precision if precision is not None else 3
        if force_two_parens:
            return f"{mean:.{p}f}(0)(0)"
        return f"{mean:.{p}f}"

    if precision is not None:
        n_dec = precision
    else:
        main_unc = max(u1, u2) if (has_1 and has_2) else (u1 if has_1 else u2)
        mag = math.floor(math.log10(abs(main_unc)))
        n_dec = -(mag - 1)
        unc_int = int(round(main_unc * 10**n_dec))
        if unc_int % 10 == 0 and unc_int >= 10 and n_dec > 1:
            n_dec -= 1
        n_dec = max(0, n_dec)

    scale = 10**n_dec
    mean_r = round(mean, n_dec)
    u1_int = int(round(u1 * scale))
    u2_int = int(round(u2 * scale))

    if force_two_parens:
        return f"{mean_r:.{n_dec}f}({u1_int})({u2_int})"

    if has_1 and has_2:
        if u1_int == 0 and u2_int == 0:
            return f"{mean_r:.{n_dec}f}"
        if u1_int == 0:
            return f"{mean_r:.{n_dec}f}({u2_int})"
        if u2_int == 0:
            return f"{mean_r:.{n_dec}f}({u1_int})"
        return f"{mean_r:.{n_dec}f}({u1_int})({u2_int})"
    elif has_1:
        if u1_int == 0:
            return f"{mean_r:.{n_dec}f}"
        return f"{mean_r:.{n_dec}f}({u1_int})"
    else:
        if u2_int == 0:
            return f"{mean_r:.{n_dec}f}"
        return f"{mean_r:.{n_dec}f}({u2_int})"


def _fmt_val_unc(mean, unc, precision=None):
    """Format mean(unc) with a single combined uncertainty."""
    if unc <= 0:
        p = precision if precision is not None else 3
        return f"{mean:.{p}f}"
    if precision is not None:
        n_dec = precision
    else:
        mag = math.floor(math.log10(abs(unc)))
        n_dec = -(mag - 1)
        unc_int = int(round(unc * 10**n_dec))
        if unc_int % 10 == 0 and unc_int >= 10 and n_dec > 1:
            n_dec -= 1
        n_dec = max(0, n_dec)
    scale = 10**n_dec
    mean_r = round(mean, n_dec)
    u_int = int(round(unc * scale))
    return f"{mean_r:.{n_dec}f}({u_int})"


def _print_table(title, col_headers, rows, col_w=22, label_w=None):
    """Print a generic table with a label column and data columns."""
    if label_w is None:
        label_w = max(len(r[0]) for r in rows) + 2
    total_w = label_w + len(col_headers) * col_w

    print(f"\n  {title}")
    print("  " + "=" * total_w)
    header = " " * label_w + "".join(h.rjust(col_w) for h in col_headers)
    print("  " + header)
    print("  " + "-" * total_w)
    for row in rows:
        line = row[0].ljust(label_w) + "".join(c.rjust(col_w) for c in row[1:])
        print("  " + line)
    print("  " + "-" * total_w)


def _onebody(d, nuc):
    """Print one-body table for a nucleus."""
    if nuc == "6Li":
        authors = [("TDA", "TDA"), ("Braun", "braun")]
    else:
        authors = [("TDA", "TDA"), ("Lenke.", "lenke"), ("Braun", "braun")]

    quantities = ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V", "E1N", "L1N"]
    col_headers = [q.replace("^S+V", " S+V").replace("^S-V", " S-V") for q in quantities]

    rows = []
    for label, key in authors:
        prec = 4 if key == "TDA" else None
        ftp = key == "TDA"
        cells = [label]
        for q in quantities:
            cells.append(_fmt(d[key][nuc][q], precision=prec, force_two_parens=ftp))
        rows.append(cells)

    _print_table(f"{nuc} One-body", col_headers, rows)


def _twobody(d, nuc):
    """Print two-body table for a nucleus."""
    if nuc == "6Li":
        authors_q3 = [("TDA", "TDA"), ("Braun*sqrt2", "braun")]
        authors_q4 = [("TDA", "TDA")]
        scale = {"braun": math.sqrt(2)}
    else:
        authors_q3 = [("TDA", "TDA"), ("Lenke.", "lenke"), ("Braun", "braun")]
        authors_q4 = [("TDA", "TDA"), ("Lenke.", "lenke")]
        scale = None

    base_q = ["F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L", "E2N", "L2N"]
    col_headers = ["dF_T", "dF_L", "E2N", "L2N"]

    rows = []
    for order_tag, author_list in [(" O(q^3)", authors_q3), (" O(q^4)", authors_q4)]:
        order_label = order_tag.strip()
        for label, key in author_list:
            has_any = any(
                not isinstance(d[key][nuc][q + order_tag], NoResult) for q in base_q
            )
            if not has_any:
                continue
            prec = 4 if key == "TDA" else None
            ftp = key == "TDA"
            cells = [f"{order_label} {label}"]
            for q in base_q:
                val = d[key][nuc][q + order_tag]
                if scale and key in scale and not isinstance(val, NoResult):
                    val = val * scale[key]
                cells.append(_fmt(val, precision=prec, force_two_parens=ftp))
            rows.append(cells)

    _print_table(f"{nuc} Two-body", col_headers, rows, label_w=22)


def _summary(d):
    """Print summary table combining 1N + 2N results."""
    twobody_orders = {"TDA": " O(q^4)", "lenke": " O(q^4)", "braun": " O(q^3)"}
    col_headers = ["E1N", "E2N", "E0+", "L1N", "L2N", "L0+", "a0"]

    rows = []
    for nuc in nuclei:
        if nuc == "6Li":
            authors = [("TDA", "TDA"), ("Braun*", "braun")]
            scale_2b = {"braun": math.sqrt(2)}
        else:
            authors = [("TDA", "TDA"), ("Lenke.", "lenke"), ("Braun", "braun")]
            scale_2b = None

        for label, key in authors:
            order = twobody_orders.get(key, " O(q^3)")
            has_data = any(
                not isinstance(d[key][nuc][q], NoResult)
                for q in ["E1N", f"E2N{order}", "L1N", f"L2N{order}"]
            )
            if not has_data:
                continue

            is_tda = key == "TDA"
            prec = 4 if is_tda else None
            ftp = is_tda

            e1n = d[key][nuc]["E1N"]
            e2n = d[key][nuc][f"E2N{order}"]
            if scale_2b and key in scale_2b and not isinstance(e2n, NoResult):
                e2n = e2n * scale_2b[key]
            l1n = d[key][nuc]["L1N"]
            l2n = d[key][nuc][f"L2N{order}"]

            cells = [f"{nuc} {label}"]

            # E1N, E2N
            cells.append("---" if isinstance(e1n, NoResult) else _fmt(e1n, precision=prec, force_two_parens=ftp))
            cells.append("---" if isinstance(e2n, NoResult) else _fmt(e2n, precision=prec, force_two_parens=ftp))

            # E0+
            if not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult):
                e_mean = e1n[0] + e2n[0]
                if is_tda:
                    e_unc1 = math.sqrt(e1n[1] ** 2 + e2n[1] ** 2)
                    e_unc2 = e1n[2]
                    cells.append(_fmt(np.array([e_mean, e_unc1, e_unc2]), precision=4, force_two_parens=True))
                else:
                    e_unc = math.sqrt(e1n[1]**2 + e1n[2]**2 + e2n[1]**2 + e2n[2]**2)
                    cells.append(_fmt_val_unc(e_mean, e_unc, precision=prec))
            else:
                cells.append("---")

            # L1N, L2N
            cells.append("---" if isinstance(l1n, NoResult) else _fmt(l1n, precision=prec, force_two_parens=ftp))
            cells.append("---" if isinstance(l2n, NoResult) else _fmt(l2n, precision=prec, force_two_parens=ftp))

            # L0+
            if not isinstance(l1n, NoResult) and not isinstance(l2n, NoResult):
                l_mean = l1n[0] + l2n[0]
                if is_tda:
                    l_unc1 = math.sqrt(l1n[1] ** 2 + l2n[1] ** 2)
                    l_unc2 = l1n[2]
                    cells.append(_fmt(np.array([l_mean, l_unc1, l_unc2]), precision=4, force_two_parens=True))
                else:
                    l_unc = math.sqrt(l1n[1]**2 + l1n[2]**2 + l2n[1]**2 + l2n[2]**2)
                    cells.append(_fmt_val_unc(l_mean, l_unc, precision=prec))
            else:
                cells.append("---")

            # a0
            if not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult):
                e_mean = e1n[0] + e2n[0]
                if is_tda:
                    e_unc1 = math.sqrt(e1n[1] ** 2 + e2n[1] ** 2)
                    e_unc2 = e1n[2]
                    e_total = math.sqrt(e_unc1**2 + e_unc2**2)
                else:
                    e_total = math.sqrt(e1n[1]**2 + e1n[2]**2 + e2n[1]**2 + e2n[2]**2)
                a0, a0_unc = a0Value(e_mean, e_total, _spin_values[nuc])
                cells.append(_fmt_val_unc(a0, a0_unc, precision=prec))
            else:
                cells.append("---")

            rows.append(cells)

    _print_table("Summary (1N + 2N)", col_headers, rows, label_w=16)


def printAll(d=None):
    if d is None:
        d = get_literature_results()

    for nuc in nuclei:
        _onebody(d, nuc)
        _twobody(d, nuc)

    _summary(d)
    print()


if __name__ == "__main__":
    printAll()
