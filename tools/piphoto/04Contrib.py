# -*- coding: utf-8 -*-
"""Print O(delta^4) two-body contributions with uncertainties.

O(delta^2) + O(delta^4) = total two-body result [O(q^4)]
So O(delta^4) = O(q^4) - O(q^3), which is the "Static" quantity in FFSlib.
"""

import math
import sys

sys.path.insert(0, "/home/alexander/OneDrive/DensityScattering/tools/package")

from FFSlib import get_literature_results, NoResult, nuclei


def _determine_precision(unc):
    """Determine decimal places from uncertainty (2 sig figs, trim trailing zero)."""
    mag = math.floor(math.log10(abs(unc)))
    n_dec = -(mag - 1)
    unc_int = int(round(unc * 10**n_dec))
    if unc_int % 10 == 0 and unc_int >= 10 and n_dec > 1:
        n_dec -= 1
    return max(0, n_dec)


def _format_split(val, precision=None, force_two_parens=False):
    """Format [mean, unc1, unc2] as mean(unc1)(unc2)."""
    if isinstance(val, NoResult):
        return "--"

    mean, unc1, unc2 = val[0], val[1], val[2]
    has_1 = unc1 > 0
    has_2 = unc2 > 0

    if not has_1 and not has_2:
        p = precision if precision is not None else 3
        if force_two_parens:
            return f"{mean:.{p}f}(0)(0)"
        return f"{mean:.{p}f}"

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


_ff_quantities = {"F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L"}


def _tda_precision(q, val, base_prec=2):
    if q in _ff_quantities and not isinstance(val, NoResult) and abs(val[0]) < 0.1:
        return 3
    return base_prec


d = get_literature_results()

quantities = [
    ("F^{(a)}_T-F^{(b)}_T", "F_T^(a)-F_T^(b)"),
    ("F^{(a)}_L-F^{(b)}_L", "F_L^(a)-F_L^(b)"),
    ("E2N", "E2N"),
    ("L2N", "L2N"),
]

authors = [("TDA", "TDA"), ("Lenkewitz", "lenke")]

# Verify O(q^3) + O(delta^4) = O(q^4) for all entries
failures = []
for author_label, author_key in authors:
    for nuc in nuclei:
        for qkey, qlabel in quantities:
            q3 = d[author_key][nuc][f"{qkey} O(q^3)"]
            q4 = d[author_key][nuc][f"{qkey} O(q^4)"]
            st = d[author_key][nuc][f"Static {qkey}"]
            if (
                isinstance(q3, NoResult)
                or isinstance(q4, NoResult)
                or isinstance(st, NoResult)
            ):
                continue
            if not math.isclose(q3[0] + st[0], q4[0], rel_tol=1e-10):
                failures.append(
                    f"  {author_label} {nuc} {qlabel}: "
                    f"q3+d4={q3[0] + st[0]:.10f} != q4={q4[0]:.10f}"
                )

w_src = 10
w_nuc = 8
w_qty = 18
w_fmt = 22
w_mean = 18
w_unc = 16

hdr = (
    f"{'Source':<{w_src}} {'Nucleus':<{w_nuc}} {'Quantity':<{w_qty}} "
    f"{'O(delta^4)':>{w_fmt}} {'mean':>{w_mean}} {'unc1':>{w_unc}} {'unc2':>{w_unc}}"
)
sep = "-" * len(hdr)

print("O(delta^4) two-body contributions")
print(sep)
print(hdr)
print(sep)

for author_label, author_key in authors:
    is_tda = author_key == "TDA"
    for nuc in nuclei:
        for qkey, qlabel in quantities:
            st = d[author_key][nuc][f"Static {qkey}"]
            if isinstance(st, NoResult):
                continue
            prec = _tda_precision(qkey, st) if is_tda else None
            formatted = _format_split(st, precision=prec, force_two_parens=is_tda)
            print(
                f"{author_label:<{w_src}} {nuc:<{w_nuc}} {qlabel:<{w_qty}} "
                f"{formatted:>{w_fmt}} {st[0]:>{w_mean}.10f} {st[1]:>{w_unc}.10f} {st[2]:>{w_unc}.10f}"
            )
        print(sep)

print()
print("Uncertainty key:")
print(
    "  TDA form factors:  unc1 = potential variation (theory spread),  unc2 = numerical precision (1.5%)"
)
print(
    "  TDA E2N/L2N:       unc1 = combined nuclear structure (theory + numeric in quadrature),  unc2 = 0"
)
print("  Lenkewitz:         unc1 = theory spread,  unc2 = numerical precision")

if failures:
    print("\nCHECK FAILURES (O(q^3) + O(delta^4) != O(q^4)):")
    for f in failures:
        print(f)
