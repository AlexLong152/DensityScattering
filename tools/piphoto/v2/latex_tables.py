# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import math
from FFSlib import get_literature_results
from FFSlib import NoResult, nuclei


def main():
    d = get_literature_results()

    onebody_rows = ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V", "E1N", "L1N"]
    twobody_rows = [
        "F^{(a)}_T-F^{(b)}_T O(q^3)",
        "F^{(a)}_L-F^{(b)}_L O(q^3)",
        "E2N O(q^3)",
        "L2N O(q^3)",
        "F^{(a)}_T-F^{(b)}_T O(q^4)",
        "F^{(a)}_L-F^{(b)}_L O(q^4)",
        "E2N O(q^4)",
        "L2N O(q^4)",
    ]

    for nuc in nuclei:
        if nuc == "6Li":
            authors = [("TDA Method", "TDA"), ("Braun", "braun")]
        else:
            authors = [
                ("TDA Method", "TDA"),
                ("Lenkewitz", "lenke"),
                ("Braun", "braun"),
            ]
        nuc_label = _nuc_labels[nuc]

        print(
            _latex_table(
                d, nuc, f"One-body results for {nuc_label}.", onebody_rows, authors
            )
        )
        print()
        print(
            _latex_table(
                d, nuc, f"Two-body results for {nuc_label}.", twobody_rows, authors
            )
        )
        print()


def _format_parens_latex(val, precision=None):
    """Format np.array([mean, unc]) as $1.551(78)$ for LaTeX.
    If precision is given, use that many decimal places.
    If None, determine decimal places from the uncertainty.
    Returns empty string for NoResult."""
    if isinstance(val, NoResult):
        return ""
    mean, unc = val
    if unc <= 0:
        p = precision if precision is not None else 3
        return f"${mean:.{p}f}$"
    if precision is not None:
        scale = 10**precision
        mean_r = round(mean, precision)
        unc_int = int(round(unc * scale))
        return f"${mean_r:.{precision}f}({unc_int})$"
    # Auto: determine decimal places from uncertainty
    mag = math.floor(math.log10(abs(unc)))
    n_dec = -(mag - 1)  # 2 sig figs of unc
    unc_r = round(unc, n_dec)
    unc_int = int(round(unc_r * 10**n_dec))
    if unc_int % 10 == 0 and unc_int >= 10:
        n_dec -= 1
        unc_r = round(unc, n_dec)
        unc_int = int(round(unc_r * 10**n_dec))
    n_dec = max(0, n_dec)
    mean_r = round(mean, n_dec)
    return f"${mean_r:.{n_dec}f}({unc_int})$"


# Map from dictionary keys to LaTeX row labels
_latex_labels = {
    "F_T^S+V": r"$F_T^{S+V}$",
    "F_T^S-V": r"$F_T^{S-V}$",
    "F_L^S+V": r"$F_L^{S+V}$",
    "F_L^S-V": r"$F_L^{S-V}$",
    "E1N": r"$E_{0+}^{1N}$",
    "L1N": r"$L_{0+}^{1N}$",
    "F^{(a)}_T-F^{(b)}_T O(q^3)": r"$F_T^{(a)}-F_T^{(b)}$ $\mathcal{O}(q^3)$",
    "F^{(a)}_L-F^{(b)}_L O(q^3)": r"$F_L^{(a)}-F_L^{(b)}$ $\mathcal{O}(q^3)$",
    "E2N O(q^3)": r"$E_{0+}^{2N}$ $\mathcal{O}(q^3)$",
    "L2N O(q^3)": r"$L_{0+}^{2N}$ $\mathcal{O}(q^3)$",
    "Static F^{(a)}_T-F^{(b)}_T": r"$F_T^{(a)}-F_T^{(b)}$ static",
    "Static F^{(a)}_L-F^{(b)}_L": r"$F_L^{(a)}-F_L^{(b)}$ static",
    "Static E2N": r"$E_{0+}^{2N}$ static",
    "Static L2N": r"$L_{0+}^{2N}$ static",
    "F^{(a)}_T-F^{(b)}_T O(q^4)": r"$F_T^{(a)}-F_T^{(b)}$ $\mathcal{O}(q^4)$",
    "F^{(a)}_L-F^{(b)}_L O(q^4)": r"$F_L^{(a)}-F_L^{(b)}$ $\mathcal{O}(q^4)$",
    "E2N O(q^4)": r"$E_{0+}^{2N}$ $\mathcal{O}(q^4)$",
    "L2N O(q^4)": r"$L_{0+}^{2N}$ $\mathcal{O}(q^4)$",
}

_nuc_labels = {
    "3H": r"${}^3\mathrm{H}$",
    "3He": r"${}^3\mathrm{He}$",
    "6Li": r"${}^6\mathrm{Li}$",
}


def _latex_table(d, nuc, caption, rows, authors):
    """Return a LaTeX table string for one nucleus.

    rows: list of quantity key strings
    authors: list of (label, author_key)
    """
    ncols = len(authors)
    col_spec = "|l||" + "c|" * ncols
    nuc_label = _nuc_labels.get(nuc, nuc)
    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")

    # Header row
    header = " & ".join([nuc_label] + [label for label, _ in authors]) + r" \\"
    lines.append(header)
    lines.append(r"\hline\hline")

    # Data rows
    for q in rows:
        label = _latex_labels.get(q, q)
        cells = [label]
        for _, author in authors:
            prec = 2 if author == "TDA" else None
            cells.append(_format_parens_latex(d[author][nuc][q], precision=prec))
        lines.append(" & ".join(cells) + r" \\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


def printLatex(d=None, showLenke=True):
    """Print LaTeX tables for all nuclei.

    Args:
        d: results dictionary (loaded automatically if None)
        showLenke: if False, omit the Lenkewitz column
    """
    if d is None:
        from FFSlib import get_literature_results

        d = get_literature_results()

    authors_full = [("TDA Method", "TDA"), ("Lenkewitz", "lenke"), ("Braun", "braun")]
    authors_no_lenke = [("TDA Method", "TDA"), ("Braun", "braun")]

    onebody_rows = ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V", "E1N", "L1N"]

    twobody_rows = [
        "F^{(a)}_T-F^{(b)}_T O(q^3)",
        "F^{(a)}_L-F^{(b)}_L O(q^3)",
        "E2N O(q^3)",
        "L2N O(q^3)",
        # "Static F^{(a)}_T-F^{(b)}_T",
        # "Static F^{(a)}_L-F^{(b)}_L",
        # "Static E2N",
        # "Static L2N",
        "F^{(a)}_T-F^{(b)}_T O(q^4)",
        "F^{(a)}_L-F^{(b)}_L O(q^4)",
        "E2N O(q^4)",
        "L2N O(q^4)",
    ]

    for nuc in nuclei:
        authors = authors_full if showLenke else authors_no_lenke
        nuc_label = _nuc_labels.get(nuc, nuc)

        print(
            _latex_table(
                d, nuc, f"One-body results for {nuc_label}.", onebody_rows, authors
            )
        )
        print()
        print(
            _latex_table(
                d, nuc, f"Two-body results for {nuc_label}.", twobody_rows, authors
            )
        )
        print()


if __name__ == "__main__":
    main()
