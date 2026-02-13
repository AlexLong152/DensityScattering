# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import math
from FFSlib import get_literature_results
from FFSlib import NoResult, nuclei


def main():
    d = get_literature_results()

    onebody_quantities = ["F_T^S+V", "F_T^S-V", "E1N", "F_L^S+V", "F_L^S-V", "L1N"]
    twobody_quantities = [
        "F^{(a)}_T-F^{(b)}_T",
        "F^{(a)}_L-F^{(b)}_L",
        "E2N",
        "L2N",
    ]
    twobody_row_groups = [
        (r"$\mathcal{O}(q^3)$", " O(q^3)"),
        (r"$\mathcal{O}(q^4)$", " O(q^4)"),
    ]

    for nuc in nuclei:
        if nuc == "6Li":
            authors = [("TDA Method", "TDA"), ("Braun", "braun")]
            authors_2b = [("TDA Method", "TDA"), (r"Braun$\times\sqrt{2}$", "braun")]
            author_scales_2b = {"braun": math.sqrt(2)}
        else:
            authors = [
                ("TDA Method", "TDA"),
                ("Lenkewitz", "lenke"),
                ("Braun", "braun"),
            ]
            authors_2b = authors
            author_scales_2b = None
        nuc_label = _nuc_labels[nuc]

        print(
            _latex_table(
                d,
                nuc,
                rf"One-body results for {nuc_label}. Form factors $F_{{T/L}}$ are unitless. $E_{{0+}}^{{1N}}$ and $L_{{0+}}^{{1N}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$.",
                onebody_quantities,
                authors,
            )
        )
        print()
        print(
            _latex_table(
                d,
                nuc,
                rf"Two-body results for {nuc_label}. $F_{{T/L}}^{{(a)}}-F_{{T/L}}^{{(b)}}$ in units of $[\mathrm{{fm}}^{{-1}}]$. $E_{{0+}}^{{2N}}$ and $L_{{0+}}^{{2N}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$. $\mathcal{{O}}(q^4)$ results omit recoil contributions.",
                twobody_quantities,
                authors_2b,
                row_groups=twobody_row_groups,
                author_scales=author_scales_2b,
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
    "F^{(a)}_T-F^{(b)}_T O(q^3)": r"$F_T^{(a)}-F_T^{(b)}$",
    "F^{(a)}_L-F^{(b)}_L O(q^3)": r"$F_L^{(a)}-F_L^{(b)}$",
    "E2N O(q^3)": r"$E_{0+}^{2N}$",
    "L2N O(q^3)": r"$L_{0+}^{2N}$",
    "Static F^{(a)}_T-F^{(b)}_T": r"$F_T^{(a)}-F_T^{(b)}$ static",
    "Static F^{(a)}_L-F^{(b)}_L": r"$F_L^{(a)}-F_L^{(b)}$ static",
    "Static E2N": r"$E_{0+}^{2N}$ static",
    "Static L2N": r"$L_{0+}^{2N}$ static",
    "F^{(a)}_T-F^{(b)}_T O(q^4)": r"$F_T^{(a)}-F_T^{(b)}$",
    "F^{(a)}_L-F^{(b)}_L O(q^4)": r"$F_L^{(a)}-F_L^{(b)}$",
    "E2N O(q^4)": r"$E_{0+}^{2N}$",
    "L2N O(q^4)": r"$L_{0+}^{2N}$",
}

_nuc_labels = {
    "3H": r"${}^3\mathrm{H}$",
    "3He": r"${}^3\mathrm{He}$",
    "6Li": r"${}^6\mathrm{Li}$",
}


def _latex_table(d, nuc, caption, quantities, authors, row_groups=None, author_scales=None):
    r"""Return a LaTeX table string for one nucleus (transposed layout).

    Columns are quantities (form factors), rows are sources (authors).

    quantities: list of quantity key strings (become columns).
        For grouped tables, these are base keys; the full data key is
        base_key + group_suffix.
    authors: list of (label, author_key) — become rows
    row_groups: optional list of (group_label, key_suffix) for grouped rows
        e.g. [(r"$\mathcal{O}(q^3)$", " O(q^3)"), ...]
        If None, quantities are used as data keys directly.
    """
    nuc_label = _nuc_labels.get(nuc, nuc)
    ncols = len(quantities)
    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")

    # Determine column labels from _latex_labels
    col_labels = []
    for q in quantities:
        if row_groups:
            lookup_key = q + row_groups[0][1]
        else:
            lookup_key = q
        col_labels.append(_latex_labels.get(lookup_key, q))

    if row_groups:
        col_spec = "|c|l||" + "c|" * ncols
        lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")
        merged = rf"\multicolumn{{2}}{{|c||}}{{{nuc_label}}}"
        header = " & ".join([merged] + col_labels) + r" \\"
    else:
        col_spec = "|l||" + "c|" * ncols
        lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")
        header = " & ".join([nuc_label] + col_labels) + r" \\"

    lines.append(header)
    lines.append(r"\thickline{0.15em}")

    # Data rows: each row is one author, columns are quantities
    if row_groups:
        cline_end = 2 + ncols
        for gi, (group_label, suffix) in enumerate(row_groups):
            # Skip authors where all values in this group are NoResult
            active = [
                (label, key)
                for label, key in authors
                if any(
                    not isinstance(d[key][nuc][q + suffix], NoResult)
                    for q in quantities
                )
            ]
            for j, (author_label, author_key) in enumerate(active):
                if j == 0:
                    group_cell = rf"\multirow{{{len(active)}}}{{*}}{{{group_label}}}"
                else:
                    group_cell = ""
                cells = [group_cell, author_label]
                for q in quantities:
                    data_key = q + suffix
                    prec = 2 if author_key == "TDA" else None
                    val = d[author_key][nuc][data_key]
                    if author_scales and author_key in author_scales:
                        if not isinstance(val, NoResult):
                            val = val * author_scales[author_key]
                    cells.append(
                        _format_parens_latex(val, precision=prec)
                    )
                lines.append(" & ".join(cells) + r" \\")
                if j < len(active) - 1:
                    lines.append(rf"\cline{{2-{cline_end}}}")
            # Between groups: double hline; after last group: single hline
            if gi < len(row_groups) - 1:
                lines.append(r"\hline\hline")
            else:
                lines.append(r"\hline")
    else:
        for author_label, author_key in authors:
            cells = [author_label]
            for q in quantities:
                prec = 2 if author_key == "TDA" else None
                val = d[author_key][nuc][q]
                if author_scales and author_key in author_scales:
                    if not isinstance(val, NoResult):
                        val = val * author_scales[author_key]
                cells.append(
                    _format_parens_latex(val, precision=prec)
                )
            lines.append(" & ".join(cells) + r" \\")
            lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(r"\end{table}")

    return "\n".join(lines).replace("$-16(2)$", "$-16.0(20)$")


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

    onebody_quantities = ["F_T^S+V", "F_T^S-V", "E1N", "F_L^S+V", "F_L^S-V", "L1N"]

    twobody_quantities = [
        "F^{(a)}_T-F^{(b)}_T",
        "F^{(a)}_L-F^{(b)}_L",
        "E2N",
        "L2N",
    ]
    twobody_row_groups = [
        (r"$\mathcal{O}(q^3)$", " O(q^3)"),
        (r"$\mathcal{O}(q^4)$", " O(q^4)"),
    ]

    for nuc in nuclei:
        authors = authors_full if showLenke else authors_no_lenke
        if nuc == "6Li":
            authors_2b = [
                (r"Braun$\times\sqrt{2}$", key) if key == "braun" else (label, key)
                for label, key in authors
            ]
            author_scales_2b = {"braun": math.sqrt(2)}
        else:
            authors_2b = authors
            author_scales_2b = None
        nuc_label = _nuc_labels.get(nuc, nuc)

        print(
            _latex_table(
                d,
                nuc,
                rf"One-body results for {nuc_label}. Form factors $F_{{T/L}}$ are unitless. $E_{{0+}}^{{1N}}$ and $L_{{0+}}^{{1N}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$.",
                onebody_quantities,
                authors,
            )
        )
        print()
        print(
            _latex_table(
                d,
                nuc,
                rf"Two-body results for {nuc_label}. $F_{{T/L}}^{{(a)}}-F_{{T/L}}^{{(b)}}$ in units of $[\mathrm{{fm}}^{{-1}}]$. $E_{{0+}}^{{2N}}$ and $L_{{0+}}^{{2N}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$. $\mathcal{{O}}(q^4)$ results omit recoil contributions.",
                twobody_quantities,
                authors_2b,
                row_groups=twobody_row_groups,
                author_scales=author_scales_2b,
            )
        )
        print()


if __name__ == "__main__":
    main()
