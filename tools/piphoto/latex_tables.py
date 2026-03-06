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
    twobody_quantities = ["F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L", "E2N", "L2N"]
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
                rf"One-body results for {nuc_label}. Form factors $F_{{T/L}}$ are unitless. $E_{{0+}}^{{1\mathrm{{N}}}}$ and $L_{{0+}}^{{1\mathrm{{N}}}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$.",
                onebody_quantities,
                authors,
            )
        )
        print()
        print(
            _latex_table(
                d,
                nuc,
                rf"Two-body results for {nuc_label}. $F_{{T/L}}^{{(a)}}-F_{{T/L}}^{{(b)}}$ in units of $[\mathrm{{fm}}^{{-1}}]$. $E_{{0+}}^{{2\mathrm{{N}}}}$ and $L_{{0+}}^{{2\mathrm{{N}}}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$. $\mathcal{{O}}(q^4)$ results omit recoil contributions.",
                twobody_quantities,
                authors_2b,
                row_groups=twobody_row_groups,
                author_scales=author_scales_2b,
            )
        )
        print()

    print(50 * "#")
    print(
        _summary_table(
            d,
            r"Summary of pion photoproduction results. $E_{0+}$ and $L_{0+}$ in units of $[10^{-3}/m_{\pi^+}]$. $a_0$ in units of $[10^{-6}/m_{\pi^+}^2]$. $E_{0+} = E_{0+}^{1\mathrm{N}} + E_{0+}^{2\mathrm{N}}$, $L_{0+} = L_{0+}^{1\mathrm{N}} + L_{0+}^{2\mathrm{N}}$. All $1\mathrm{N}$ values to $\mathcal{O}(q^3)$. TDA and Lenkewitz $2\mathrm{N}$ results at $\mathcal{O}(q^4)$; Braun at $\mathcal{O}(q^3)$.",
            twobody_orders={"TDA": " O(q^4)", "lenke": " O(q^4)", "braun": " O(q^3)"},
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
    "E1N": r"$E_{0+}^{1\mathrm{N}}$",
    "L1N": r"$L_{0+}^{1\mathrm{N}}$",
    "F^{(a)}_T-F^{(b)}_T O(q^3)": r"$F_T^{(a)}-F_T^{(b)}$",
    "F^{(a)}_L-F^{(b)}_L O(q^3)": r"$F_L^{(a)}-F_L^{(b)}$",
    "E2N O(q^3)": r"$E_{0+}^{2\mathrm{N}}$",
    "L2N O(q^3)": r"$L_{0+}^{2\mathrm{N}}$",
    "Static F^{(a)}_T-F^{(b)}_T": r"$F_T^{(a)}-F_T^{(b)}$ static",
    "Static F^{(a)}_L-F^{(b)}_L": r"$F_L^{(a)}-F_L^{(b)}$ static",
    "Static E2N": r"$E_{0+}^{2\mathrm{N}}$ static",
    "Static L2N": r"$L_{0+}^{2\mathrm{N}}$ static",
    "F^{(a)}_T-F^{(b)}_T O(q^4)": r"$F_T^{(a)}-F_T^{(b)}$",
    "F^{(a)}_L-F^{(b)}_L O(q^4)": r"$F_L^{(a)}-F_L^{(b)}$",
    "E2N O(q^4)": r"$E_{0+}^{2\mathrm{N}}$",
    "L2N O(q^4)": r"$L_{0+}^{2\mathrm{N}}$",
}

_nuc_labels = {
    "3H": r"${}^3\mathrm{H}$",
    "3He": r"${}^3\mathrm{He}$",
    "6Li": r"${}^6\mathrm{Li}$",
}

_spins = {"3H": r"$1/2$", "3He": r"$1/2$", "6Li": r"$1$"}
_spin_values = {"3H": 0.5, "3He": 0.5, "6Li": 1.0}


def _latex_table(
    d, nuc, caption, quantities, authors, row_groups=None, author_scales=None
):
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
        col_spec = "|l|l||" + "l|" * ncols
        lines.append(r"{\setlength{\extrarowheight}{3pt}")
        lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")
        merged = (
            rf"\multicolumn{{2}}{{|c||}}{{{nuc_label} \rule[-1.1ex]{{0pt}}{{2.8ex}}}}"
        )
        header = " & ".join([merged] + col_labels) + r" \\"
    else:
        col_spec = "|l||" + "l|" * ncols
        lines.append(r"{\setlength{\extrarowheight}{3pt}")
        lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")
        header = (
            " & ".join([rf"{nuc_label} \rule[-1.1ex]{{0pt}}{{2.8ex}}"] + col_labels)
            + r" \\"
        )

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
                    cells.append(_format_parens_latex(val, precision=prec))
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
                cells.append(_format_parens_latex(val, precision=prec))
            lines.append(" & ".join(cells) + r" \\")
            lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append("}")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(r"\end{table}")

    return "\n".join(lines).replace("$-16(2)$", "$-16.0(20)$")


def _summary_table(d, caption, twobody_orders=None):
    r"""Return a LaTeX summary table combining 1N and 2N results for all nuclei.

    Rows grouped by nucleus (left column, multirow), with author rows within each.
    Columns: E_{0+}^{1N}, E_{0+}^{2N}, E_{0+}, L_{0+}^{1N}, L_{0+}^{2N}, L_{0+}, a_0.
    E_{0+} = E_{0+}^{1N} + E_{0+}^{2N}; similarly for L.

    twobody_orders: dict mapping author_key to order suffix, e.g.
        {"TDA": " O(q^4)", "lenke": " O(q^4)", "braun": " O(q^3)"}.
        Defaults to " O(q^3)" for all authors.
    """
    if twobody_orders is None:
        twobody_orders = {}
    default_order = " O(q^3)"
    col_labels = [
        r"$E_{0+}^{1\mathrm{N}}$",
        r"$E_{0+}^{2\mathrm{N}}$",
        r"$E_{0+}$",
        r"$L_{0+}^{1\mathrm{N}}$",
        r"$L_{0+}^{2\mathrm{N}}$",
        r"$L_{0+}$",
        r"$a_0$",
    ]
    ncols = len(col_labels)

    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    col_spec = "|l|l||" + "l|" * ncols
    lines.append(r"{\setlength{\extrarowheight}{3pt}")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")

    merged = r"\multicolumn{2}{|c||}{\rule[-1.1ex]{0pt}{2.8ex}}"
    header = " & ".join([merged] + col_labels) + r" \\"
    lines.append(header)
    lines.append(r"\thickline{0.15em}")

    cline_end = 2 + ncols

    for ni, nuc in enumerate(nuclei):
        nuc_label = _nuc_labels[nuc]
        if nuc == "6Li":
            authors = [("TDA", "TDA"), ("Braun*", "braun")]
            author_scales_2b = {"braun": math.sqrt(2)}
        else:
            authors = [("TDA", "TDA"), ("Lenkewitz", "lenke"), ("Braun", "braun")]
            author_scales_2b = None

        # Filter to authors that have at least one non-NoResult value
        active = [
            (label, key)
            for label, key in authors
            if any(
                not isinstance(d[key][nuc][q], NoResult)
                for q in [
                    "E1N",
                    f"E2N{twobody_orders.get(key, default_order)}",
                    "L1N",
                    f"L2N{twobody_orders.get(key, default_order)}",
                ]
            )
        ]

        for j, (author_label, author_key) in enumerate(active):
            if j == 0:
                group_cell = rf"\multirow{{{len(active)}}}{{*}}{{{nuc_label}}}"
            else:
                group_cell = ""

            cells = [group_cell, author_label]
            order = twobody_orders.get(author_key, default_order)

            # E1N
            e1n = d[author_key][nuc]["E1N"]
            cells.append(
                ""
                if isinstance(e1n, NoResult)
                else _format_parens_latex([e1n[0], 0], precision=2)
            )

            # E2N
            e2n = d[author_key][nuc][f"E2N{order}"]
            if author_scales_2b and author_key in author_scales_2b:
                if not isinstance(e2n, NoResult):
                    e2n = e2n * author_scales_2b[author_key]
            cells.append(
                ""
                if isinstance(e2n, NoResult)
                else _format_parens_latex([e2n[0], 0], precision=2)
            )

            # E0+ = E1N + E2N
            if not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult):
                e_total = [e1n[0] + e2n[0], math.sqrt(e1n[1] ** 2 + e2n[1] ** 2)]
                cells.append(_format_parens_latex(e_total, precision=2))
            else:
                cells.append("")

            # L1N
            l1n = d[author_key][nuc]["L1N"]
            cells.append(
                ""
                if isinstance(l1n, NoResult)
                else _format_parens_latex([l1n[0], 0], precision=2)
            )

            # L2N
            l2n = d[author_key][nuc][f"L2N{order}"]
            cells.append(
                ""
                if isinstance(l2n, NoResult)
                else _format_parens_latex([l2n[0], 0], precision=2)
            )

            # L0+ = L1N + L2N
            if not isinstance(l1n, NoResult) and not isinstance(l2n, NoResult):
                l_total = [l1n[0] + l2n[0], math.sqrt(l1n[1] ** 2 + l2n[1] ** 2)]
                cells.append(_format_parens_latex(l_total, precision=2))
            else:
                cells.append("")

            # a_0
            if not isinstance(e1n, NoResult) and not isinstance(e2n, NoResult):
                a0, a0_unc = a0Value(e_total[0], e_total[1], _spin_values[nuc])
                cells.append(_format_parens_latex([a0, a0_unc], precision=2))
            else:
                cells.append("")

            lines.append(" & ".join(cells) + r" \\")
            if j < len(active) - 1:
                lines.append(rf"\cline{{2-{cline_end}}}")

        if ni < len(nuclei) - 1:
            lines.append(r"\hline\hline")
        else:
            lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append("}")
    note = r" The value of $E_{0+}^{2\mathrm{N}}$ for ${}^6\mathrm{Li}$ in the work by Braun has been multiplied by $\sqrt{2}$ in accordance with discussion."
    lines.append(r"\caption{" + caption + note + "}")
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

    onebody_quantities = ["F_T^S+V", "F_T^S-V", "E1N", "F_L^S+V", "F_L^S-V", "L1N"]

    twobody_quantities = ["F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L", "E2N", "L2N"]
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
                rf"One-body results for {nuc_label}. Form factors $F_{{T/L}}$ are unitless. $E_{{0+}}^{{1\mathrm{{N}}}}$ and $L_{{0+}}^{{1\mathrm{{N}}}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$.",
                onebody_quantities,
                authors,
            )
        )
        print()
        print(
            _latex_table(
                d,
                nuc,
                rf"Two-body results for {nuc_label}. $F_{{T/L}}^{{(a)}}-F_{{T/L}}^{{(b)}}$ in units of $[\mathrm{{fm}}^{{-1}}]$. $E_{{0+}}^{{2\mathrm{{N}}}}$ and $L_{{0+}}^{{2\mathrm{{N}}}}$ in units of $[10^{{-3}}/m_{{\pi^+}}]$. $\mathcal{{O}}(q^4)$ results omit recoil contributions.",
                twobody_quantities,
                authors_2b,
                row_groups=twobody_row_groups,
                author_scales=author_scales_2b,
            )
        )
        print()


def a0Value(E0, uncer, j):
    """
    Given E0=E1N+E2N, and the uncertainty of E0, along with the spin j of the nucleus,
    calculates the scattering length, and the uncertainty of the scattering length using error propagation.
    """
    if isinstance(j, str):
        j = float(j)
    a0 = (4 / 3) * j * (j + 1) * (E0**2)
    a0_uncer = (4 / 3) * j * (j + 1) * 2 * abs(E0) * uncer
    return a0, a0_uncer


if __name__ == "__main__":
    main()
