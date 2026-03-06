import math
import numpy as np
from DirectCalc import (
    aPlusVals,
    aMinusVals,
    aPlus,
    aMinus,
    mpiArr,
    MAve,
)
import importlib
TwoN = importlib.import_module("2N")


charge = np.array([-1, 0, 1])
chargeLabels = ["pi-", "pi0", "pi+"]
chargeLabelsLatex = [r"$\pi^-$", r"$\pi^0$", r"$\pi^+$"]

nuclei = {
    "3He": {"Z": 2, "N": 1},
    "4He": {"Z": 2, "N": 2},
    "6Li": {"Z": 3, "N": 3},
}

_nuc_labels = {
    "3He": r"${}^3\mathrm{He}$",
    "4He": r"${}^4\mathrm{He}$",
    "6Li": r"${}^6\mathrm{Li}$",
}


def get_1N(nucleus):
    """
    Compute 1N scattering lengths from the direct calculation.

    Returns
    -------
    numpy.ndarray
        Shape (3, 2) with columns [central, uncertainty] for pi-, pi0, pi+.
    """
    Z = nuclei[nucleus]["Z"]
    N = nuclei[nucleus]["N"]
    A = Z + N

    out = np.zeros((3, 2))
    for ic, Qpi in enumerate(charge):
        mpi_ic = mpiArr[ic]
        numer = 1 + mpi_ic / MAve
        denom = 1 + mpi_ic / (A * MAve)
        kinFactor = numer / denom
        vals = []
        for ap in aPlusVals:
            for am in aMinusVals:
                vals.append(kinFactor * (A * ap - Qpi * (Z - N) * am))
        central = kinFactor * (A * aPlus - Qpi * (Z - N) * aMinus)
        unc = (max(vals) - min(vals)) / 2.0
        out[ic][0] = central
        out[ic][1] = unc
    return out


def get_2N(nucleus):
    """
    Compute 2N scattering lengths from density files.

    Returns
    -------
    numpy.ndarray
        Shape (3, 2) with columns [mean, uncertainty] for pi-, pi0, pi+.
    """
    arr = TwoN.get_scattering_length_array(nucleus)
    return TwoN.getMeanUncer(arr)


def _format_parens_latex(val, precision=None):
    """Format [mean, unc] as $1.551(78)$ for LaTeX.
    If precision is given, use that many decimal places.
    If None, determine decimal places from the uncertainty."""
    mean, unc = val
    if unc <= 0:
        p = precision if precision is not None else 3
        return f"${mean:.{p}f}$"
    if precision is not None:
        scale = 10**precision
        mean_r = round(mean, precision)
        unc_int = int(round(unc * scale))
        return f"${mean_r:.{precision}f}({unc_int})$"
    mag = math.floor(math.log10(abs(unc)))
    n_dec = max(0, -(mag - 1))  # 2 sig figs of uncertainty
    unc_int = int(round(unc * 10**n_dec))
    mean_r = round(mean, n_dec)
    return f"${mean_r:.{n_dec}f}({unc_int})$"


def _align_column(formatted):
    r"""Add \phantom padding so decimal points align within a column.

    Pads both the integer part (left of decimal) and the fractional+suffix
    part (right of decimal) so that entries with different numbers of
    decimal places still line up at the decimal point.

    Parameters
    ----------
    formatted : list of str
        LaTeX-formatted strings like ``$88.5(54)$`` or ``$-24.37(78)$``.

    Returns
    -------
    list of str
        Padded versions of the input strings.
    """
    int_parts = []
    frac_parts = []  # everything after the decimal point, before closing $
    for s in formatted:
        inner = s.strip("$")
        if "." in inner:
            ip, fp = inner.split(".", 1)
        else:
            ip, fp = inner, ""
        int_parts.append(ip)
        frac_parts.append(fp)

    max_int_len = max(len(p) for p in int_parts)
    max_frac_len = max(len(p) for p in frac_parts)

    out = []
    for s, ip, fp in zip(formatted, int_parts, frac_parts):
        left_pad = max_int_len - len(ip)
        right_pad = max_frac_len - len(fp)
        result = "$"
        if left_pad > 0:
            result += r"\phantom{" + "0" * left_pad + "}"
        # strip the leading $
        result += s[1:-1]
        if right_pad > 0:
            result += r"\phantom{" + "0" * right_pad + "}"
        result += "$"
        out.append(result)
    return out


def latex_table():
    r"""Return a LaTeX table of pion-nucleus scattering lengths.

    Rows are grouped by nucleus, with sub-rows for each pion charge.
    Columns: Nucleus | pion | a^{1N} | a^{2N} | a^{1N+2N}.
    """
    col_labels = [
        r"$a^{1\mathrm{N}}$",
        r"$a^{2\mathrm{N}}$",
        r"$a=a^{1\mathrm{N}}+a^{2\mathrm{N}}$",
    ]
    ncols = len(col_labels)

    # Pre-compute all values
    nuc_list = list(nuclei.keys())
    col_vals = [[], [], []]  # collect [mean, unc] per data column

    for nuc in nuc_list:
        a1N = get_1N(nuc)
        a2N = get_2N(nuc)
        for ic in range(3):
            total_mean = a1N[ic, 0] + a2N[ic, 0]
            total_unc = math.sqrt(a1N[ic, 1] ** 2 + a2N[ic, 1] ** 2)
            col_vals[0].append(a1N[ic])
            col_vals[1].append(a2N[ic])
            col_vals[2].append([total_mean, total_unc])

    # Format with auto-precision (2 sig figs of uncertainty per entry)
    col_strs = [[], [], []]
    for k in range(3):
        for v in col_vals[k]:
            col_strs[k].append(_format_parens_latex(v))

    # Align each column at the decimal point
    for k in range(3):
        col_strs[k] = _align_column(col_strs[k])

    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    col_spec = "|l|l||" + "r|" * ncols
    lines.append(r"{\setlength{\extrarowheight}{3pt}")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}\hline")

    merged = r"\multicolumn{2}{|c||}{\rule[-1.1ex]{0pt}{2.8ex}}"
    header = " & ".join([merged] + col_labels) + r" \\"
    lines.append(header)
    lines.append(r"\thickline{0.15em}")

    row_idx = 0
    for ni, nuc in enumerate(nuc_list):
        nuc_label = _nuc_labels[nuc]

        for ic in range(3):
            if ic == 0:
                group_cell = rf"\multirow{{3}}{{*}}{{{nuc_label}}}"
            else:
                group_cell = ""

            cells = [
                group_cell,
                chargeLabelsLatex[ic],
                col_strs[0][row_idx],
                col_strs[1][row_idx],
                col_strs[2][row_idx],
            ]
            lines.append(" & ".join(cells) + r" \\")
            if ic < 2:
                lines.append(rf"\cline{{2-{2 + ncols}}}")
            row_idx += 1

        if ni < len(nuc_list) - 1:
            lines.append(r"\hline\hline")
        else:
            lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append("}")
    lines.append(
        r"\caption{Pion-nucleus scattering lengths in units of $[10^{-3}/m_{\pi}]$.}"
    )
    lines.append(r"\end{table}")

    return "\n".join(lines)


if __name__ == "__main__":
    for nuc in ["3He", "4He", "6Li"]:
        a1N = get_1N(nuc)
        a2N = get_2N(nuc)
        total_central = a1N[:, 0] + a2N[:, 0]
        total_unc = np.sqrt(a1N[:, 1] ** 2 + a2N[:, 1] ** 2)

        print(f"{nuc}:")
        for ic, label in enumerate(chargeLabels):
            print(
                f"  {label}:  a1N = {a1N[ic][0]:8.2f} +/- {a1N[ic][1]:.2f}"
                f"   a2N = {a2N[ic][0]:8.2f} +/- {a2N[ic][1]:.2f}"
                f"   total = {total_central[ic]:8.2f} +/- {total_unc[ic]:.2f}"
            )
        print()

    print()
    print(latex_table())
