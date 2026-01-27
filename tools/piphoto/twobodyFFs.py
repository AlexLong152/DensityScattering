# -*- coding: utf-8 -*-
"""Two-body form factor analysis for pion photoproduction."""

import numpy as np
from os import listdir
from os.path import isfile, join, dirname, basename
import array_to_latex as a2l

from FFSLib import (
    He3,
    H3,
    He4,
    Li6,
    percent,
    getFormFactorsFromFile,
    getScatMatFromFile,
    calculateSpread,
    _format_parens,
    _format_value,
    _get_files_from_folder,
)

import math

breakLine = "\n" + 150 * "%"


def _auto_precision(mean, max_sig_figs=4):
    """Compute decimal precision for at most max_sig_figs significant digits."""
    abs_mean = abs(mean)
    if abs_mean < 1e-10:
        return 3
    int_digits = max(1, int(math.floor(math.log10(abs_mean))) + 1)
    return max(0, min(3, max_sig_figs - int_digits))


# ============================================================================
# Two-body file pairing (Odelta2 vs Odelta4)
# ============================================================================


def getCombinedResults(path):
    """Extract [F_T, F_L, E_0+, L_0+] from output file."""
    ffs = getFormFactorsFromFile(path)
    scatMat = getScatMatFromFile(path)
    return np.concatenate([ffs, scatMat])


def getOrderedPaths(path):
    """Get (path_o2, path_o4) for a given Odelta2 or Odelta4 file."""
    file = basename(path)
    folder = dirname(path)

    if "twobody" not in file:
        raise ValueError(f"Expected twobody file, got: {file}")

    # Determine order and find matching file
    is_o2 = "Odelta2" in file
    is_o4 = "Odelta4" in file

    if not (is_o2 or is_o4):
        raise ValueError(f"File doesn't contain Odelta2 or Odelta4: {file}")

    target_order = "Odelta4" if is_o2 else "Odelta2"
    match_file = file.replace("Odelta2" if is_o2 else "Odelta4", target_order)

    # Find matching file
    matches = [
        f for f in listdir(folder) if isfile(join(folder, f)) and match_file in f
    ]

    if len(matches) != 1:
        raise AssertionError(
            f"Expected 1 match, found {len(matches)}:\n" + "\n".join(matches)
        )

    path_match = join(folder, matches[0])

    # Return in order: (Odelta2, Odelta4)
    return (path, path_match) if is_o2 else (path_match, path)


def calculateStaticContribution(path, extraction_func=getFormFactorsFromFile):
    """Calculate static contribution: Odelta4 - Odelta2.

    Returns (static_contrib, o2_result, o4_result).
    """
    path_o2, path_o4 = getOrderedPaths(path)
    o2_result = extraction_func(path_o2)
    o4_result = extraction_func(path_o4)
    return o4_result - o2_result, o2_result, o4_result


# ============================================================================
# Statistical analysis (two-body)
# ============================================================================


def analyzeTwoBodyFiles(files, extraction_func=getCombinedResults):
    """Analyze two-body files for Odelta2, Odelta4, and static contributions.

    Returns dict with 'o2', 'o4', 'static' keys, each with 'mean' and 'spread'.
    """
    results = {"o2": [], "o4": [], "static": []}

    for path in files:
        static, o2, o4 = calculateStaticContribution(path, extraction_func)
        results["o2"].append(o2)
        results["o4"].append(o4)
        results["static"].append(static)

    # Convert to arrays and calculate statistics
    output = {}
    for key, vals in results.items():
        arr = np.array(vals)
        spread = calculateSpread(arr)
        mean = np.mean(arr, axis=0)
        spread = np.sqrt(((percent * mean / 100) ** 2) + (spread**2))
        output[key] = {"mean": mean, "spread": spread}

    return output


# ============================================================================
# Output formatting (two-body)
# ============================================================================


def printTwoBodyTable(results, labels, precision=3):
    """Print formatted table of two-body results."""
    col_width = 25
    col_diff = 35
    label_width = max(len(label + "   [fm^-1]") for label in labels[:2])
    total_width = label_width + 2 * col_width + col_diff + 15

    full_labels = [
        labels[0] + "   [fm^-1]",
        labels[1] + "   [fm^-1]",
        labels[2] + "  [10^-3/m_π]",
        labels[3] + "  [10^-3/m_π]",
    ]

    print("-" * total_width)
    print(
        f"{'Two Body Results':<{label_width}} | "
        f"{'Oδ2 (mean ± σ)':^{col_width}} | "
        f"{'Oδ4 (mean ± σ)':^{col_width}} | "
        f"{'Diff (mean ± σ),  Oδ2 + Diff = Oδ4':^{col_diff}}"
    )
    print("-" * total_width)

    o2, o4, static = results["o2"], results["o4"], results["static"]

    for i, label in enumerate(full_labels):
        if i == 2:
            print("-" * total_width)
        print(
            f"{label:<{label_width}} | "
            f"{_format_value(o2['mean'][i], o2['spread'][i], col_width, precision)} | "
            f"{_format_value(o4['mean'][i], o4['spread'][i], col_width, precision)} | "
            f"{_format_value(static['mean'][i], static['spread'][i], col_diff, precision)}"
        )

    print("-" * total_width)
    print()


def createLatexTable(results, nuc_name):
    """Create LaTeX table for two-body results."""
    header_map = {
        "Lithium 6": r"\LiS",
        "Helium 4": r"\HeF",
        "Helium 3": r"\HeT",
        "Hydrogen 3": r"\HThree",
    }

    # Extract data and format as latex strings
    data_cols = [results[key] for key in ["o2", "o4", "static"]]
    out = np.empty((4, 3), dtype="U50")
    for i, data in enumerate(data_cols):
        out[:, i] = [_format_parens(m, u) for m, u in zip(data["mean"], data["spread"])]

    # Build full table with headers
    row_labels = [
        r"\hline$F_T  \;\;\;[\mathrm{fm}^{-1}]$",
        r"\hline$F_L  \;\;\;[\mathrm{fm}^{-1}]$",
        r"\hline$E^{2N}_{0+}\ [10^{-3}/m_\pi]$",
        r"\hline$L^{2N}_{0+}\ [10^{-3}/m_\pi]$",
    ]
    column_labels = [
        r"$\calO(\delta^2)$",
        r"$\calO(\delta^4)$ ",
        r"Static contribution to $\calO(\delta^4)$",
    ]

    out2 = np.empty((5, 4), dtype="U50")
    out2[:, :] = ""
    out2[0, 0] = "$" + header_map.get(nuc_name, nuc_name) + "$"
    out2[1:, 0] = row_labels
    out2[0, 1:] = column_labels
    out2[1:, 1:] = out

    latex_str = a2l.to_ltx(out2, arraytype="tabular", print_out=False)
    replaceStr = (
        r"\begin{table}[H]"
        + "\n"
        + r"\centering"
        + "\n"
        + r"\begin{tabular}{|l||c|c|c|}\hline"
    )
    latex_str = latex_str.replace(r"\begin{tabular}", replaceStr)
    latex_str = latex_str.replace(r"\end{tabular}", r"\hline\end{tabular}")
    latex_str = r"\renewcommand{\arraystretch}{1.4}" + "\n" + latex_str

    latex_str += "\n" + r"\end{table}"

    latex_str = latex_str.replace(
        "\n" + r"\hline\end{tabular}", r"\\" + "\n" + r"\hline\end{tabular}"
    )
    return latex_str


# ============================================================================
# High-level interface (two-body)
# ============================================================================


def createCombinedTwoBodyLatexTable(all_results):
    """Create combined LaTeX table of two-body Odelta2/Odelta4 results for all nuclei.

    all_results: list of (nuc_name, results_dict) excluding Helium 4.
    Table has 7 columns (1 label + 2 per nucleus) and 6 rows
    (nucleus names, order labels, F_T, F_L, E_{0+}, L_{0+}).
    """
    header_map = {"Lithium 6": r"\LiS", "Helium 3": r"\HeT", "Hydrogen 3": r"\HThree"}

    row_labels = [
        r"$F_T  \;\;\;[\mathrm{fm}^{-1}]$",
        r"$F_L  \;\;\;[\mathrm{fm}^{-1}]$",
        r"$E^{2N}_{0+}\ [10^{-3}/m_\pi]$",
        r"$L^{2N}_{0+}\ [10^{-3}/m_\pi]$",
    ]

    n_nuc = len(all_results)
    lines = []
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")

    # Column spec: |l|| then c|c|| repeated per nucleus (last pair ends with |)
    col_pairs = "c|c||" * (n_nuc - 1) + "c|c|"
    lines.append(r"\begin{tabular}{|l||" + col_pairs + r"}\hline")

    # Row 1: nucleus names spanning 2 columns each
    nuc_headers = []
    for i, (nuc_name, _) in enumerate(all_results):
        latex_name = header_map.get(nuc_name, nuc_name)
        sep = "||" if i < n_nuc - 1 else "|"
        nuc_headers.append(r"\multicolumn{2}{c" + sep + "}{$" + latex_name + "$}")
    lines.append(" & " + " & ".join(nuc_headers) + r" \\")
    lines.append(r"\cline{2-7}")

    # Row 2: Odelta2 / Odelta4 labels
    order_cells = []
    for _ in all_results:
        order_cells.append(r"$\calO(\delta^2)$")
        order_cells.append(r"$\calO(\delta^4)$")
    lines.append(" & " + " & ".join(order_cells) + r" \\")
    lines.append(r"\hline")

    # Rows 3-6: data
    for i, label in enumerate(row_labels):
        cells = [label]
        for _, results in all_results:
            o2_m, o2_s = results["o2"]["mean"][i], results["o2"]["spread"][i]
            o4_m, o4_s = results["o4"]["mean"][i], results["o4"]["spread"][i]
            cells.append(_format_parens(o2_m, o2_s, _auto_precision(o2_m)))
            cells.append(_format_parens(o4_m, o4_s, _auto_precision(o4_m)))
        lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


def createLiteratureTwoBodyLatexTable():
    """Create combined LaTeX table of Braun and Lenkewitz two-body results. Excludes 2H."""
    header_map = {"Helium 3": r"\HeT", "Hydrogen 3": r"\HThree", "Lithium 6": r"\LiS"}

    braun_nucs = ["Helium 3", "Hydrogen 3", "Lithium 6"]
    lenke_nucs = ["Helium 3", "Hydrogen 3"]

    braun_row_labels = [
        r"$F_T^a - F_T^b\ [\mathrm{fm}^{-1}]$",
        r"$E^{2N}_{0+}\ [10^{-3}/M_{\pi^+}]$",
    ]
    lenke_row_labels = [
        r"$F_T^a - F_T^b\ [\mathrm{fm}^{-1}]$",
        r"$F_L^a - F_L^b\ [\mathrm{fm}^{-1}]$",
        r"$E^{2N}_{0+}\ [10^{-3}/M_{\pi^+}]$",
        r"$L^{2N}_{0+}\ [10^{-3}/M_{\pi^+}]$",
    ]

    braun = {
        "Hydrogen 3": ("$-27.2(33)$", "$-3.55(43)$"),
        "Helium 3": ("$-27.1(33)$", "$-3.53(42)$"),
        "Lithium 6": ("$-11.4(14)$", "$-1.52(18)$"),
    }
    H3F_T = -29.7
    H3F_L = -23.2

    He3F_T = -29.3
    He3F_L = -22.9
    H32N = -4.01
    He32N = -3.95

    lenkewitz = {
        "Hydrogen 3": (f"${H3F_T}(3)$", f"${H3F_L}(1)$", f"${H32N}(3)$", "$-3.13$"),
        "Helium 3": (f"${He3F_T}(3)$", f"${He3F_L}(2)$", f"${He32N}(3)$", "$-3.09$"),
    }

    H3F_TO4 = H3F_T - 0.15
    H3F_LO4 = H3F_L - 0.49

    He3F_TO4 = He3F_T - 0.134
    He3F_LO4 = He3F_L - 0.542
    H32NO4 = H32N - 0.2
    He32NO4 = He32N - 0.2

    lenkewitz_o4 = {
        "Hydrogen 3": (f"${H3F_TO4:.2f}(3)$", f"${H3F_LO4:.2f}(1)$", f"${H32NO4:.2f}(3)$", "$-3.20$"),
        "Helium 3": (f"${He3F_TO4:.3f}(3)$", f"${He3F_LO4:.3f}(2)$", f"${He32NO4:.2f}(3)$", "$-3.15$"),
    }

    n_braun = len(braun_nucs)
    n_lenke = len(lenke_nucs)

    lines = []
    lines.append(r"\renewcommand{\arraystretch}{1.2}")
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")

    # --- Braun tabular ---
    braun_col_spec = "|l||" + "c|" * (n_braun - 1) + "c|"
    lines.append(r"\begin{tabular}{" + braun_col_spec + r"}\hline")

    # Header: title + order on same line
    n_braun_total = 1 + n_braun
    lines.append(
        r" & \multicolumn{" + str(n_braun)
        + r"}{c|}{Braun~\cite{BraunThesis} $\calO(\delta^2)$} \\"
    )
    lines.append(r"\cline{2-" + str(n_braun_total) + "}")

    # Nucleus names
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    nuc_headers = ["$" + header_map[nuc] + "$" for nuc in braun_nucs]
    lines.append(" & " + " & ".join(nuc_headers) + r" \\")
    lines.append(r"\hline")

    # Data rows
    for i, label in enumerate(braun_row_labels):
        cells = [label]
        for nuc in braun_nucs:
            cells.append(braun[nuc][i])
        lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(r"\vspace{1em}")

    # --- Lenkewitz tabular ---
    lenke_col_spec = "|l||" + "c|" * (n_lenke - 1) + "c||" + "c|" * (n_lenke - 1) + "c|"
    lines.append(r"\begin{tabular}{" + lenke_col_spec + r"}\hline")

    # Header: source name spanning all data columns
    n_lenke_total = 1 + 2 * n_lenke
    lines.append(
        r" & \multicolumn{" + str(2 * n_lenke)
        + r"}{c|}{Lenkewitz~\cite{LenkeThesis}} \\"
    )
    lines.append(r"\cline{2-" + str(n_lenke_total) + "}")

    # Order labels
    lines.append(
        r" & \multicolumn{" + str(n_lenke) + r"}{c||}{$\calO(\delta^2)$}"
        r" & \multicolumn{" + str(n_lenke) + r"}{c|}{$\calO(\delta^4)$} \\"
    )
    lines.append(r"\cline{2-" + str(n_lenke_total) + "}")

    # Nucleus names
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    nuc_headers = ["$" + header_map[nuc] + "$" for nuc in lenke_nucs] + [
        "$" + header_map[nuc] + "$" for nuc in lenke_nucs
    ]
    lines.append(" & " + " & ".join(nuc_headers) + r" \\")
    lines.append(r"\hline")

    # Data rows
    for i, label in enumerate(lenke_row_labels):
        cells = [label]
        for nuc in lenke_nucs:
            cells.append(lenkewitz[nuc][i])
        for nuc in lenke_nucs:
            cells.append(lenkewitz_o4[nuc][i])
        lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


def processTwoBodyFolder(folder, name):
    """Process all two-body files in folder and print results.

    Returns the analysis results dict.
    """
    all_files = _get_files_from_folder(folder)
    files = [f for f in all_files if "Odelta2" in f and "j12max=1" in f]

    print(f"\n{name} Two Body Results")
    labels = ["  F_T  ", "  F_L  ", "E_0+", "L_0+"]

    results = analyzeTwoBodyFiles(files)
    printTwoBodyTable(results, labels, precision=6)
    print(createLatexTable(results, name))
    return results


def runTwoBodyAnalysis():
    """Run analysis for all two-body calculations."""
    folders = {
        He3: r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/",
        H3: "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod/",
        He4: r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/",
        Li6: "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/",
    }

    print("2 Body Threshold Pion Photoproduction Results")
    all_results = []
    for name, folder in folders.items():
        results = processTwoBodyFolder(folder, name)
        if name != He4:
            all_results.append((name, results))
        print(breakLine)

    print("\n\nNow printing combined data")
    print(createCombinedTwoBodyLatexTable(all_results))
    print("\n\nNow printing literature Results\n\n")
    print(createLiteratureTwoBodyLatexTable())


if __name__ == "__main__":
    runTwoBodyAnalysis()
