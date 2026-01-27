# -*- coding: utf-8 -*-
"""Form factor and scattering amplitude utilities for pion photoproduction.

Extracts form factors (F_T, F_L) from nuclear matrix elements, and scattering
amplitudes (E_0+^1N, L_0+^1N) from output files of one-body and two-body
calculations (Odelta2 and Odelta4 orders).

One-body scattering amplitudes are read directly from the output text.
Two-body scattering amplitudes are computed from the Matrix section divided
by spin operators.

Author: alexl
"""

import numpy as np
from math import isinf, isnan
import readDensity as rd
from os import listdir
from os.path import isfile, join, dirname, basename
import array_to_latex as a2l
import re


# ============================================================================
# Core extraction functions
# ============================================================================

percent = 2  # assumed percent error from densities


def getSpinOperFromFile(path):
    """Extract spin operator matrices from output file.

    For one-body files: reads SpinVec from file.
    For two-body files: uses rd.spin4Nuc() as fallback.
    """
    with open(path, "r") as f:
        contents = f.read()

    nucName = rd.getNucName(path)

    expected = rd.spin4Nuc(nucName)

    # Try to read SpinVec (one-body files)
    if "onebody" in path:
        try:
            block = rd.getBlock(contents, "SpinVec")
            if block.strip():
                parsed_data = rd.parseBlock(block)
                if parsed_data:
                    result = rd.vals2matrix(nucName, parsed_data)
                    passed = np.allclose(result, expected)
                    if not passed:
                        print(f"SpinVec from file doesn't match spin4Nuc for {nucName}")
                        # print("result=\n", result)
                        # print("expected=\n", expected, "\n")
                        for i in range(3):
                            res = result[i]
                            exp = expected[i]
                            if not np.allclose(res, exp):
                                print("i=", i)
                                print(res, "\n")
                                print(exp, "\n")
                                print(res - exp)

                        print("\n\n")
                        # print("result/expected=", result / expected)

                    return result
        except AssertionError:
            raise  # checks for the spin vectors
        except Exception:
            pass

    # Fallback to rd.spin4Nuc() for two-body files
    return expected


def extractMatrixFromFile(path, kind="ScatMat"):
    """Extract matrix elements from output file."""
    matOut = rd.getQuantNums(path, kind=kind)
    return matOut["MatVals"]


def extractScalarFromFile(path, pattern, label):
    """Extract scalar value using regex pattern."""
    with open(path, "r") as f:
        match = re.search(pattern, f.read())

    if not match:
        raise ValueError(f"Could not find {label} in file {path}")

    return float(match.group(1))


def calculateFormFactors(mat, oper, divide=True):
    """Calculate form factors from matrix elements and spin operators.

    Returns [F_T, F_L] where F_T is average of x,y components and F_L is z.
    Multiplies by -1 to match literature conventions.
    """
    if divide:
        with np.errstate(divide="ignore", invalid="ignore"):
            fullMat = mat / oper
            fullMat = np.nan_to_num(fullMat, nan=0.0, posinf=0.0, neginf=0.0)
            tmpMat = fullMat.real
            # assert np.allclose(fullMat, tmpMat, atol=1e-5, rtol=1e-3), (
            #     f"fullMat\n={fullMat},\ntmpMat=\ntmpMat"
            # )
    else:
        tmpMat = mat.real

    formfacts = np.zeros(3)
    for i in range(3):
        flattened = tmpMat[i].flatten()
        valid = [f for f in flattened if not (isinf(f) or isnan(f)) and abs(f) > 0.0001]
        formfacts[i] = np.mean(valid) if valid else 0.0

    # print("formfacts\n=", formfacts)
    # Combine x,y into transverse and return [F_T, F_L]
    return np.array([(formfacts[0] + formfacts[1]) / 2, formfacts[2]]) * -1


def getFormFactorsFromFile(path, kind="FormFactors", divide=True):
    """Extract form factors [F_T, F_L] from output file."""
    mat = extractMatrixFromFile(path, kind=kind)
    oper = getSpinOperFromFile(path)
    return calculateFormFactors(mat, oper, divide=divide)


def getScatMatFromFile(path):
    """Extract E_0+^1N and L_0+^1N from output file.

    One-body: reads directly from text.
    Two-body: reads Matrix section and divides by spin operator.
    """
    if "onebody" in path:
        with open(path, "r") as f:
            contents = f.read()

        e_match = re.search(r"E_0\+\^1N=\s*([-\d.]+)", contents)
        l_match = re.search(r"L_0\+\^1N=\s*([-\d.]+)", contents)

        if not e_match or not l_match:
            raise ValueError(f"Could not find E_0+^1N or L_0+^1N in file {path}")

        return np.array([float(e_match.group(1)), float(l_match.group(1))])
    else:
        assert "twobody" in path
        mat = extractMatrixFromFile(path, kind="Matrix")
        oper = getSpinOperFromFile(path)

        return calculateFormFactors(mat, oper, divide=True)


def getCombinedResults(path):
    """Extract [F_T, F_L, E_0+, L_0+] from output file."""
    ffs = getFormFactorsFromFile(path)
    scatMat = getScatMatFromFile(path)
    return np.concatenate([ffs, scatMat])


# ============================================================================
# Two-body file pairing (Odelta2 vs Odelta4)
# ============================================================================


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
# Statistical analysis
# ============================================================================


def calculateSpread(values):
    """Calculate spread (max - min) / 2 for each column."""
    return (np.max(values, axis=0) - np.min(values, axis=0)) / 2


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


def analyzeOneBodyFiles(files):
    """Analyze one-body files for F^{S+V}, F^{S-V}, and ScatMat.

    Returns array of shape (6, 2) with [mean, spread] for each quantity.
    """
    kinds = ["F^{S+V}", "F^{S-V}", "ScatMat"]
    results = []

    for kind in kinds:
        # Extract values for this kind from all files
        if kind == "ScatMat":
            values = [getScatMatFromFile(path) for path in files]
        else:
            # Cancel out -1 factor in calculateFormFactors
            values = [
                getFormFactorsFromFile(path, kind=kind, divide=False) * -1
                for path in files
            ]

        values = np.array(values)
        mean = np.mean(values, axis=0)
        spread = calculateSpread(values)

        spread = np.sqrt(((percent * mean / 100) ** 2) + (spread**2))
        # Flatten to [[mean1, spread1], [mean2, spread2]]
        results.extend([[m, s] for m, s in zip(mean, spread)])

    return np.array(results)


# ============================================================================
# Output formatting
# ============================================================================


def _format_parens(mean, spread, precision=3):
    """Format mean(uncertainty) string for LaTeX, e.g. -0.047(5) or -0.047(13)."""
    scale = 10**precision
    unc = int(round(abs(spread) * scale))
    return rf"${mean:.{precision}f}({unc})$"


def _format_value(mean, spread, width, precision):
    """Format mean ± spread string centered in given width."""
    return f"{mean: .{precision}f} ± {spread:.{precision}f}".center(width)


def printOneBodyTable(data):
    """Print formatted table of one-body results."""
    labels = [
        "F_T^{S+V}",
        "F_L^{S+V}",
        "F_T^{S-V}",
        "F_L^{S-V}",
        "E_{0+}^{1N}  [10^-3/m_π]",
        "L_{0+}^{1N}  [10^-3/m_π]",
    ]

    col_width = 25
    prec = 3
    label_width = max(len(label) for label in labels)
    total_width = label_width + col_width + 10

    print("-" * total_width)
    print(f"{'One Body Results':<{label_width}} | {' (mean ± σ)':^{col_width}} | ")
    print("-" * total_width)

    for i, label in enumerate(labels):
        if i in [2, 4]:
            print("-" * total_width)
        print(
            f"{label:<{label_width}} | {_format_value(data[i, 0], data[i, 1], col_width, prec)} | "
        )

    print("-" * total_width)
    print()


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


def createOneBodyLatexTable(all_data):
    """Create LaTeX table for one-body results across all nuclei.

    all_data: dict mapping nucleus name to (6, 2) array from analyzeOneBodyFiles.
    Rows are nuclei, columns are quantities.
    """
    header_map = {
        "Lithium 6": r"\LiS",
        "Helium 4": r"\HeF",
        "Helium 3": r"\HeT",
        "Hydrogen 3": r"\HThree",
    }

    nuc_order = ["Hydrogen 3", "Helium 3", "Lithium 6"]

    col_labels = [
        r"$F_T^{S+V}$",
        r"$F_T^{S-V}$",
        r"$E^{1N}_{0+}\ [10^{-3}/m_\pi]$",
        r"$F_L^{S+V}$",
        r"$F_L^{S-V}$",
        r"$L^{1N}_{0+}\ [10^{-3}/m_\pi]$",
    ]
    # Indices into the (6,2) data array for the column order above
    col_data_idx = [0, 2, 4, 1, 3, 5]

    num_rows = 1 + len(nuc_order)  # 1 header + nuclei
    num_cols = 1 + len(col_labels)  # 1 label + quantities

    out = np.empty((num_rows, num_cols), dtype="U50")
    out[:, :] = ""

    # Header row: quantity labels
    for i, label in enumerate(col_labels):
        out[0, i + 1] = label

    # Data rows: one per nucleus
    for j, nuc in enumerate(nuc_order):
        out[j + 1, 0] = r"\hline$" + header_map.get(nuc, nuc) + "$"
        if nuc in all_data:
            data = all_data[nuc]
            for i, di in enumerate(col_data_idx):
                m, s = data[di, 0], data[di, 1]
                out[j + 1, i + 1] = _format_parens(m, s)

    col_spec = "l|" + "c|" * (len(col_labels) - 1) + "c"
    latex_str = a2l.to_ltx(out, arraytype="tabular", print_out=False)
    replaceStr = (
        r"\begin{table}[H]"
        + "\n"
        + r"\centering"
        + "\n"
        + r"\begin{tabular}{"
        + col_spec
        + "}"
    )
    latex_str = latex_str.replace(r"\begin{tabular}", replaceStr)
    latex_str = r"\renewcommand{\arraystretch}{1.4}" + "\n" + latex_str
    latex_str += "\n" + r"\end{table}"
    return latex_str


def createComparisonLatexTable():
    """Create LaTeX comparison table of F_T results: Braun vs Lenkewitz.

    Excludes 2H. Shows F_T^{S+V}, F_T^{S-V}, E_{0+}^{1N} with hardcoded
    values from Braun and Lenkewitz.
    """
    # Rows: [nucleus_label, Braun (F_T^{S+V}, F_T^{S-V}, E_{0+}^{1N}),
    #                       Lenkewitz (F_T^{S+V}, F_T^{S-V}, E_{0+}^{1N})]
    data = [
        # ³H
        [
            r"$\HThree$",
            "$1.551(78)$",
            "$0.039(2)$",
            "$-0.94(5)(5)$",
            "$1.493(25)$",
            "$0.012(13)$",
            "$-0.93(3)(5)$",
        ],
        # ³He
        [
            r"$\HeT$",
            "$0.041(2)$",
            "$1.544(77)$",
            "$1.77(9)(9)$",
            "$0.017(13)$",
            "$1.480(26)$",
            "$1.71(4)(9)$",
        ],
        # ⁶Li
        [r"$\LiS$", "$0.476(24)$", "$0.479(24)$", "$0.26(3)(3)$", "", "", ""],
    ]

    lines = []
    lines.append(r"\renewcommand{\arraystretch}{1.2}")
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{|l||c|c|c||c|c|c|}\hline")
    lines.append(
        r" & \multicolumn{3}{c||}{Braun}"
        r" & \multicolumn{3}{c|}{Lenkewitz} \\"
    )
    lines.append(r"\hline")
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    lines.append(
        r" & $F_T^{S+V}$ & $F_T^{S-V}$"
        r" & $E^{1N}_{0+}\ [10^{-3}/M_{\pi^+}]$"
        r" & $F_T^{S+V}$ & $F_T^{S-V}$"
        r" & $E^{1N}_{0+}\ [10^{-3}/M_{\pi^+}]$ \\"
    )
    lines.append(r"\hline")

    for row in data:
        lines.append(" & ".join(row) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


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
# High-level interface
# ============================================================================


def _get_files_from_folder(folder):
    """Get all files from a folder."""
    return [join(folder, f) for f in listdir(folder) if isfile(join(folder, f))]


def processOneBodyFolder(folder, name):
    """Process all one-body files in folder and print results."""
    files = _get_files_from_folder(folder)
    print(f"\n{name} One Body Results")
    data = analyzeOneBodyFiles(files)
    printOneBodyTable(data)


def processTwoBodyFolder(folder, name):
    """Process all two-body files in folder and print results."""
    all_files = _get_files_from_folder(folder)
    files = [f for f in all_files if "Odelta2" in f and "j12max=1" in f]

    print(f"\n{name} Two Body Results")
    labels = ["  F_T  ", "  F_L  ", "E_0+", "L_0+"]

    results = analyzeTwoBodyFiles(files)
    printTwoBodyTable(results, labels, precision=6)
    print(createLatexTable(results, name))


# ============================================================================
# Main execution
# ============================================================================

# Nucleus names
He3 = "Helium 3"
H3 = "Hydrogen 3"
He4 = "Helium 4"
Li6 = "Lithium 6"
breakLine = "\n" + 150 * "%"
breakLine2 = "\n" + 100 * "%"


def runTwoBodyAnalysis():
    """Run analysis for all two-body calculations."""
    folders = {
        He3: r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/",
        H3: "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod/",
        He4: r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/",
        Li6: "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/",
    }

    print("2 Body Threshold Pion Photoproduction Results")
    for name, folder in folders.items():
        processTwoBodyFolder(folder, name)
        print(breakLine)


def runOneBodyAnalysis():
    """Run analysis for all one-body calculations."""
    folders = {
        He3: r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/thresh/132MeV/",
        H3: "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/1bod/thresh",
        He4: "/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/1bod/thresh",
        Li6: "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/1bod/thresh",
    }

    all_data = {}
    for name, folder in folders.items():
        processOneBodyFolder(folder, name)
        files = _get_files_from_folder(folder)
        all_data[name] = analyzeOneBodyFiles(files)
        print(breakLine2)

    print(createOneBodyLatexTable(all_data))
    print()
    print(createComparisonLatexTable())


if __name__ == "__main__":
    runTwoBodyAnalysis()
    print("\n\n")
    print(150 * "#")
    print(150 * "#")
    print("\n\n")
    runOneBodyAnalysis()
