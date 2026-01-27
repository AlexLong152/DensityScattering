# -*- coding: utf-8 -*-
"""One-body form factor analysis for pion photoproduction."""

import numpy as np
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

breakLine2 = "\n" + 100 * "%"


# ============================================================================
# Statistical analysis (one-body)
# ============================================================================


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
# Output formatting (one-body)
# ============================================================================


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


def createBraunLatexTable():
    """Create LaTeX table of Braun one-body results. Excludes 2H."""
    data = [
        [r"$\HThree$", "$1.551(78)$", "$0.039(2)$", "$-0.94(5)(5)$"],
        [r"$\HeT$", "$0.041(2)$", "$1.544(77)$", "$1.77(9)(9)$"],
        [r"$\LiS$", "$0.476(24)$", "$0.479(24)$", "$0.26(3)(3)$"],
    ]

    lines = []
    lines.append(r"\renewcommand{\arraystretch}{1.2}")
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{|l||c|c|c|}\hline")
    lines.append(r" & \multicolumn{3}{c|}{Braun~\cite{BraunThesis}} \\")
    lines.append(r"\cline{2-4}")
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    lines.append(
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


def createLenkewitzLatexTable():
    """Create LaTeX table of Lenkewitz one-body results. Excludes 2H."""
    data = [
        [
            r"$\HThree$",
            "$1.493(25)$",
            "$0.012(13)$",
            "$-0.93(3)(5)$",
            "$1.487(27)(8)$",
            "$-0.083(14)(8)$",
        ],
        [
            r"$\HeT$",
            "$0.017(13)$",
            "$1.480(26)$",
            "$1.71(4)(9)$",
            "$-0.079(14)(8)$",
            "$1.479(26)(8)$",
        ],
    ]

    lines = []
    lines.append(r"\renewcommand{\arraystretch}{1.2}")
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{|l||c|c|c|c|c|}\hline")
    lines.append(r" & \multicolumn{5}{c|}{Lenkewitz~\cite{LenkeThesis}} \\")
    lines.append(r"\cline{2-6}")
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    lines.append(
        r" & $F_T^{S+V}$ & $F_T^{S-V}$"
        r" & $E^{1N}_{0+}\ [10^{-3}/M_{\pi^+}]$"
        r" & $F_L^{S+V}$ & $F_L^{S-V}$ \\"
    )
    lines.append(r"\hline")

    for row in data:
        lines.append(" & ".join(row) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(r"\caption{")
    lines.append(
        r"  Numerical form factors from the literature. The first error is estimation of theory uncertainty,"
    )
    lines.append(
        r"  and the second error (where present) comes from the uncertainty in single nucleon amplitudes."
    )
    lines.append(r"}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ============================================================================
# High-level interface (one-body)
# ============================================================================


def processOneBodyFolder(folder, name):
    """Process all one-body files in folder and print results."""
    files = _get_files_from_folder(folder)
    print(f"\n{name} One Body Results")
    data = analyzeOneBodyFiles(files)
    printOneBodyTable(data)


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
    print(createBraunLatexTable())
    print()
    print(createLenkewitzLatexTable())


if __name__ == "__main__":
    runOneBodyAnalysis()
