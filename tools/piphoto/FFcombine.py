# -*- coding: utf-8 -*-
"""Combined form factor results for E1N, E2N, L1N, L2N across nuclei."""

import numpy as np

from FFSLib import He3, H3, Li6, _get_files_from_folder, _format_parens
from onebodyFFs import analyzeOneBodyFiles
from twobodyFFs import analyzeTwoBodyFiles

# Folder paths for each nucleus
onebody_folders = {
    H3: "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/1bod/thresh",
    He3: "/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/thresh/132MeV/",
    Li6: "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/1bod/thresh",
}

twobody_folders = {
    H3: "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod/",
    He3: "/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/",
    Li6: "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/",
}

# Nucleus order for the array (height 3)
nuc_order = [H3, He3, Li6]

# Comes from averaging
LiSfactor = 1.0


def getCombinedFFArray():
    """Get combined array of [E1N, E2N, E1N+E2N, a_0, L1N, L2N, L1N+L2N] for each nucleus.

    Returns:
        np.ndarray: Shape (3, 7) array where:
            - Rows are nuclei: [H3, He3, Li6]
            - Columns are: [E1N, E2N, E1N+E2N, a_0, L1N, L2N, L1N+L2N]
            where a_0 = (E1N + E2N)^2.
            Each value is the mean.
    """
    results = np.zeros((3, 7))

    for i, nuc in enumerate(nuc_order):
        # One-body: E1N is index 4, L1N is index 5 (mean is column 0)
        onebody_files = _get_files_from_folder(onebody_folders[nuc])
        onebody_data = analyzeOneBodyFiles(onebody_files)
        E1N = onebody_data[4, 0]
        L1N = onebody_data[5, 0]

        # Two-body: E2N is index 2, L2N is index 3 in o4 mean
        twobody_files = _get_files_from_folder(twobody_folders[nuc])
        twobody_files = [f for f in twobody_files if "Odelta2" in f and "j12max=1" in f]
        twobody_data = analyzeTwoBodyFiles(twobody_files)
        E2N = twobody_data["o4"]["mean"][2]
        L2N = twobody_data["o4"]["mean"][3]

        # a_0 = (E1N + E2N)^2
        a_0 = (E1N + E2N) ** 2
        if nuc == "6Li":
            a_0 = a_0 * LiSfactor
        results[i] = [E1N, E2N, E1N + E2N, a_0, L1N, L2N, L1N + L2N]

    return results


def getUncerArray():
    """Get uncertainties for [E1N, E2N, E1N+E2N, a_0, L1N, L2N, L1N+L2N] for each nucleus.

    Returns:
        np.ndarray: Shape (3, 7) array where:
            - Rows are nuclei: [H3, He3, Li6]
            - Columns are: [σ_E1N, σ_E2N, σ_(E1N+E2N), σ_a0, σ_L1N, σ_L2N, σ_(L1N+L2N)]
            Combined uncertainties use σ = sqrt(σ_1^2 + σ_2^2).
            For squared terms: σ_z = 2|x+y| sqrt(σ_x^2 + σ_y^2).
    """
    results = np.zeros((3, 7))

    for i, nuc in enumerate(nuc_order):
        # One-body: E1N is index 4, L1N is index 5 (spread is column 1)
        onebody_files = _get_files_from_folder(onebody_folders[nuc])
        onebody_data = analyzeOneBodyFiles(onebody_files)
        sigma_E1N = onebody_data[4, 1]
        sigma_L1N = onebody_data[5, 1]

        # Two-body: E2N is index 2, L2N is index 3 in o4 spread
        twobody_files = _get_files_from_folder(twobody_folders[nuc])
        twobody_files = [f for f in twobody_files if "Odelta2" in f and "j12max=1" in f]
        twobody_data = analyzeTwoBodyFiles(twobody_files)
        sigma_E2N = twobody_data["o4"]["spread"][2]
        sigma_L2N = twobody_data["o4"]["spread"][3]

        # Combined uncertainties: sqrt(σ_1^2 + σ_2^2)
        sigma_E_combined = np.sqrt(sigma_E1N**2 + sigma_E2N**2)
        sigma_L_combined = np.sqrt(sigma_L1N**2 + sigma_L2N**2)

        """
        Error propogation for z=f(x,y)=(x+y)^2 is
        σ_z= sqrt( (∂z/∂x)^2 σ_x^2 + (∂z/∂y)^2 σ_y^2)
           =2 sqrt[ (x+y)^2 (σ_x^2+σ_y^2) ]
        """
        # a_0 = (E1N + E2N)^2
        E1N = onebody_data[4, 0]
        E2N = twobody_data["o4"]["mean"][2]
        a_0 = (E1N + E2N) ** 2
        if nuc == "6Li":  # spin averaging
            a_0 = a_0 * LiSfactor
            sigma_E1N = sigma_E1N * LiSfactor
            sigma_E2N = sigma_E2N * LiSfactor

        arg_a = (a_0) * (sigma_E1N**2 + sigma_E2N**2)
        sigma_a0 = 2 * np.sqrt(arg_a)

        results[i] = [
            sigma_E1N,
            sigma_E2N,
            sigma_E_combined,
            sigma_a0,
            sigma_L1N,
            sigma_L2N,
            sigma_L_combined,
        ]

    return results


def createCombinedLatexTable(means, uncertainties):
    """Create LaTeX table for combined results.

    Args:
        means: Shape (3, 7) array of mean values.
        uncertainties: Shape (3, 7) array of uncertainties.

    Returns:
        str: LaTeX table string.
    """
    header_map = {H3: r"\HThree", He3: r"\HeT", Li6: r"\LiS"}

    col_labels = [
        r"$E^{1N}_{0+}$",
        r"$E^{2N}_{0+}$",
        r"$E_{0+}$",
        r"$a_0$",
        r"$L^{1N}_{0+}$",
        r"$L^{2N}_{0+}$",
        r"$L_{0+}$",
    ]

    lines = []
    lines.append(r"\renewcommand{\arraystretch}{1.4}")
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{|l||c|c|c|c||c|c|c|}\hline")

    # Header row
    lines.append(" & " + " & ".join(col_labels) + r" \\")
    lines.append(r"\hline")

    # Data rows
    for i, nuc in enumerate(nuc_order):
        nuc_label = "$" + header_map[nuc] + "$"
        cells = [nuc_label]
        for j in range(7):
            m, s = means[i, j], uncertainties[i, j]
            cells.append(_format_parens(m, s, precision=2))
        lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(
        r"\caption{Combined one-body and two-body results. Units are $[10^{-3}/m_\pi]$ for $E_{0+}$ and $L_{0+}$, and $[10^{-6}/m_\pi^2]$ for $a_0$.}"
    )
    lines.append(r"\end{table}")

    return "\n".join(lines)


# Array of combined results: shape (3, 7)
# Rows: [Hydrogen 3, Helium 3, Lithium 6]
# Columns: [E1N, E2N, E1N+E2N, a_0, L1N, L2N, L1N+L2N]
FFcombine = getCombinedFFArray()

# Array of uncertainties: shape (3, 7)
FFuncertainties = getUncerArray()


if __name__ == "__main__":
    print("Combined FF Results")
    print("Rows: Hydrogen 3, Helium 3, Lithium 6")
    print("Columns: E1N, E2N, E1N+E2N, a_0, L1N, L2N, L1N+L2N")
    print()
    print(FFcombine)
    print()
    print("Uncertainties")
    print()
    print(FFuncertainties)
    print()
    print("a_0 values [value, uncertainty] for each nucleus:")
    print()
    print("LaTeX Table:")
    print()
    print(createCombinedLatexTable(FFcombine, FFuncertainties))

    print(
        "Do not trust the values of scattering length, especially for Li6, there is spin averaging stuff that is going on"
    )
