# -*- coding: utf-8 -*-

"""
@author: alexl
"""

# import os
# import sys
import readDensity as rd
import numpy as np
from os import listdir
from os.path import isfile, join
from math import isinf, isnan


class NoResult:
    pass


quantities_1bod = ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V", "E1N", "L1N"]
quantities_2bod = [
    "F^{(a)}_T-F^{(b)}_T O(q^3)",
    "F^{(a)}_L-F^{(b)}_L O(q^3)",
    "E2N O(q^3)",
    "L2N O(q^3)",
    "Static F^{(a)}_T-F^{(b)}_T",
    "Static F^{(a)}_L-F^{(b)}_L",
    "Static E2N",
    "Static L2N",
    "F^{(a)}_T-F^{(b)}_T O(q^4)",
    "F^{(a)}_L-F^{(b)}_L O(q^4)",
    "E2N O(q^4)",
    "L2N O(q^4)",
]
quantities = quantities_1bod + quantities_2bod
nuclei = ["3H", "3He", "6Li"]

# ============================================================
# Folder paths for TDA (this work)
# ============================================================

onebody_folders = {
    "3H": "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/1bod/thresh",
    "3He": "/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/thresh/132MeV/",
    "6Li": "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/1bod/thresh",
}

above_folders = {
    "3H": "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/1bod/above",
    "3He": "/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/above",
    "6Li": "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/1bod/above",
}

twobody_folders = {
    "3H": "/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod/",
    "3He": "/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/",
    "6Li": "/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/",
}


def getFormFactorsFromFile(path, kind="FormFactors", divide=True):
    """Extract form factors [F_T, F_L] from output file."""
    mat = extractMatrixFromFile(path, kind=kind)
    oper = getSpinOperFromFile(path)
    return calculateFormFactors(mat, oper, divide=divide)


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


def extractMatrixFromFile(path, kind="ScatMat"):
    """Extract matrix elements from output file."""
    matOut = rd.getQuantNums(path, kind=kind)
    return matOut["MatVals"]


def _analyze_onebody_files(files):
    """Extract F_T and F_L form factors from output files, return mean and spread.

    For each file, getFormFactorsFromFile (with divide=False) returns [F_T, F_L]
    using the average of non-zero values from extQnum=1,2 for F_T
    and extQnum=3 for F_L.  The *-1 cancels the sign convention inside FFSLib.

    Returns dict: {quantity_name: np.array([mean, spread])}
    """
    ff_keys = ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V"]
    all_vals = {k: [] for k in ff_keys}

    for path in files:
        for kind, sv in [("F^{S+V}", "S+V"), ("F^{S-V}", "S-V")]:
            ft, fl = getFormFactorsFromFile(path, kind=kind, divide=False) * -1
            all_vals[f"F_T^{sv}"].append(ft)
            all_vals[f"F_L^{sv}"].append(fl)

    result = {}
    for k in ff_keys:
        arr = np.array(all_vals[k])
        mean = np.mean(arr)
        spread = (np.max(arr) - np.min(arr)) / 2
        result[k] = np.array([mean, spread])

    return result


def _analyze_above_files(files):
    """Extract E1N and L1N from above-threshold output files, return mean and spread.

    For each file, extracts the FormFacts matrix and uses calculateFormFactors
    with divide=False. The *-1 cancels the sign convention inside calculateFormFactors.
    Returns [E1N, L1N] where E1N is from extQnum=1,2 and L1N is from extQnum=3.

    Returns dict: {"E1N": np.array([mean, spread]), "L1N": np.array([mean, spread])}
    """
    all_E1N = []
    all_L1N = []

    for path in files:
        mat = extractMatrixFromFile(path, kind="FormFacts")
        e1n, l1n = calculateFormFactors(mat, None, divide=False) * -1
        all_E1N.append(e1n)
        all_L1N.append(l1n)

    result = {}
    for k, vals in [("E1N", all_E1N), ("L1N", all_L1N)]:
        arr = np.array(vals)
        mean = np.mean(arr)
        spread = (np.max(arr) - np.min(arr)) / 2
        result[k] = np.array([mean, spread])

    return result


def _get_combined_twobody(path):
    """Extract [F_T, F_L, E2N, L2N] from a twobody output file.

    Uses divide=True (divides by spin operator), matching twobodyFFs.py.
    """
    ffs = getFormFactorsFromFile(path, kind="FormFactors", divide=True)
    scatmat = getFormFactorsFromFile(path, kind="Matrix", divide=True)
    return np.concatenate([ffs, scatmat])


def _find_odelta4_match(path_o2):
    """Given an Odelta2 file path, find the matching Odelta4 file."""
    from os.path import basename, dirname

    folder = dirname(path_o2)
    name_o4 = basename(path_o2).replace("Odelta2", "Odelta4")
    path_o4 = join(folder, name_o4)
    if isfile(path_o4):
        return path_o4
    raise FileNotFoundError(f"No Odelta4 match for {path_o2}")


def _analyze_twobody_files(files):
    """Extract F_T, F_L, E2N, L2N from twobody output files.

    Filters for Odelta2 files with j12max=1, pairs each with its Odelta4 match.
    Uses divide=True (divides by spin operator) matching twobodyFFs.py.

    Returns dict with display-name keys, each np.array([mean, spread]).
    """
    o2_files = [f for f in files if "Odelta2" in f and "j12max=1" in f]

    o2_all = []
    o4_all = []
    for path_o2 in o2_files:
        path_o4 = _find_odelta4_match(path_o2)
        o2_all.append(_get_combined_twobody(path_o2))
        o4_all.append(_get_combined_twobody(path_o4))

    o2_arr = np.array(o2_all)  # shape (n_files, 4): [F_T, F_L, E2N, L2N]
    o4_arr = np.array(o4_all)

    # Internal keys and display name mapping
    keys = ["F_T", "F_L", "E2N", "L2N"]
    name = {
        "F_T": "F^{(a)}_T-F^{(b)}_T",
        "F_L": "F^{(a)}_L-F^{(b)}_L",
        "E2N": "E2N",
        "L2N": "L2N",
    }

    result = {}
    for i, k in enumerate(keys):
        for tag, arr in [(" O(q^3)", o2_arr), (" O(q^4)", o4_arr)]:
            col = arr[:, i]
            mean = np.mean(col)
            spread = (np.max(col) - np.min(col)) / 2
            result[f"{name[k]}{tag}"] = np.array([mean, spread])

    # Static = Total - base (uncertainty 0)
    for k in keys:
        mean = result[f"{name[k]} O(q^4)"][0] - result[f"{name[k]} O(q^3)"][0]
        result[f"Static {name[k]}"] = np.array([mean, 0])

    return result


def get_literature_results():
    """Return nested dictionary of results from Braun (Table 5.2), Lenkewitz (2013), and TDA.

    Structure: d[author][nucleus][quantity] = np.array([mean, combined_uncertainty])
    Entries without results are NoResult instances.

    Combined uncertainty: sqrt(uncertainty_1^2 + uncertainty_2^2) when two are given.
    """
    d = {}
    for author in ["braun", "lenke", "TDA"]:
        d[author] = {}
        for nuc in nuclei:
            d[author][nuc] = {}
            for q in quantities:
                d[author][nuc][q] = NoResult()

    # === Braun Table 5.2 "our result" ===
    # 3H
    d["braun"]["3H"]["F_T^S+V"] = np.array([1.551, 0.078])
    d["braun"]["3H"]["F_T^S-V"] = np.array([0.039, 0.002])
    d["braun"]["3H"]["E1N"] = np.array([-0.94, np.sqrt(0.05**2 + 0.05**2)])

    # 3He
    d["braun"]["3He"]["F_T^S+V"] = np.array([0.041, 0.002])
    d["braun"]["3He"]["F_T^S-V"] = np.array([1.544, 0.077])
    d["braun"]["3He"]["E1N"] = np.array([1.77, np.sqrt(0.09**2 + 0.09**2)])

    # 6Li
    d["braun"]["6Li"]["F_T^S+V"] = np.array([0.476, 0.024])
    d["braun"]["6Li"]["F_T^S-V"] = np.array([0.479, 0.024])
    d["braun"]["6Li"]["E1N"] = np.array([0.26, np.sqrt(0.03**2 + 0.03**2)])

    # === Lenkewitz results ===
    # F_T from Braun Table 5.2 "literature" column [Ref 26]
    d["lenke"]["3He"]["F_T^S+V"] = np.array([0.017, 0.013])
    d["lenke"]["3He"]["F_T^S-V"] = np.array([1.480, 0.026])
    d["lenke"]["3H"]["F_T^S+V"] = np.array([1.493, 0.025])
    d["lenke"]["3H"]["F_T^S-V"] = np.array([0.012, 0.013])

    # F_L from Lenkewitz 2013 Table 1
    d["lenke"]["3He"]["F_L^S+V"] = np.array([-0.079, np.sqrt(0.014**2 + 0.008**2)])
    d["lenke"]["3He"]["F_L^S-V"] = np.array([1.479, np.sqrt(0.026**2 + 0.008**2)])
    d["lenke"]["3H"]["F_L^S+V"] = np.array([1.487, np.sqrt(0.027**2 + 0.008**2)])
    d["lenke"]["3H"]["F_L^S-V"] = np.array([-0.083, np.sqrt(0.014**2 + 0.008**2)])

    # E1N from Braun Table 5.2 "literature" (= Lenkewitz Table 4 first column)
    d["lenke"]["3He"]["E1N"] = np.array([1.71, np.sqrt(0.04**2 + 0.09**2)])
    d["lenke"]["3H"]["E1N"] = np.array([-0.93, np.sqrt(0.03**2 + 0.05**2)])

    # L1N from Lenkewitz 2013 Table 4 first column (1N at O(q^4))
    d["lenke"]["3He"]["L1N"] = np.array([-1.89, np.sqrt(0.04**2 + 0.09**2)])
    d["lenke"]["3H"]["L1N"] = np.array([-0.99, np.sqrt(0.04**2 + 0.05**2)])

    # 6Li: no Lenkewitz results (all remain NoResult)
    # Braun: no F_L or L1N results (all remain NoResult)

    # ============================================================
    # Two-body results
    # ============================================================

    # === Braun Table 5.4 "our result" ===
    # 3H
    d["braun"]["3H"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-27.2, 3.3])
    d["braun"]["3H"]["E2N O(q^3)"] = np.array([-3.55, 0.43])

    # 3He
    d["braun"]["3He"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-27.1, 3.3])
    d["braun"]["3He"]["E2N O(q^3)"] = np.array([-3.53, 0.42])

    # 6Li
    d["braun"]["6Li"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-11.4, 1.4])
    d["braun"]["6Li"]["E2N O(q^3)"] = np.array([-1.52, 0.18])

    # === Lenkewitz two-body results ===
    # F_T from Braun Table 5.4 "literature" column [Ref 26]
    d["lenke"]["3H"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-29.7, 0.3])
    d["lenke"]["3He"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-29.3, 0.3])

    # F_L from Lenkewitz 2013 Table 2
    d["lenke"]["3He"]["F^{(a)}_L-F^{(b)}_L O(q^3)"] = np.array(
        [-22.9, np.sqrt(0.2**2 + 0.1**2)]
    )
    d["lenke"]["3H"]["F^{(a)}_L-F^{(b)}_L O(q^3)"] = np.array(
        [-23.2, np.sqrt(0.1**2 + 0.1**2)]
    )

    # E2N from Braun Table 5.4 "literature" column [Ref 26]
    d["lenke"]["3H"]["E2N O(q^3)"] = np.array([-4.01, 0.03])
    d["lenke"]["3He"]["E2N O(q^3)"] = np.array([-3.95, 0.03])

    # L2N from Lenkewitz 2013 Table 4 second column (2N at O(q^3))
    d["lenke"]["3He"]["L2N O(q^3)"] = np.array([-3.09, 0.02])
    d["lenke"]["3H"]["L2N O(q^3)"] = np.array([-3.13, 0.01])

    # === Lenkewitz 2N-static form factors from Thesis Table 5.8 ===
    d["lenke"]["3He"]["Static F^{(a)}_T-F^{(b)}_T"] = np.array(
        [-0.134, np.sqrt(0.006**2 + 0.072**2)]
    )
    d["lenke"]["3H"]["Static F^{(a)}_T-F^{(b)}_T"] = np.array(
        [-0.150, np.sqrt(0.050**2 + 0.084**2)]
    )
    d["lenke"]["3He"]["Static F^{(a)}_L-F^{(b)}_L"] = np.array(
        [-0.542, np.sqrt(0.056**2 + 0.110**2)]
    )
    d["lenke"]["3H"]["Static F^{(a)}_L-F^{(b)}_L"] = np.array(
        [-0.490, np.sqrt(0.038**2 + 0.099**2)]
    )

    # O(q^4) total = O(q^3) + Static, for form factors
    for nuc in ["3H", "3He"]:
        for ff in ["F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L"]:
            base = d["lenke"][nuc][f"{ff} O(q^3)"]
            static = d["lenke"][nuc][f"Static {ff}"]
            mean = base[0] + static[0]
            sigma = np.array([base[1], static[1]])
            sigma = np.sqrt(np.dot(sigma, sigma))
            d["lenke"][nuc][f"{ff} O(q^4)"] = np.array([mean, sigma])

    # === Lenkewitz 2N-static (q^4) from Table 4 column 5 ===
    d["lenke"]["3He"]["Static E2N"] = np.array([-0.02, np.sqrt(0.00**2 + 0.01**2)])
    d["lenke"]["3He"]["Static L2N"] = np.array([-0.07, np.sqrt(0.01**2 + 0.01**2)])
    d["lenke"]["3H"]["Static E2N"] = np.array([-0.02, np.sqrt(0.01**2 + 0.01**2)])
    d["lenke"]["3H"]["Static L2N"] = np.array([-0.07, np.sqrt(0.00**2 + 0.01**2)])

    # O(q^4) total = O(q^3) + Static, for E2N/L2N
    for nuc in ["3H", "3He"]:
        for prefix in ["E2N", "L2N"]:
            base = d["lenke"][nuc][f"{prefix} O(q^3)"]
            static = d["lenke"][nuc][f"Static {prefix}"]
            mean = base[0] + static[0]
            sigma = np.array([base[1], static[1]])
            sigma = np.sqrt(np.dot(sigma, sigma))
            d["lenke"][nuc][f"{prefix} O(q^4)"] = np.array([mean, sigma])

    # 6Li: no Lenkewitz two-body results (all remain NoResult)
    # Braun: no 2bodF_L or L2N results (all remain NoResult)

    # ============================================================
    # TDA results (this work) — one-body form factors from files
    # ============================================================

    for nuc in nuclei:
        # Form factors from thresh/ folders
        files = _get_files_from_folder(onebody_folders[nuc])
        onebody_data = _analyze_onebody_files(files)
        for q in ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V"]:
            d["TDA"][nuc][q] = onebody_data[q]

        # E1N, L1N from above/ folders
        above_files = _get_files_from_folder(above_folders[nuc])
        above_data = _analyze_above_files(above_files)
        d["TDA"][nuc]["E1N"] = above_data["E1N"]
        d["TDA"][nuc]["L1N"] = above_data["L1N"]

        # Two-body from 2bod/ folders
        twobody_files = _get_files_from_folder(twobody_folders[nuc])
        twobody_data = _analyze_twobody_files(twobody_files)
        for q in twobody_data:
            d["TDA"][nuc][q] = twobody_data[q]

    # Add 3% uncertainty to all TDA results
    for nuc in nuclei:
        for q in d["TDA"][nuc]:
            val = d["TDA"][nuc][q]
            if not isinstance(val, NoResult):
                old_unc = val[1]
                mean = val[0]
                val[1] = np.sqrt(old_unc**2 + (0.04 * mean) ** 2)

    return d


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


def _get_files_from_folder(folder):
    """Get all files from a folder."""
    return [join(folder, f) for f in listdir(folder) if isfile(join(folder, f))]
