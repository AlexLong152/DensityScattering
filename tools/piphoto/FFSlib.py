# -*- coding: utf-8 -*-

"""
@author: alexl

v2: Standardized uncertainty slots.
  For all quantities: unc1 = nuclear structure / form factor, unc2 = elementary amplitude.
  Form factors: [mean, theory_spread, numeric_unc] (no elementary amplitude dependence).
  E1N/L1N: [mean, sqrt(theory^2 + numeric^2), multipole_unc].
  E2N/L2N: [mean, sqrt(theory^2 + numeric^2), 0] (no elementary amplitude dependence).
"""

import readDensity as rd
import numpy as np
from os import listdir
from os.path import isfile, join
from math import isinf, isnan

# 3% numeric uncertainty for all TDA results
err = 3 / 100

# Elementary amplitudes in 10^{-3}/m_{pi+}
E0p_pi0p = -1.16
E0p_pi0n = +2.13
L0p_pi0p = -1.35
L0p_pi0n = -2.41
amp_frac = 0.05  # 5% uncertainty on elementary amplitudes

# Masses in MeV for K1N calculation
m_N = 938.918
m_pi = 139.571
M_nuc = {"3H": 2808.921, "3He": 2808.391, "6Li": 5601.518}


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


def K1N(nuc):
    """Kinematical factor for 1N to 3N phase space (Eq. 6, Lenkewitz 2013)."""
    Mnuc = M_nuc[nuc]
    return (m_N + m_pi) / (Mnuc + m_pi) * Mnuc / m_N


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
    else:
        tmpMat = mat.real

    formfacts = np.zeros(3)
    for i in range(3):
        flattened = tmpMat[i].flatten()
        valid = [f for f in flattened if not (isinf(f) or isnan(f)) and abs(f) > 0.0001]
        formfacts[i] = np.mean(valid) if valid else 0.0

    # Combine x,y into transverse and return [F_T, F_L]
    return np.array([(formfacts[0] + formfacts[1]) / 2, formfacts[2]]) * -1


def extractMatrixFromFile(path, kind="ScatMat"):
    """Extract matrix elements from output file."""
    matOut = rd.getQuantNums(path, kind=kind)
    return matOut["MatVals"]


def _analyze_onebody_files(files):
    """Extract F_T and F_L form factors from output files.

    Returns dict: {quantity_name: np.array([mean, theory_unc, numeric_unc])}
    where theory_unc is half the spread and numeric_unc is 1.5% of |mean|.
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
        theory_unc = (np.max(arr) - np.min(arr)) / 2
        numeric_unc = err * abs(mean)
        result[k] = np.array([mean, theory_unc, numeric_unc])

    return result


def _analyze_above_files(files):
    """Extract E1N and L1N from above-threshold output files.

    Returns dict: {"E1N": np.array([mean, theory_unc, numeric_unc]),
                    "L1N": np.array([mean, theory_unc, numeric_unc])}
    where theory_unc is half the spread and numeric_unc is 1.5% of |mean|.
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
        theory_unc = (np.max(arr) - np.min(arr)) / 2
        numeric_unc = err * abs(mean)
        result[k] = np.array([mean, theory_unc, numeric_unc])

    return result


def _calc_multipole_unc(nuc, onebody_data):
    """Compute theory uncertainty on E1N/L1N from 5% elementary amplitude uncertainty.

    Uses E1N = K/2 * (Ep * F_T^{S+V} + En * F_T^{S-V}) (Eq. 5, Lenkewitz 2013)
    and varies only the elementary amplitudes, holding form factors fixed.

    Returns dict: {"E1N": unc_value, "L1N": unc_value}
    """
    k = K1N(nuc)

    ft_spv = onebody_data["F_T^S+V"][0]
    ft_smv = onebody_data["F_T^S-V"][0]
    fl_spv = onebody_data["F_L^S+V"][0]
    fl_smv = onebody_data["F_L^S-V"][0]

    ep, en = E0p_pi0p, E0p_pi0n
    dep, den = abs(ep) * amp_frac, abs(en) * amp_frac
    lp, ln = L0p_pi0p, L0p_pi0n
    dlp, dln = abs(lp) * amp_frac, abs(ln) * amp_frac

    # E1N: vary elementary E amplitudes over all sign combinations
    e_corners = []
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            val = k / 2 * ((ep + s1 * dep) * ft_spv + (en + s2 * den) * ft_smv)
            e_corners.append(val)
    e_unc = (max(e_corners) - min(e_corners)) / 2

    # L1N: vary elementary L amplitudes over all sign combinations
    l_corners = []
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            val = k / 2 * ((lp + s1 * dlp) * fl_spv + (ln + s2 * dln) * fl_smv)
            l_corners.append(val)
    l_unc = (max(l_corners) - min(l_corners)) / 2

    return {"E1N": e_unc, "L1N": l_unc}


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

    Returns dict with display-name keys:
      Form factors: np.array([mean, theory_unc, numeric_unc])
      E2N/L2N: np.array([mean, combined_nuclear_unc, 0])
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
            theory_unc = (np.max(col) - np.min(col)) / 2
            numeric_unc = err * abs(mean)
            if k in ("E2N", "L2N"):
                # Combined nuclear structure uncertainty; no elementary amplitude dependence
                result[f"{name[k]}{tag}"] = np.array(
                    [mean, np.sqrt(theory_unc**2 + numeric_unc**2), 0]
                )
            else:
                result[f"{name[k]}{tag}"] = np.array([mean, theory_unc, numeric_unc])

    # Static = Total - base
    for k in keys:
        mean = result[f"{name[k]} O(q^4)"][0] - result[f"{name[k]} O(q^3)"][0]
        numeric_unc = err * abs(mean)
        if k in ("E2N", "L2N"):
            result[f"Static {name[k]}"] = np.array([mean, numeric_unc, 0])
        else:
            result[f"Static {name[k]}"] = np.array([mean, 0, numeric_unc])

    return result


def get_literature_results():
    """Return nested dictionary of results.

    Structure: d[author][nucleus][quantity] = np.array([mean, unc1, unc2])
    For TDA results:
      Form factors: unc1 = potential variation, unc2 = numerical precision (1.5%)
      E1N/L1N: unc1 = total nuclear structure, unc2 = elementary amplitude uncertainty
      E2N/L2N: unc1 = combined nuclear structure, unc2 = 0
    Literature entries are kept as originally published.
    Entries without results are NoResult instances.
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
    d["braun"]["3H"]["F_T^S+V"] = np.array([1.551, 0.078, 0])
    d["braun"]["3H"]["F_T^S-V"] = np.array([0.039, 0.002, 0])
    d["braun"]["3H"]["E1N"] = np.array([-0.94, 0.05, 0.05])

    # 3He
    d["braun"]["3He"]["F_T^S+V"] = np.array([0.041, 0.002, 0])
    d["braun"]["3He"]["F_T^S-V"] = np.array([1.544, 0.077, 0])
    d["braun"]["3He"]["E1N"] = np.array([1.77, 0.09, 0.09])

    # 6Li
    d["braun"]["6Li"]["F_T^S+V"] = np.array([0.476, 0.024, 0])
    d["braun"]["6Li"]["F_T^S-V"] = np.array([0.479, 0.024, 0])
    d["braun"]["6Li"]["E1N"] = np.array([0.26, 0.03, 0.03])

    # === Lenkewitz results ===
    # F_T from Braun Table 5.2 "literature" column [Ref 26]
    d["lenke"]["3He"]["F_T^S+V"] = np.array([0.017, 0.013, 0])
    d["lenke"]["3He"]["F_T^S-V"] = np.array([1.480, 0.026, 0])
    d["lenke"]["3H"]["F_T^S+V"] = np.array([1.493, 0.025, 0])
    d["lenke"]["3H"]["F_T^S-V"] = np.array([0.012, 0.013, 0])

    # F_L from Lenkewitz 2013 Table 1
    d["lenke"]["3He"]["F_L^S+V"] = np.array([-0.079, 0.014, 0.008])
    d["lenke"]["3He"]["F_L^S-V"] = np.array([1.479, 0.026, 0.008])
    d["lenke"]["3H"]["F_L^S+V"] = np.array([1.487, 0.027, 0.008])
    d["lenke"]["3H"]["F_L^S-V"] = np.array([-0.083, 0.014, 0.008])

    # E1N from Braun Table 5.2 "literature" (= Lenkewitz Table 4 first column)
    d["lenke"]["3He"]["E1N"] = np.array([1.71, 0.04, 0.09])
    d["lenke"]["3H"]["E1N"] = np.array([-0.93, 0.03, 0.05])

    # L1N from Lenkewitz 2013 Table 4 first column (1N at O(q^4))
    d["lenke"]["3He"]["L1N"] = np.array([-1.89, 0.04, 0.09])
    d["lenke"]["3H"]["L1N"] = np.array([-0.99, 0.04, 0.05])

    # 6Li: no Lenkewitz results (all remain NoResult)
    # Braun: no F_L or L1N results (all remain NoResult)

    # ============================================================
    # Two-body results
    # ============================================================

    # === Braun Table 5.4 "our result" ===
    # 3H
    d["braun"]["3H"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-27.2, 3.3, 0])
    d["braun"]["3H"]["E2N O(q^3)"] = np.array([-3.55, 0.43, 0])

    # 3He
    d["braun"]["3He"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-27.1, 3.3, 0])
    d["braun"]["3He"]["E2N O(q^3)"] = np.array([-3.53, 0.42, 0])

    # 6Li
    d["braun"]["6Li"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-11.4, 1.4, 0])
    d["braun"]["6Li"]["E2N O(q^3)"] = np.array([-1.52, 0.18, 0])

    # === Lenkewitz two-body results ===
    # F_T from Braun Table 5.4 "literature" column [Ref 26]
    d["lenke"]["3H"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-29.7, 0.3, 0])
    d["lenke"]["3He"]["F^{(a)}_T-F^{(b)}_T O(q^3)"] = np.array([-29.3, 0.3, 0])

    # F_L from Lenkewitz 2013 Table 2
    d["lenke"]["3He"]["F^{(a)}_L-F^{(b)}_L O(q^3)"] = np.array([-22.9, 0.2, 0.1])
    d["lenke"]["3H"]["F^{(a)}_L-F^{(b)}_L O(q^3)"] = np.array([-23.2, 0.1, 0.1])

    # E2N from Braun Table 5.4 "literature" column [Ref 26]
    d["lenke"]["3H"]["E2N O(q^3)"] = np.array([-4.01, 0.03, 0])
    d["lenke"]["3He"]["E2N O(q^3)"] = np.array([-3.95, 0.03, 0])

    # L2N from Lenkewitz 2013 Table 4 second column (2N at O(q^3))
    d["lenke"]["3He"]["L2N O(q^3)"] = np.array([-3.09, 0.02, 0])
    d["lenke"]["3H"]["L2N O(q^3)"] = np.array([-3.13, 0.01, 0])

    # === Lenkewitz 2N-static form factors from Thesis Table 5.8 ===
    d["lenke"]["3He"]["Static F^{(a)}_T-F^{(b)}_T"] = np.array([-0.134, 0.006, 0.072])
    d["lenke"]["3H"]["Static F^{(a)}_T-F^{(b)}_T"] = np.array([-0.150, 0.050, 0.084])
    d["lenke"]["3He"]["Static F^{(a)}_L-F^{(b)}_L"] = np.array([-0.542, 0.056, 0.110])
    d["lenke"]["3H"]["Static F^{(a)}_L-F^{(b)}_L"] = np.array([-0.490, 0.038, 0.099])

    # O(q^4) total = O(q^3) + Static, for form factors
    for nuc in ["3H", "3He"]:
        for ff in ["F^{(a)}_T-F^{(b)}_T", "F^{(a)}_L-F^{(b)}_L"]:
            base = d["lenke"][nuc][f"{ff} O(q^3)"]
            static = d["lenke"][nuc][f"Static {ff}"]
            mean = base[0] + static[0]
            theory = np.sqrt(base[1] ** 2 + static[1] ** 2)
            numeric = np.sqrt(base[2] ** 2 + static[2] ** 2)
            d["lenke"][nuc][f"{ff} O(q^4)"] = np.array([mean, theory, numeric])

    # === Lenkewitz 2N-static (q^4) from Table 4 column 5 ===
    d["lenke"]["3He"]["Static E2N"] = np.array([-0.02, 0.00, 0.01])
    d["lenke"]["3He"]["Static L2N"] = np.array([-0.07, 0.01, 0.01])
    d["lenke"]["3H"]["Static E2N"] = np.array([-0.02, 0.01, 0.01])
    d["lenke"]["3H"]["Static L2N"] = np.array([-0.07, 0.00, 0.01])

    # O(q^4) total = O(q^3) + Static, for E2N/L2N
    for nuc in ["3H", "3He"]:
        for prefix in ["E2N", "L2N"]:
            base = d["lenke"][nuc][f"{prefix} O(q^3)"]
            static = d["lenke"][nuc][f"Static {prefix}"]
            mean = base[0] + static[0]
            theory = np.sqrt(base[1] ** 2 + static[1] ** 2)
            numeric = np.sqrt(base[2] ** 2 + static[2] ** 2)
            d["lenke"][nuc][f"{prefix} O(q^4)"] = np.array([mean, theory, numeric])

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

        # Combine nuclear structure uncertainties; keep elementary amplitude separate
        multipole_unc = _calc_multipole_unc(nuc, onebody_data)
        for amp in ["E1N", "L1N"]:
            nuclear_total = np.sqrt(above_data[amp][1] ** 2 + above_data[amp][2] ** 2)
            d["TDA"][nuc][amp] = np.array(
                [above_data[amp][0], nuclear_total, multipole_unc[amp]]
            )

        # Two-body from 2bod/ folders
        twobody_files = _get_files_from_folder(twobody_folders[nuc])
        twobody_data = _analyze_twobody_files(twobody_files)
        for q in twobody_data:
            d["TDA"][nuc][q] = twobody_data[q]

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
                        for i in range(3):
                            res = result[i]
                            exp = expected[i]
                            if not np.allclose(res, exp):
                                print("i=", i)
                                print(res, "\n")
                                print(exp, "\n")
                                print(res - exp)

                        print("\n\n")

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
