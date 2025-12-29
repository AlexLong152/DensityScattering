import re
import numpy as np
from typing import Dict
from copy import copy


def main():
    testfile = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/132MeV/onebody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.Odelta3-.denshash=fac31e8f17bc2b1a6d84411b408e05cca71aae63c04be8f6fbf3ba311acb45a3.v2.0.dat"
    d = parseOnebody(testfile)

    for k, v in d.items():
        print(f"{k}: {v}")
    d = He3Parse(testfile)
    print(10 * "\n")
    for k, v in d.items():
        print(f"{k}: {v}")


def He3Parse(filepath):
    myDict = parseOnebody(filepath)
    one = myDict[1]
    two = myDict[2]
    Transverse = np.vstack((one, two))
    vals = np.mean(Transverse, axis=0)
    three = myDict[3]
    del myDict[2]
    myDict[1] = vals
    tmp = copy(myDict[1])
    myDict[1] = np.array([f.real for f in copy(myDict[1])])
    myDict[3] = np.array([f.real for f in copy(myDict[3])])
    return myDict


def parseOnebody(filepath: str) -> Dict[object, object]:
    """
    Given a path to an output-*.dat file, return a dictionary with:
      - metadata parsed from the densityFileName substring:
          myDict["lambdaCut"], myDict["lambdaSRG"], myDict["Ntotmax"],
          myDict["omegaH"], myDict["denshash"]
      - per-eps arrays stored as myDict[1], myDict[2], ... (dtype=complex)
        Each array contains the values encountered in that eps block, in order.
        Length can be from 0 to 8.

    for 3He
    myDict[1] contains averages eps=1,0,0 and eps=0,1,0 results
    myDict[1] contains eps=0,0,1 results
    """
    myDict: Dict[object, object] = {}

    # Regex for the densityFileName line
    re_dens_line = re.compile(r"^\s*densityFileName=(?P<path>.+?)\s*$")

    # Extract fields and denshash
    re_meta = re.compile(
        r"lambda(?P<cut>\d+)-lambdaSRG(?P<srg>[0-9.]+)setNtotmax(?P<nmax>\d+)omegaH(?P<omega>\d+)\.denshash=(?P<hash>[^.\s]+)"
    )

    # eps header like: "eps=1,0,0 Result"
    re_eps = re.compile(
        r"^\s*eps\s*=\s*([-0-9]+)\s*,\s*([-0-9]+)\s*,\s*([-0-9]+)\s*Result", re.I
    )

    # F lines like: "FMinus(-1,1) =   1.5013 +0.0000j"
    re_value = re.compile(r"=\s*([+-]?\d+(?:\.\d+)?)")

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # ---- Parse metadata ----
    for line in lines:
        m = re_dens_line.match(line)
        if m:
            path = m.group("path")
            mm = re_meta.search(path)
            if mm:
                myDict["lambdaCut"] = int(mm.group("cut"))
                myDict["lambdaSRG"] = float(mm.group("srg"))
                myDict["Ntotmax"] = int(mm.group("nmax"))
                myDict["omegaH"] = int(mm.group("omega"))
                myDict["denshash"] = mm.group("hash")
            break

    # ---- Parse eps blocks ----
    idx = 1
    current_values = []

    def commit_block():
        nonlocal idx, current_values
        if current_values:
            myDict[idx] = np.unique(np.array(current_values, dtype="complex"))
            idx += 1
            current_values = []

    for line in lines:
        if line.strip().startswith("Lenkewitz"):
            commit_block()
            break

        if re_eps.match(line):
            commit_block()  # close previous block if any
            current_values = []  # start new block
            continue

        if "FPlus" in line or "FMinus" in line:
            mv = re_value.search(line)
            if mv:
                val = float(mv.group(1))
                current_values.append(complex(val, 0.0))

    commit_block()  # commit last block if file ends cleanly

    return myDict


if __name__ == "__main__":
    main()
