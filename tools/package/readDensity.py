# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import re
from copy import copy
from defs import pauliVec, spin1Vec, getSpinOper


def main():
    file1 = r"/home/alexander/OneDrive/DensityScattering/tools/package/twobody-3He.132MeV-060deg.dens-chiralsmsN4LO+3nfN2LO-lambda450-lambdaSRG0.000setNtotmax00omegaH00.Odelta2-j12max=1-.denshash=e8f1ba5e7fafac8431e05ae0e6b63381f668664f52d39abfcea0830f420f4706.v2.0.dat"
    for f in [file1]:
        tmp = getQuantNums(f, kind="FormFactors")["MatVals"]
        print(tmp)
        print(75 * "%")
        tmp = getQuantNums(f, kind="Matrix")["MatVals"]
        print(tmp)
        print(75 * "%")


def getBlock(filetext, kind):
    """
    Extracts the text block corresponding to 'kind' from the file contents.

    Parameters
    ----------
    filetext : str
        The full text contents of the file
    kind : str
        The block type to extract (e.g., "FormFactors" or "Matrix")

    Returns
    -------
    str
        The text block containing lines matching the kind pattern
    """
    lines = filetext.splitlines()
    block_lines = []

    for line in lines:
        # Check if line contains the kind pattern
        if kind + "(extQnum=" in line:
            block_lines.append(line)

    if not block_lines:
        raise BadDataError(f"No block found for kind='{kind}'")

    return "\n".join(block_lines)


def parseBlock(block):
    """
    Parses a block of text and extracts quantum numbers and complex values.

    Parameters
    ----------
    block : str
        Text block containing lines with quantum numbers and values

    Returns
    -------
    list
        List of tuples, each containing (extQnum, twoMzp, twoMz, complex_number)
    """

    lines = block.splitlines()
    result = []

    for line in lines:
        # Parse pattern: kind(extQnum=X,twoMzp=Y, twoMz=Z): real_part +/- imag_part i
        # Example: FormFactors(extQnum=   1,twoMzp=   1, twoMz=   1):   -0.053... +0.000... i

        # Extract quantum numbers
        match = re.search(
            r"extQnum=\s*(\d+),\s*twoMzp=\s*(-?\d+),\s*twoMz=\s*(-?\d+)", line
        )
        if not match:
            BadDataError(f"Bad data, line={line}")

        extQnum = int(match.group(1))
        twoMzp = int(match.group(2))
        twoMz = int(match.group(3))

        # Extract complex number (after the colon)
        value_part = line.split(":", 1)[1].strip()

        # Replace 'i' with 'j' and remove spaces, then cast to complex
        value_part = value_part.replace(" ", "").replace("i", "j")
        complex_number = complex(value_part)

        result.append((extQnum, twoMzp, twoMz, complex_number))

    return result


def getQuantNums(filename, returnMat=True, kind="Result"):
    """
    Reads in quantum numbers as real numbers from filename
    """
    out = {}
    nucName = getNucName(filename)
    out["nuc"] = nucName
    if returnMat:
        with open(filename, "r") as f:
            contents = f.read()
            if isinstance(kind, (list, np.ndarray)):
                for tmp in kind:
                    try:
                        block = getBlock(contents, kind)
                    except BadDataError:
                        pass
            else:
                block = getBlock(contents, kind)
            parsed_data = parseBlock(block)

            out["MatVals"] = vals2matrix(nucName, parsed_data)
    try:
        densityName = densityFileName(filename)
    except UnboundLocalError:
        densityName = filename

    out["twoSpin"] = getTwoSpin(nucName)
    densityName = densityFileName(filename)
    out["name"] = filename
    out["file"] = filename
    out["hash"] = getHash(densityName)

    out["omega"] = getOmega(densityName)
    out["energy"] = getOmega(densityName)

    out["theta"] = getTheta(densityName)
    out["angle"] = getTheta(densityName)

    out["Ntotmax"] = getNtotmax(densityName)
    out["omegaH"] = getOmegaH(densityName)
    out["lambda"] = getLambda(densityName)
    out["lambdaCut"] = getLambda(densityName)
    # print("filename=", filename)
    out["lambdaSRG"] = getLambdaSRG(densityName)
    out["numBodies"] = getNumBodies(densityName)
    out["Odelta"] = getOdelta(densityName)
    # print("len(vals)=", len(vals))
    # print(np.shape(out["MatVals"]))
    return out


def getTwoSpin(NucName):
    out = {"3He": 1, "3H": 1, "4He": 0, "6Li": 2}
    return out[NucName]


# def data2array(parsed_data):
class BadDataError(Exception):
    def __init__(self, message):
        super().__init__(message)


def densityFileName(filename):
    """
    Gets the density file name from input file
    if the input file exists at the end of the output file
    """
    substr = r"**************************************************"  # gets contents in input file
    with open(filename, "r") as f:
        contents = f.read()
        contents = contents.split(substr)[-1]
        contents = contents.split("\n")
        # magic number comes from input file density being on 6th line of input file
        # for i, line in enumerate(contents):
        #     print(i, ":", line)
        found = False
        for line in contents:
            if "denshash" in line:
                inFile = line
                found = True
                break
        if found:
            inFile = inFile.split(r"/")[-1]
        else:
            inFile = filename
        # print("boop")
        # print("inFile=", inFile)
    return inFile


def spin4Nuc(nuc):
    """
     Given the nucleus `nuc` as a string, returns the spin vector operator.

       - 6Li  : spin 1   → 3 spin-1 matrices (3x3)
       - 3He  : spin 1/2 → Pauli matrices / 2 (2x2)
       - 4He  : spin 0   →


     Parameters
     ----------
     nuc : str
         The nucleus name as a string. Supported values are:
         - "6Li": spin 1, returns the three 3×3 spin-1 matrices.
         - "3He": spin 1/2, returns Pauli matrices divided by 2.
         - "4He": spin 0, Spin operator isn't a thing return (1,1,1)

    Returns
     -------
     numpy.ndarray
         A length-3 array whose elements are the spin operator matrices
         appropriate for the nuclear spin. For spin-0 (4He) the operator
         is returned as the vector [1.0, 1.0, 1.0].
    Raises
    ------
    ValueError
        If the nucleus name is not recognized.
    """

    match nuc:
        case "4He":
            return getSpinOper(0)

        case "3He":
            return getSpinOper(1)

        case "3H":
            return getSpinOper(1)
        case "6Li":
            return getSpinOper(2)

        case _:
            raise ValueError(f"Unknown nucleus '{nuc}'")


def vals2matrix(nucName, parsed_data):
    """
    Converts parsed data to a 3D matrix using proper quantum number to index mapping.

    Parameters
    ----------
    nucName : str
        Nucleus name ("3He", "4He", or "6Li")
    parsed_data : list
        List of tuples (extQnum, twoMzp, twoMz, complex_number)

    Returns
    -------
    numpy.ndarray
        3D array with shape (numextQnums, numStates, numStates)

    Notes
    -----
    The mapping from quantum numbers to array indices follows the convention
    from IsospinMap.py and the Fortran code: index = (twoSpin - twoMz) // 2
    This ensures consistency with the spin matrix definitions.
    """
    # Determine twoSpin from nucleus name
    match nucName:
        case "6Li":
            twoSpin = 2
        case "4He":
            twoSpin = 0
        case "3He":
            twoSpin = 1
        case "3H":
            twoSpin = 1
        case _:
            raise ValueError(f"Unknown nucleus: {nucName}")

    numStates = twoSpin + 1  # 2S+1 quantum states

    # Find number of external quantum numbers
    extQnums = set(item[0] for item in parsed_data)
    numextQnums = len(extQnums)

    # Initialize output array
    out = np.zeros((numextQnums, numStates, numStates), dtype=np.complex128)

    # Mapping function: twoMz -> array index
    # Follows the convention from IsospinMap.py: index = (twoSpin - twoMz) // 2
    def get_index(twoMz):
        return (twoSpin - twoMz) // 2

    # Fill the array using quantum numbers from parsed_data
    allZeros = True
    for extQnum, twoMzp, twoMz, complex_number in parsed_data:
        # extQnum starts at 1, array index starts at 0
        ext_idx = extQnum - 1

        # Map quantum numbers to array indices
        mzp_idx = get_index(twoMzp)
        mz_idx = get_index(twoMz)

        # Store the value
        out[ext_idx, mzp_idx, mz_idx] = complex_number

        if complex_number != 0:
            allZeros = False

    if allZeros:
        raise ValueError("All zero entries in parsed data")

    return out


def getStringbtwn(myString, prefix, suffix):
    assert prefix in myString
    assert suffix in myString
    myString = copy(myString)
    myString = " " + myString
    myString = myString.split(prefix)[1]
    myString = myString.split(suffix)[0]
    return myString


def getNucName(filepath):
    path = copy(filepath)
    path = path.split(r"/")[-1]
    name = getStringbtwn(path, "body", "MeV")

    nuclei = ["6Li", "4He", "3He", "3H"]

    for nuc in nuclei:
        if nuc in name:
            if "3H" in name:
                assert "3He" not in filepath[:20]
            return nuc

    print("readDensity.py:getNucName path=", path)
    print("name=", name)
    raise ValueError("Nucleus not found")


def getOdelta(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    btwn = getStringbtwn(filename, "omegaH", "-.")
    return btwn[3:]


def getOmega(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    omega = filename.split("MeV")[0][-3:]
    omega = int(omega)
    return omega


def getNumBodies(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    if "twobody" in filename:
        return 1
    elif "onebody" in filename:
        return 2
    else:
        raise ValueError(f"Num bodies not found for {filename}")


def getTheta(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    theta = filename.split("deg")[0][-3:]
    theta = int(theta)
    return theta


def getNtotmax(filename):
    """
    'onebody-6Li.060MeV-180deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax014omegaH24.denshash=b7a1bdaca9e2d2ec9bedc928a210baa4fc442b27027a63fbb2512686170f40f9.v2.0.h5'
    """
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    Ntotmax = getStringbtwn(filename, "setNtotmax", "omegaH")
    if Ntotmax == "":
        Ntotmax = 0
    else:
        Ntotmax = int(Ntotmax)
    return Ntotmax


def getOmegaH(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    # print("filename=", filename)
    # omegaH = filename.split("omegaH")[-1][:2]
    omegaH = getStringbtwn(filename, "omegaH", r".")
    if omegaH == "":
        omegaH = 0
    else:
        omegaH = int(omegaH)
    return omegaH


def getLambdaSRG(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    lambdaSRG = filename.split("lambdaSRG")[1]
    lambdaSRG = lambdaSRG[:5]
    return float(lambdaSRG)


def getLambda(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    lambdaVal = filename.split("lambda")[1]
    lambdaVal = lambdaVal[:-1]
    return int(lambdaVal)


def getHash(filename):
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
    if "denshash" in filename:
        hash = filename.split("denshash=")[-1]
        hash = hash.split(".v2")[0]
        return hash
    else:
        return ""


if __name__ == "__main__":
    main()
