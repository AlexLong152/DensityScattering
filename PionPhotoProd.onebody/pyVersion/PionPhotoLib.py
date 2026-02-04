# -*- coding: utf-8 -*-

"""
@author: alexl

TODO: implement beyond l=1 in sum using experimental data
Note for extrapolation: E_l, M_l behave like |vec{q}|^l in the threshold region
"""

# import sys
import os
from pdb import set_trace
import numpy as np
from scipy.integrate import quad
import os
from os import listdir
from os.path import isfile, join
import re

# from copy import copy


def passFoo():
    """
    :return: None
    Trivial function for turning off trace when desired
    """
    pass


DEBUG = True
trace = set_trace if DEBUG else passFoo


sigx = np.array([[0, 1], [1, 0]])
sigy = np.array([[0, -1j], [1j, 0]])
sigz = np.array([[1, 0], [0, -1]])
sigVec = np.array([sigx, sigy, sigz])

# TODO: change out mpi for mpi plus
MeVtofm = 197.3  # 1fm =(1/197.3) MeV^-1
# TODO: mpi masses appear to matter a lot, so add sig figs, and check for mpi+ vs mpi0
mpi = 134.97
mpiPlus = 139.57
mProton = 938.272
mNeutron = 939.565
mN = (mProton + mNeutron) / 2

"""
Based on N. Rijneveen thesis. Inverting equations in 3.12 and solving for 3.11 gives

T_gamma p -> pi^0 p= T_gamma p ^(1/2) + (2/3) T_gamma N^(3/2)
T_gamma n -> pi^0 n= -T_gamma n ^(1/2) + (2/3) T_gamma N^(3/2)
"""


def main():
    """
    Run a basic example calculation.
    """
    poleData = SaidPoles()
    x = 0.5  # cos(theta)
    nucs = "pp0"  # reaction type
    sqrtS = 1200  # Mandelstam S
    # print(calcCrossSection(S, x, nucs, poleData))
    # getRawM(sqrtS, x, nucs, poleData)

    print(calcCrossSection(sqrtS**2, x, nucs, poleData))
    # print(calcCrossSectionFromRaw(sqrtS**2, x, nucs, poleData))


def calcCCFromRaw(S, x, nucs, data, roundThresh=False):
    """
    calculates the crossSection from `getRawM`

    Parameters
    ----------
    roundThresh: bool, optional
        If True, round Epi to threshold when within 2.5 MeV (matching Fortran behavior).
        Default is False.
    """
    epsVecs = np.array([[1, 0, 0], [0, 1, 0]])
    sqrtS = np.sqrt(S)
    S, kVec, qVec = getKinematics(sqrtS, x, nucs, roundThresh=roundThresh)
    # prefactor = 8 * np.pi * sqrtS
    MSquare = 0
    for epsVec in epsVecs:
        MMat = getRawM(sqrtS, x, nucs, data, epsVec=epsVec, roundThresh=roundThresh)
        M_Dot_MDag = MMat @ MMat.conjugate().T
        MSquare += np.trace(M_Dot_MDag).real

    crossSec = (vecAbs(qVec) / vecAbs(kVec)) * MSquare
    crossSec = 0.25 * crossSec * (1 / (64 * S * np.pi**2))
    # crossSec = crossSec * S * np.pi
    crossSec = crossSec * MeVtofm**2
    crossSec = crossSec / 100  # Convert to barns
    crossSec = crossSec / (10**-6)  # Convert to microbarns
    return crossSec


def calcCrossSection(S, x, nucs, data, roundThresh=False):
    """
    Calculate the cross section for a given reaction.
    Result is in microbarns

    Parameters
    -----------
    S: float
        mandalstam S
    x: float
        cos(theta) where theta is the scattering angle
    nucs: str
        options are: "pp0" "nn0" "pn+" "np-"
        specifies the reaction,
        string order is nucleus before, nucleus after, and outgoing pion charge
        from getKinematics
        ```
        if nucs == "pp0":
            mNucl = mProton
            mPion = mpi
        elif nucs == "nn0":
            mNucl = mNeutron
            mPion = mpi
        elif nucs == "pn+":
            mNucl = mProton
            mPion = mpiPlus
        elif nucs == "np-":
            mNucl = mNeutron
            mPion = mpiPlus
        ```
    data: SaidPoles
        Object containing pole data for the calculation
    roundThresh: bool, optional
        If True, round Epi to threshold when within 2.5 MeV (matching Fortran behavior).
        Default is False.

    Returns
    -------
    float
        The cross section in microbarns
    """
    sqrtS = np.sqrt(S)
    S, kVec, qVec = getKinematics(sqrtS, x, nucs, roundThresh=roundThresh)
    MSquare = getMSquare(sqrtS, x, nucs, data, roundThresh=roundThresh)
    crossSec = (vecAbs(qVec) / vecAbs(kVec)) * MSquare
    crossSec = 0.25 * crossSec * (1 / (64 * S * np.pi**2))
    # crossSec = crossSec * S * np.pi
    crossSec = crossSec * MeVtofm**2
    crossSec = crossSec / 100  # Convert to barns
    crossSec = crossSec / (10**-6)  # Convert to microbarns
    return crossSec


def getRawM(sqrtS, x, nucs, data, epsVec=None, roundThresh=False):
    """
    Calculate the matrix (not matrix squared), unitless

    Parameters
    ----------
    sqrtS: float
        The square root of the Mandelstam variable S
    x: float
        cos(theta) where theta is the scattering angle
    nucs: str
        Specifies the reaction ("pp0", "nn0", "pn+", "np-")
    data: SaidPoles
        Object containing pole data for the calculation
    epsVec: ndarray, optional, length 3
        The polarization vector of the incoming photon
        if None, assigns epsVec=[1,0,0]
    roundThresh: bool, optional
        If True, round Epi to threshold when within 2.5 MeV (matching Fortran behavior).
        Default is False.

    Returns
    -------
    Mat: ndarray 2x2 complex
        The matrix element
    """

    _, kVec, qVec = getKinematics(sqrtS, x, nucs, roundThresh=roundThresh)
    targets = ["p12", "n12", "32q"]

    if nucs == "pp0":
        coefs = np.array([1, 0, 2 / 3])
    elif nucs == "nn0":
        coefs = np.array([0, -1, 2 / 3])
    elif nucs == "pn+":
        coefs = np.array([1, 0, -1 / 3]) * np.sqrt(2)
    elif nucs == "np-":
        coefs = np.array([0, 1, 1 / 3]) * np.sqrt(2)
    else:
        raise ValueError(f"nucs value '{nucs}' not supported")

    # Initialize the total squared amplitude
    # Prefactor for the amplitude
    if type(epsVec) is type(None):
        epsVec = np.array([1, 0, 0])  # WLOG set to x polarization if not set

    assert len(epsVec) == 3
    assert np.dot(epsVec, epsVec) == 1

    # Remember to square then sum not sum then square
    prefactor = 8 * np.pi * sqrtS
    # MSquare = np.zeros((2, 2), dtype=complex)
    Mmat = np.zeros((2, 2), dtype=complex)
    for i, target in enumerate(targets):
        Mmat += coefs[i] * F(x, sqrtS, qVec, kVec, epsVec, data, target)
    Mmat *= prefactor

    return Mmat


def getMSquare(sqrtS, x, nucs, data, roundThresh=False):
    """
    Calculate the squared matrix element.

    Parameters
    ----------
    sqrtS: float
        The square root of the Mandelstam variable S
    x: float
        cos(theta) where theta is the scattering angle
    nucs: str
        Specifies the reaction ("pp0", "nn0", "pn+", "np-")
    data: SaidPoles
        Object containing pole data for the calculation
    roundThresh: bool, optional
        If True, round Epi to threshold when within 2.5 MeV (matching Fortran behavior).
        Default is False.

    Returns
    -------
    float
        The squared matrix element
    """

    _, kVec, qVec = getKinematics(sqrtS, x, nucs, roundThresh=roundThresh)

    # Use 3 linear polarizations (x, y, z) to match Fortran implementation
    epsVecs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=complex)
    epsVecs = np.array([[1, 0, 0], [0, 1, 0]], dtype=complex)
    targets = ["p12", "n12", "32q"]

    if nucs == "pp0":
        coefs = np.array([1, 0, 2 / 3])
    elif nucs == "nn0":
        coefs = np.array([0, -1, 2 / 3])
    elif nucs == "pn+":
        coefs = np.array([1, 0, -1 / 3]) * np.sqrt(2)
    elif nucs == "np-":
        coefs = np.array([0, 1, 1 / 3]) * np.sqrt(2)
    else:
        raise ValueError(f"nucs value '{nucs}' not supported")

    # Initialize the total squared amplitude
    # Prefactor for the amplitude

    # Remember to square then sum not sum then square
    prefactor = 8 * np.pi * sqrtS
    MSquare = np.zeros((2, 2), dtype=complex)
    sumMSq = 0.0
    for epsVec in epsVecs:
        # build the 2×2 amplitude for this polarization only
        Mpol = np.zeros((2, 2), dtype=complex)
        for i, target in enumerate(targets):
            Mpol += coefs[i] * F(x, sqrtS, qVec, kVec, epsVec, data, target)
        Mpol *= prefactor

        # now square _that_ one, sum over final spins, then add to sumMSq
        Mpol_ms = Mpol @ Mpol.conjugate().T
        trace = np.trace(Mpol_ms).real
        sumMSq += trace
    return sumMSq
    # MSquare = np.sum(Mtmp_MS).real  # both are equivalent, trace is slightly more accurate with floating point
    # MSquare now contains the sum over the two photon polarizations of |M|^2.


def printMat(mat, val=None):
    assert np.shape(mat) == (2, 2)
    if val is not None:
        for i in range(len(mat[0])):
            for j in range(len(mat[:, 0])):
                print(f"{val}[{i},{j}]=", f"{mat[i, j]:.6E}")
    else:
        for i in range(len(mat[0])):
            for j in range(len(mat[:, 0])):
                print(f"Mat[{i},{j}]=", f"{mat[i, j]:.6E}")

    print("")


def getKinematics(sqrtS, x, nucs, roundThresh=False):
    """
    Calculate kinematic variables for a given reaction.

    Parameters
    ----------
    sqrtS: float
        The square root of the Mandelstam variable S
    x: float
        cos(theta) where theta is the scattering angle
    nucs: str
        Specifies the reaction ("pp0", "nn0", "pn+", "np-")
    roundThresh: bool, optional
        If True, round Epi to threshold when within 2.5 MeV (matching Fortran behavior).
        Default is False.

    Returns
    -------
    S, kVec, qVec - Mandelstam S, photon momentum vector, pion momentum vector
    """
    if nucs == "pp0":
        mNucl = mProton
        mPion = mpi
    elif nucs == "nn0":
        mNucl = mNeutron
        mPion = mpi
    elif nucs == "pn+":
        mNucl = mProton
        mPion = mpiPlus
    elif nucs == "np-":
        mNucl = mNeutron
        mPion = mpiPlus
    else:
        print()
        errStr = "Given value of `nucs` has length " + str(len(nucs))
        errStr += f"\n nucs value '{nucs}' not supported"
        raise ValueError(errStr)
    S = sqrtS**2
    omega = (S - mNucl**2) / (2 * sqrtS)  # cm frame
    kVec = np.array([0, 0, omega])

    # stu = 2 * mN**2 + mpi**2
    Epi = (S + mPion**2 - mNucl**2) / (2 * sqrtS)

    # Match Fortran behavior: set Epi to threshold if within 2.5 MeV
    if roundThresh and abs(Epi - mPion) <= 2.5:
        Epi = mPion

    # Enucl = (S - mpi**2 + mN**2) / (2 * sqrtS)
    # print("Epi=", Epi)
    # print("Enucl=", Enucl)
    absQ = np.sqrt(Epi**2 - mPion**2)
    qVec = np.array([0, np.sqrt(1 - x**2), x]) * absQ

    # abspp = np.sqrt(Enucl**2 - mN**2)
    # qVec = np.array([0, np.sin(x), np.cos(x)]) * absQ

    return S, kVec, qVec


def getFullKinematics(sqrtS, x):
    """
    Calculate full set of kinematic variables.

    Parameters
    ----------
    sqrtS: float
        The square root of the Mandelstam variable S
    x: float
        cos(theta) where theta is the scattering angle

    Returns
    -------
    tuple
        Multiple kinematic variables including S, momentum vectors, and energies
    """
    # theta = np.arccos(x)
    S = sqrtS**2
    omega = (S - mN**2) / (2 * sqrtS)
    kVec = np.array([0, 0, omega])

    # stu = 2 * mN**2 + mpi**2

    Epi = (S + mpi**2 - mN**2) / (2 * sqrtS)
    Enucl = (S - mpi**2 + mN**2) / (2 * sqrtS)
    # print("Epi=", Epi)
    # print("Enucl=", Enucl)
    absQ = np.sqrt(Epi**2 - mpi**2)
    abspp = np.sqrt(Enucl**2 - mN**2)
    qVec = np.array([0, np.sqrt(1 - x**2), x]) * absQ
    # qVec = np.array([0, np.sin(x), np.cos(x)]) * absQ

    omegaLab = (S - mN**2) / (2 * mN)  # lab frame
    # print("omegaLab=", omegaLab)

    # In the center of mass frame with kVec=k zHat, p=-k zHat
    pVec = -1 * kVec  # cm frame
    k0 = omega
    p0 = np.sqrt(mN**2 + omega**2)

    ppVec = -1 * qVec
    # ppVec = np.array([0, np.sin(theta), np.cos(theta)]) * -1 * abspp

    # k4 = np.array([omega, 0, 0, omega])
    # p4 = np.array([np.sqrt(omega**2 + mN**2), 0, 0, omega])
    q0 = np.sqrt(mpi**2 + np.dot(qVec, qVec))
    pp0 = np.sqrt(mN**2 + np.dot(ppVec, ppVec))
    assert abs(Epi == q0) < 0.0001
    assert abs(Enucl - pp0) < 0.0001
    assert abs(abspp - vecAbs(ppVec)) < 0.0001
    # print("Assertions passed: Epi==q0, Enucl==pp0, abspp==|ppVec|\n")
    return S, kVec, qVec, pVec, pVec, ppVec, k0, p0, q0, pp0, omega


def dag(M):
    """
    Calculate the Hermitian conjugate (dagger) of a matrix.

    Parameters
    ----------
    M: numpy.ndarray
        Input matrix

    Returns
    -------
    numpy.ndarray
        Hermitian conjugate of the input matrix
    """
    return np.conjugate(np.transpose(M))


def F(x, sqrtS, qVec, kVec, epsVec, data, target):
    """
    Calculate the F amplitude.

    Parameters
    ----------
    x: float
        cos(theta) where theta is the scattering angle
    sqrtS: float
        The square root of the Mandelstam variable S
    qVec: numpy.ndarray
        The pion momentum vector
    kVec: numpy.ndarray
        The photon momentum vector
    epsVec: numpy.ndarray
        The photon polarization vector
    data: SaidPoles
        Object containing pole data for the calculation
    target: str
        Valid targets are "p12", "n12", "32q"

    Returns
    -------
    numpy.ndarray
        2x2 complex matrix for the F amplitude
    """
    # in units of fm^-1
    # xs = data.Epole["plus"][target][0][0]

    F1 = getF(x, sqrtS, data, 1, target=target)
    F1Term = matDotVec(sigVec, epsVec) * F1 * 1j

    qAbs = vecAbs(qVec)
    kAbs = vecAbs(kVec)

    # Match Fortran: only compute F2-F4 terms if qAbs != 0
    if qAbs != 0.0:
        F2 = getF(x, sqrtS, data, 2, target=target)
        F2Term = (
            matDotVec(sigVec, qVec) @ matDotVec(sigVec, np.cross(kVec, epsVec)) * F2
        )  # operator @ is dot product
        F2Term = F2Term / (qAbs * kAbs)

        F3 = getF(x, sqrtS, data, 3, target=target)
        F3Term = 1j * matDotVec(sigVec, kVec) * np.dot(qVec, epsVec) * F3
        F3Term = F3Term / (qAbs * kAbs)

        F4 = getF(x, sqrtS, data, 4, target=target)
        F4Term = 1j * matDotVec(sigVec, qVec) * np.dot(qVec, epsVec) * F4
        F4Term = F4Term / (qAbs * qAbs)

        out = F1Term + F2Term + F3Term + F4Term
    else:
        out = F1Term
    # print("in python F function")
    # print("x", x)
    # print("sqrtS", sqrtS)
    # print("qVec", qVec)
    # print("kVec", kVec)
    # print("epsVec", epsVec)
    # print("target", target)
    #
    # printMat(F1Term, "F1Term")
    # printMat(F2Term, "F2Term")
    # printMat(F3Term, "F3Term")
    # printMat(F4Term, "F4Term")
    # print("end of python F function prints")
    return out


def getF(x, sqrtS, data, Fi, target):
    """
    Get the F_i component of the amplitude.

    Parameters
    ----------
    x: float
        cos(theta) where theta is the scattering angle
    sqrtS: float
        The square root of the Mandelstam variable S
    data: SaidPoles
        Object containing pole data for the calculation
    Fi: int
        Which F value is being used, i=1,2,3,4
    target: str
        Target identifier ("p12", "n12", "32q")

    Returns
    -------
    complex
        The F_i value for the given parameters
    """
    out = 0
    for ell in range(6):  # Fixed to iterate through all ell
        Eplus, Mplus, Eminus, Mminus = getPoles(data, target, ell, sqrtS)

        # Debug prints for all ell values
        match Fi:
            case 1:
                legP1 = legP(x, ell + 1, deriv=1)
                tmpF = (ell * Mplus + Eplus) * legP1
                legP2 = legP(x, ell - 1, deriv=1)
                tmpF += ((ell + 1) * Mminus + Eminus) * legP2
                # print("Fi=",Fi,"ell=",ell,"tmpF=",tmpF)
                # print("legP1=",legP1,"legP2=",legP2)
                # print("Fi=1")
                # print("tmpF=", tmpF)
                out += tmpF
            case 2:
                tmpF = ((ell + 1) * Mplus + (ell * Mminus)) * legP(x, ell, deriv=1)
                # print("Fi=",Fi,"tmpF=",tmpF)
                # print("Fi=2")
                # print("tmpF=", tmpF)
                out += tmpF
            case 3:
                tmpF = (Eplus - Mplus) * legP(x, ell + 1, deriv=2)
                tmpF += (Eminus + Mminus) * legP(x, ell - 1, deriv=2)
                # print("Fi=",Fi,"tmpF=",tmpF)
                # print("Fi=3")
                # print("tmpF=", tmpF)
                out += tmpF
            case 4:
                tmpF = (Mplus - Eplus - Mminus - Eminus) * legP(x, ell, deriv=2)
                # print("Fi=",Fi,"tmpF=",tmpF)
                # print("Fi=4")
                # print("tmpF=", tmpF)
                out += tmpF
    return out


def getPoles(data, target, ell, sqrtS):
    """
    Retrieve the E and M amplitudes for a given target, angular momentum, and energy.

    Parameters
    ----------
    data : SaidPoles
        Object containing the amplitude data
    target : str
        One of 'p12', 'n12', '32q'
    ell : int
        Angular momentum (0 for S-wave, 1 for P-wave, etc.)
    sqrtS : float
        Center-of-mass energy in MeV

    Returns
    -------
    tuple
        (Eplus, Mplus, Eminus, Mminus) - Tuple of complex amplitude values
    """
    # print(f"Looking for {target}, ell={ell}, sqrtS={sqrtS}")
    Eplus = data.Epole["plus"][target]
    Mplus = data.Mpole["plus"][target]
    Eminus = data.Epole["minus"][target]
    Mminus = data.Mpole["minus"][target]

    poles = [Eplus, Mplus, Eminus, Mminus]
    names = ["Eplus", "Mplus", "Eminus", "Mminus"]
    tmp = np.zeros(4, dtype=complex)

    for i, pole in enumerate(poles):
        # print(f"Looking for {names[i]} amplitude...")
        arr = pole.get(ell, None)
        if isinstance(arr, type(None)):
            # print(f"  {names[i]}: No data found for ell={ell}")
            tmp[i] = 0.0j
        else:
            sqrtSs = np.real(arr[0])  # Energy values
            amps = arr[1]  # Amplitude values
            if len(sqrtSs) == 0:
                # print(f"  {names[i]}: Empty array for ell={ell}")
                tmp[i] = 0.0j
            elif len(sqrtSs) == 1:
                # Only one point, use it directly
                tmp[i] = amps[0]
            else:
                # Linear interpolation matching Fortran behavior
                if sqrtS <= sqrtSs[0]:
                    # Below minimum: extrapolate using first two points
                    idx1, idx2 = 0, 1
                elif sqrtS >= sqrtSs[-1]:
                    # Above maximum: use last point (no extrapolation)
                    tmp[i] = amps[-1]
                    continue
                else:
                    # Find bracketing points
                    idx2 = np.searchsorted(sqrtSs, sqrtS)
                    idx1 = idx2 - 1

                # Linear interpolation
                sqrtS1, sqrtS2 = sqrtSs[idx1], sqrtSs[idx2]
                amp1, amp2 = amps[idx1], amps[idx2]

                if abs(sqrtS2 - sqrtS1) < 1e-10:
                    tmp[i] = amp1
                else:
                    t = (amp1 - amp2) / (sqrtS1 - sqrtS2)
                    tmp[i] = amp1 + t * (sqrtS - sqrtS1)
    Eplus, Mplus, Eminus, Mminus = tmp

    return Eplus, Mplus, Eminus, Mminus


def sqrtS_to_labE(sqrtS):
    """
    Convert center-of-mass energy to lab-frame photon energy.

    Parameters
    ----------
    sqrtS : float
        Center-of-mass energy in MeV

    Returns
    -------
    float
        Lab-frame photon energy E_lab in MeV
    """
    return (sqrtS**2 - mN**2) / (2.0 * mN)


def viewData(data, target, ell, sqrtS):
    """
    Print diagnostic information about amplitude data for given parameters.

    Shows the sqrtS, target, EG (lab energy), spin values, and relevant
    lines from the said-SM22.txt file.

    Parameters
    ----------
    data : SaidPoles
        Object containing the amplitude data
    target : str
        One of 'p12', 'n12', '32q'
    ell : int
        Angular momentum (0 for S-wave, 1 for P-wave, etc.)
    sqrtS : float
        Center-of-mass energy in MeV
    """

    # Convert sqrtS to lab energy
    EG = sqrtS_to_labE(sqrtS)

    print("=" * 60)
    print(f"viewData: Amplitude lookup information")
    print("=" * 60)
    print(f"  sqrtS  = {sqrtS:.4f} MeV")
    print(f"  EG     = {EG:.4f} MeV (lab frame)")
    print(f"  target = {target}")
    print(f"  ell    = {ell}")
    print()

    # Determine isospin from target
    if target == "32q":
        I = 1.5
    else:
        I = 0.5

    # Build spin strings for all 4 amplitude types
    amp_info = []
    for plusMinus in ["plus", "minus"]:
        for ampType in ["E", "M"]:
            # For target p12/n12, use p/n prefix; for 32q, use p prefix
            if target == "p12":
                subChan = "p" + ampType
            elif target == "n12":
                subChan = "n" + ampType
            else:  # 32q
                subChan = "p" + ampType  # Use 'p' for I=3/2

            spinStr = buildspinstring(plusMinus, ell, I, subChan)
            if spinStr is not None:
                spinStr = spinStr.upper()
                # Format as "[D13 nM]" style
                # Extract wave label and subchan
                wave_label = spinStr[:-2].upper()
                sub = spinStr[-2:].lower()
                display = f"[{wave_label} {sub}]"
            else:
                display = "[N/A]"

            amp_info.append((plusMinus, ampType, display, spinStr))

    # Print amplitude types being looked up
    print("Amplitude types:")
    for plusMinus, ampType, display, spinStr in amp_info:
        # Get the data array
        if ampType == "E":
            arr = data.Epole[plusMinus][target].get(ell, None)
        else:
            arr = data.Mpole[plusMinus][target].get(ell, None)

        if arr is None or arr.shape[1] == 0:
            status = "NO DATA"
            range_str = ""
        else:
            sqrtSs = np.real(arr[0])
            EGs = [sqrtS_to_labE(s) for s in sqrtSs]
            min_EG, max_EG = min(EGs), max(EGs)
            if sqrtS < sqrtSs[0]:
                status = "BELOW RANGE"
            elif sqrtS > sqrtSs[-1]:
                status = "ABOVE RANGE"
            else:
                status = "IN RANGE"
            range_str = f"  (EG: {min_EG:.1f} - {max_EG:.1f})"

        name = f"{ampType}{plusMinus}"
        print(f"  {name:12s} {display:12s} {status}{range_str}")

    print()

    # Read the said-SM22.txt file and find relevant lines
    said_file = data.filename
    if not os.path.isfile(said_file):
        # Try in parent directory
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        said_file = os.path.join(parent_dir, data.filename)

    if os.path.isfile(said_file):
        print(f"Source file: {said_file}")
        print()

        with open(said_file, "r") as f:
            lines = f.readlines()

        # Find sections matching our parameters
        for plusMinus, ampType, display, spinStr in amp_info:
            if spinStr is None:
                continue

            # Search for the header line
            # Format in file: "PI0P S11  pE" or "PI0N S11  nE"
            wave_label = spinStr[:-2].upper()
            sub = spinStr[-2:].lower()  # Lowercase for comparison

            found_section = False
            section_start = None
            header_line = None

            for i, line in enumerate(lines):
                stripped = line.strip()
                if stripped.startswith("PI0P") or stripped.startswith("PI0N"):
                    parts = stripped.split()
                    if len(parts) >= 3:
                        file_wave = parts[1].upper()
                        file_sub = parts[2].lower()  # Lowercase for comparison
                        if file_wave == wave_label and file_sub == sub:
                            found_section = True
                            section_start = i
                            header_line = i + 1  # Next line is header
                            break

            if found_section:
                print(f"{display} found at line {section_start + 1}:")

                # Find the data lines closest to our EG
                data_lines = []
                for j in range(header_line + 1, min(header_line + 60, len(lines))):
                    dline = lines[j].strip()
                    if not dline:
                        break
                    if dline.startswith("GWU") or dline.startswith("PI0"):
                        break
                    parts = dline.split()
                    if len(parts) >= 3:
                        try:
                            file_EG = float(parts[0])
                            data_lines.append((j + 1, file_EG, lines[j].rstrip()))
                        except ValueError:
                            break

                if data_lines:
                    # Find closest lines to our EG
                    closest_idx = min(
                        range(len(data_lines)), key=lambda i: abs(data_lines[i][1] - EG)
                    )

                    # Print header
                    if header_line < len(lines):
                        print(f"  {lines[header_line].rstrip()}")

                    # Print a few lines around the closest match
                    start_idx = max(0, closest_idx - 1)
                    end_idx = min(len(data_lines), closest_idx + 2)

                    for idx in range(start_idx, end_idx):
                        line_num, file_EG, content = data_lines[idx]
                        marker = " <--" if idx == closest_idx else ""
                        print(f"  (line {line_num}) {content}{marker}")
                else:
                    print("  No data lines found in section")
                print()
            else:
                print(f"{display}: Section not found in file")
                print()

    # Also print the interpolated values
    print("Interpolated amplitude values:")
    Eplus, Mplus, Eminus, Mminus = getPoles(data, target, ell, sqrtS)
    Eplus = Eplus / data.unitsFactor
    Mplus = Mplus / data.unitsFactor
    Eminus = Eminus / data.unitsFactor
    Mminus = Mminus / data.unitsFactor
    print(f"  Eplus  = {Eplus}")
    print(f"  Mplus  = {Mplus}")
    print(f"  Eminus = {Eminus}")
    print(f"  Mminus = {Mminus}")
    print("\n\n")


def parseSpinString(spinString):
    """
    Parse a spin string into its components.

    Parameters
    ----------
    spinString: str
        String like 'S11pE', 'P33nM', 'H111pE', or 'H311nE' to be parsed

    Returns
    -------
    tuple or None
        (plusMinus, ell, I, subChan='pE'/'nM' etc.) or None if parse fails

        plusMinus ∈ {'plus','minus'}
        ell ∈ {0,1,2,3,...} from S,P,D,F,...
        I = 0.5 or 1.5 (usually)
        subChan = the substring after the partial-wave block, e.g. 'pE' or 'nE'
    """
    # Extract the angular momentum letter
    letter = spinString[0].upper()  # 'S','P','D','F',...

    # Map letter -> L
    L_map = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5}
    if letter not in L_map:
        return None
    ell = L_map[letter]

    # Extract the numerical part (could be 11, 111, 31, 311, etc.)
    match = re.match(r"[SDPFGH](\d+)([a-zA-Z][a-zA-Z])", spinString)
    if not match:
        return None

    # The numerical part could now be like '11', '111', '31', '311'
    num_part = match.group(1)
    subChan = match.group(2)  # e.g. 'pE', 'nE', etc.

    # Handle cases with 2-digit J values (like H111, H311)
    if len(num_part) == 3:
        # For formats like H111, H311 - first digit is 2*I, last two digits are 2*J
        twoI = int(num_part[0])
        twoJ = int(num_part[1:3])
    elif len(num_part) == 2:
        # For standard format like S11, P33 - first digit is 2*I, second digit is 2*J
        twoI = int(num_part[0])
        twoJ = int(num_part[1])
    else:
        return None

    I = 0.5 * twoI
    J = 0.5 * twoJ

    # plus minus if J - L = ± 1/2
    diff = J - ell
    if abs(diff - 0.5) < 1e-9:
        plusMinus = "plus"
    elif abs(diff + 0.5) < 1e-9:
        plusMinus = "minus"
    else:
        return None

    return (plusMinus, ell, I, subChan)


def buildspinstring(plusminus, ell, i, subchan):
    """
    Convert component values to a spin string (inverse of parseSpinString).

    Parameters
    ----------
    plusminus: str
        Either "plus" or "minus"
    ell: int
        The partial wave (0 for S-wave, 1 for P-wave, etc.)
    i: float
        The isospin, typically 0.5 or 1.5
    subchan: str
        The subchannel identifier, e.g. 'pm', 'ne'

    Returns
    -------
    str or None
        Spin string like 'p33nm', 's11pe', 'h111pe', or None if inputs are unphysical
    """
    # invert ell-> letter
    l_map_inv = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}

    if ell not in l_map_inv:
        return None
    letter = l_map_inv[ell]

    # 2i
    twoi = int(round(2 * i))
    # if we only allow i=0.5 or 1.5

    # plus => j = ell+0.5 => twoj= 2ell+1
    # minus => j= ell-0.5 => twoj= 2ell-1
    if plusminus == "plus":
        twoj = 2 * ell + 1
    elif plusminus == "minus":
        twoj = 2 * ell - 1
    else:
        return None

    if twoj < 1:
        return None  # e.g. 'minus' for ell=0 => unphysical

    # partial wave label e.g. 'p33' or 'h111'
    # Use special case for double-digit J values
    if twoj >= 10:
        pw = f"{letter}{twoi}{twoj}"
    else:
        pw = f"{letter}{twoi}{twoj}"

    # subchan appended
    out = pw + subchan

    badstrings = [
        "s11pm",
        "s11nm",
        "s31nm",
        "s31pm",
        "p11pe",
        "p11ne",
        "p31pe",
        "p31ne",
    ]
    if out in badstrings:
        print("you are attempting to parse an unphysical pole")
    return out


def lab_E_sqrtS(E_lab):
    """
    Convert lab-frame photon energy to center-of-mass energy.

    Parameters
    ----------
    E_lab : float
        Photon energy in laboratory frame (MeV)

    Returns
    -------
    float
        Center-of-mass energy sqrt(s) = sqrt(mN^2 + 2*mN*E_lab),
        assuming the target nucleon is at rest (proton or neutron ~ same mass)
    """
    return np.sqrt(mN**2 + 2.0 * mN * E_lab)


class Poles:
    """
    A class for loading E and M polarization data from Rijneveen

    Attributes
    ----------
    method : str
        Method used for calculations (e.g., "covariant")
    order : str
        Order of the calculation (e.g., "qTo3epsTo2")
    folderPath : str
        Path to the data folder
    unitsFactor : float
        Factor to convert units
    signs : list
        List of sign values ["plus", "minus"]
    targets : list
        List of target values ["p12", "n12", "32q"]
    Epole : dict
        Dictionary storing E polarization data
    Mpole : dict
        Dictionary storing M polarization data

    Notes
    -----
    After loading:
    - pol_dict[sign][target][i] will be a 2xN complex128 array: arr[0] = x-values, arr[1] = y-values.
    - pol_dict["filename"][sign][target][i] will be a list of filenames used to produce that array.
      If imaginary data is appended, its filename will be appended to this list.
    """

    def __init__(self, method="covariant", order="qTo3epsTo2", unitsFactor=None):
        """
        Initialize the Poles class.

        Parameters
        ----------
        method : str, optional
            Method used for calculations, default is "covariant"
        order : str, optional
            Order of the calculation, default is "qTo3epsTo2"
        unitsFactor : float, optional
            Factor to convert units, default is (10**-3)/mpiPlus
        """
        self.method = method
        self.order = order
        self.folderPath = os.path.join(self.method, "data_" + self.order)
        if unitsFactor is None:
            self.unitsFactor = (10**-3) / mpiPlus
        else:
            self.unitsFactor = unitsFactor

        self.signs = ["plus", "minus"]
        self.targets = ["p12", "n12", "32q"]

        # Initialize data dictionaries
        self.Epole = {sign: {t: {} for t in self.targets} for sign in self.signs}
        # self.Epole["filename"] = {
        #     sign: {t: {} for t in self.targets} for sign in self.signs
        # }

        self.Epole["filenames"] = []
        self.Mpole = {sign: {t: {} for t in self.targets} for sign in self.signs}
        # self.Mpole["filename"] = {
        #     sign: {t: {} for t in self.targets} for sign in self.signs
        # }
        self.Mpole["filenames"] = []

        self.get_files()

    def get_files(self):
        """
        Get all data files and load polarization data.

        Looks for files in the folderPath directory, excluding those with 'err' in the name,
        and loads the E and M polarization data.
        """
        # Get all files excluding those with 'err'
        onlyfiles = [
            f for f in listdir(self.folderPath) if isfile(join(self.folderPath, f))
        ]
        onlyfiles = [f for f in onlyfiles if "err" not in f]
        onlyfiles.sort()

        self.load_polarization_data("E", self.Epole, onlyfiles)
        self.apply_im_data("E", self.Epole)

        self.load_polarization_data("M", self.Mpole, onlyfiles)
        self.apply_im_data("M", self.Mpole)

    def load_polarization_data(self, letter, pol_dict, onlyfiles):
        """
        Load polarization data (E or M) for all sign-target-(0,1) combinations.

        Parameters
        ----------
        letter : str
            'E' or 'M' indicating the type of polarization data
        pol_dict : dict
            Dictionary to store the polarization data
        onlyfiles : list
            List of filenames to process

        Notes
        -----
        If a file is found, store the array. If not, store an empty array.
        Also record filenames in pol_dict["filenames"].
        """
        for target in self.targets:
            for sign in self.signs:
                for i in [0, 1]:
                    filename = f"{letter}{i}{sign}{target}{self.order}.dat".replace(
                        "qq", "q"
                    )
                    fullpath = join(self.folderPath, filename)
                    if filename in onlyfiles:
                        arr = self.parse_data_file(fullpath)
                        self.assert_spacing(arr[0])
                        pol_dict[sign][target][i] = arr
                        # pol_dict["filename"][sign][target][i] = [fullpath]
                        pol_dict["filenames"].append(
                            (sign + "_" + target + "_" + str(i), fullpath)
                        )
                    else:
                        # File not found: store empty arrays and empty filename list
                        pol_dict[sign][target][i] = np.zeros(
                            (2, 0), dtype=np.complex128
                        )
                        # pol_dict["filename"][sign][target][i] = []
                        pol_dict["filenames"].append(
                            (sign + "_" + target + "_" + str(i), None)
                        )

    def apply_im_data(self, letter, pol_dict):
        """
        Apply imaginary component data from the separate folder.

        Parameters
        ----------
        letter : str
            'E' or 'M' indicating the type of polarization data
        pol_dict : dict
            Dictionary containing the polarization data to update

        Notes
        -----
        Applies imaginary component data for plus/32q/1 combination.
        Before applying, ensures x-values match.
        Appends the imaginary data filename to pol_dict["filenames"].
        """
        folderPath2 = os.path.join(self.method, "data_im_epsTo3")
        fileName = f"im{letter}1plus32epsTo3.dat"
        fullpath = join(folderPath2, fileName)

        if not os.path.isfile(fullpath):
            return

        arr_im = self.parse_data_file(fullpath)
        self.assert_spacing(arr_im[0])

        # Target combination to apply imaginary data
        sign, target, i = "plus", "32q", 1
        if (
            sign in pol_dict
            and target in pol_dict[sign]
            and i in pol_dict[sign][target]
        ):
            original_arr = pol_dict[sign][target][i]
            if original_arr.shape[1] == arr_im.shape[1]:
                # Check that the x-values are the same before adding imaginary part
                if np.allclose(original_arr[0], arr_im[0]):
                    # Add imaginary part
                    original_arr[1] += arr_im[1] * 1j
                    # Append the imaginary filename
                    pol_dict["filenames"].append((sign + target + str(i), fullpath))
                else:
                    # If x-values don't match, handle as needed, here we do nothing
                    pass

    def parse_data_file(self, fullpath):
        """
        Parse a data file into a formatted array.

        Parameters
        ----------
        fullpath : str
            Full path to the data file

        Returns
        -------
        numpy.ndarray
            2xN array where:
            arr[0] = x-values
            arr[1] = y-values*unitsFactor
        """
        with open(fullpath, "r") as f:
            lines = f.read().strip().split("\n")

        x_values = []
        y_values = []
        for line in lines:
            parts = line.split()
            if len(parts) == 2:
                x_val = np.complex128(parts[0])
                y_val = np.complex128(parts[1])
                x_values.append(x_val)
                y_values.append(y_val)

        x_values = np.array(x_values, dtype=np.complex128)
        y_values = np.array(y_values, dtype=np.complex128) * self.unitsFactor

        arr = np.zeros((2, len(x_values)), dtype=np.complex128)
        arr[0] = x_values
        arr[1] = y_values
        return arr

    @staticmethod
    def assert_spacing(x_values):
        """
        Check that x-values have the expected spacing.

        Parameters
        ----------
        x_values : numpy.ndarray
            Array of x-values to check

        Raises
        ------
        AssertionError
            If x_values spacing is not 0.5 (when more than one point exists)
        """
        if len(x_values) > 1:
            assert np.isclose(x_values[1] - x_values[0], 0.5), (
                "Unexpected spacing in x_values"
            )


class SaidPoles:
    """
    A SAID-specific data handler for pole information.

    Attributes
    ----------
    filename : str
        Name of the SAID data file
    unitsFactor : float
        Factor to convert units
    signs : list
        List of sign values ["plus", "minus"]
    targets : list
        List of target values ["p12", "n12", "32q"]
    Epole : dict
        Dictionary storing E polarization data
    Mpole : dict
        Dictionary storing M polarization data

    Notes
    -----
    This needs to be fixed to make it work with different reactions... probably
    currently the curves are wrong.

    A SAID-specific variant of Poles that:
     - Reads the big 'said-SM22.txt' file once in __init__.
     - Stores the real, imaginary E- or M-amplitudes in self.Epole or self.Mpole
       in exactly the same style as the Poles class (dictionaries by sign/target/ell).
     - A separate function getPartialWaveData(spinString) returns the same 2xN array
       (rows = [ x-values, y-values ]) for that partial wave, with x = sqrt(s) in MeV
       and y = (real + i * imag) amplitude in MeV^-1 (converted from the data's mfm).
    """

    def __init__(self, filename="said-SM22.txt", unitsFactor=None):
        """
        Initialize the SaidPoles class.

        Parameters
        ----------
        filename : str, optional
            Name of the SAID data file, default is "said-SM22.txt"
        unitsFactor : float, optional
            Factor to convert units, default is 1/(MeVtofm*1000)

        Notes
        -----
        Parses the entire 'said-SM22.txt' file and stores data in:
          self.Epole[plus|minus][p12|n12|32q][ell]
          self.Mpole[plus|minus][p12|n12|32q][ell]
        Each entry is a 2 x N array [ [ sqrtS... ],
                                     [ amplitude(real+imag) ... ] ].

        If a partial wave is missing in the file, it becomes a zero-sized array.
        """
        self.filename = filename
        if unitsFactor is None:
            self.unitsFactor = 1 / (MeVtofm * 1000)
        else:
            self.unitsFactor = unitsFactor

        # internal storages just like Poles
        self.signs = ["plus", "minus"]
        self.targets = ["p12", "n12", "32q"]

        # Build nested dictionaries for Epole, Mpole
        # Epole[sign][target][ell] => 2 x N array
        # Mpole[sign][target][ell] => 2 x N array
        self.Epole = {s: {t: {} for t in self.targets} for s in self.signs}
        self.Mpole = {s: {t: {} for t in self.targets} for s in self.signs}

        # We'll fill these from parseFile()
        # but first, initialize them with empty 2x0 arrays for all sign/target/ell up to e.g. ell=3 if you like:
        for s in self.signs:
            for t in self.targets:
                # We'll typically only see ell=0 or 1 for S/P waves in the snippet, but let's go up to 3
                for ell in range(0, 5):
                    self.Epole[s][t][ell] = np.zeros((2, 0), dtype=np.complex128)
                    self.Mpole[s][t][ell] = np.zeros((2, 0), dtype=np.complex128)

        # Parse the text file in one pass
        self._parseFile()

    def _parseFile(self):
        """
        Parse the SAID data file and organize the data into dictionaries.

        Reads said-SM22.txt line by line, looks for blocks like:
           PI0P S11  pE   ...
             EG       EMreal   EMimag  ...
           ...
        and then stops at 'GWU Single energy values' or next partial wave.
        Stores all these lines as one partial-wave set.

        The lab-energy 'EG' is converted to sqrtS using lab_E_sqrtS().
        EMreal, EMimag are multiplied by 'unitsFactor' to yield amplitude in MeV^-1.

        Data is stored in self.Epole[...] or self.Mpole[...] depending on whether
        subChan is 'pE','nE','pM','nM', etc.
        The 'p' or 'n' letter + the isospin determines target='p12','n12',
        or '32q' (for I=3/2).

        The processing flow is:
        spinString 'S11pE' → parseSpinString → (plusMinus, ell, I, subChan='pE').
        """
        with open(self.filename, "r") as f:
            lines = f.readlines()

        idx = 0
        nLines = len(lines)

        while idx < nLines:
            line = lines[idx].strip()
            idx += 1
            if not line:
                # blank line, skip
                continue

            # Look for partial-wave header like:
            # "PI0P S11  pE        1/21/25"
            if line.startswith("PI0P") or line.startswith("PI0N"):
                # e.g. line might be "PI0P S11  pE  1/21/25"
                # split by whitespace
                parts = line.split()
                # typical structure: [ 'PI0P', 'S11', 'pE', '1/21/25' ]
                if len(parts) < 3:
                    continue

                # build spinString = 'S11pE' or 'H111pE'
                waveLabel = parts[1].strip()  # e.g. 'S11' or 'H111'
                subChan = parts[2].strip()  # e.g. 'pE'
                spinString = waveLabel + subChan  # => 'S11pE' or 'H111pE'

                # parse it
                parsed = parseSpinString(spinString)
                if parsed is None:
                    # e.g. it might be S31 repeated or something not standard
                    # we'll try to keep going
                    continue

                (plusMinus, ell, I, subChan) = parsed

                # figure out target from I and from subChan[0] in {p,n}
                # If I=1/2 => target = p12 if subChan[0]=='p' else n12
                # If I=1.5 => target = '32q'
                # (the data "S31 pE" would have I=1.5 so we do '32q')
                if abs(I - 1.5) < 1e-9:
                    # I=3/2
                    target = "32q"
                else:
                    # I=1/2
                    if subChan[0] == "p":
                        target = "p12"
                    else:
                        target = "n12"

                # The second character in subChan is E or M
                # e.g. 'pE' => we use Epole, 'nM' => Mpole, etc.
                ampType = subChan[1].upper()  # 'E' or 'M'

                # Prepare arrays to read in
                xvals = []
                yvals = []

                # Next lines: read until we see "GWU" or next partial wave or end-of-file
                while idx < nLines:
                    peek = lines[idx].strip()
                    if (
                        (not peek)
                        or peek.startswith("GWU")
                        or peek.startswith("PI0P")
                        or peek.startswith("PI0N")
                    ):
                        # we've hit the next partial wave or the single-energy block
                        break

                    # otherwise, parse the amplitude line, e.g.
                    #  "200.00      13.974     1.624    13.371     6.629    16.670"
                    cols = peek.split()
                    # typically we need at least 3 columns: EG, EMreal, EMimag
                    if len(cols) >= 3:
                        try:
                            labE = float(cols[0])  # MeV
                            re_val = float(cols[1])  # in mfm
                            im_val = float(cols[2])  # in mfm

                            sqrtS = lab_E_sqrtS(labE)
                            # amplitude in MeV^-1
                            amp = (re_val + 1j * im_val) * self.unitsFactor

                            xvals.append(sqrtS)
                            yvals.append(amp)
                        except ValueError:
                            pass

                    idx += 1

                # Now we have xvals,yvals for this partial wave
                if len(xvals) > 0:
                    arr = np.zeros((2, len(xvals)), dtype=np.complex128)
                    arr[0] = xvals  # row 0 is sqrtS
                    arr[1] = yvals  # row 1 is amplitude

                    # store in self.Epole or self.Mpole
                    if ampType == "E":
                        self.Epole[plusMinus][target][ell] = arr
                    else:
                        self.Mpole[plusMinus][target][ell] = arr

                # loop back around for more lines
                continue

            # otherwise, not a partial-wave header, keep scanning
            # e.g. lines that say "SM22 2700 MEV P(200) CHI/DP=..." or "GWU single energy"
            # we just skip them.
            continue

    def getPartialWaveData(self, spinString):
        """
        Get amplitude data for a specific partial wave.

        Parameters
        ----------
        spinString : str
            String specifying the partial wave, e.g., 'S11pE', 'S31pE', 'S11nE'

        Returns
        -------
        numpy.ndarray
            2xN array where:
            arr[0] = sqrtS values in MeV
            arr[1] = amplitude values (complex)

            Returns an empty array if no data was found for the specified partial wave.
        """
        parsed = parseSpinString(spinString)
        if parsed is None:
            # e.g. invalid wave specification
            return np.zeros((2, 0), dtype=np.complex128)

        (plusMinus, ell, I, subChan) = parsed

        # figure out target from I and from subChan[0]
        if abs(I - 1.5) < 1e-9:
            target = "32q"
        else:
            # I=0.5
            if subChan[0] == "p":
                target = "p12"
            else:
                target = "n12"

        # amplitude type from subChan[1] in {E,M}
        ampType = subChan[1].upper()

        if ampType == "E":
            return self.Epole[plusMinus][target][ell]
        else:
            return self.Mpole[plusMinus][target][ell]


def legP(x, n, deriv=0):
    """
    Calculate Legendre polynomials or their derivatives.

    Parameters
    ----------
    x: float or numpy.ndarray
        Input value(s) where -1 <= x <= 1, typically cos(theta)
    n: int
        Order of the Legendre polynomial
    deriv: int, optional
        Derivative order (0=P_n, 1=P_n', 2=P_n''), default is 0

    Returns
    -------
    float or numpy.ndarray
        Value of the Legendre polynomial or its derivative
    """
    # assert isinstance(x, np.ndarray)
    if deriv == 0:
        match n:
            # case -1:
            #     tmp = np.zero(len(x))
            #     tmp[:] = 1.0
            #     return tmp
            case 0 | -1:
                return 1
            case 1:
                return x
            case 2:
                return 1.5 * (x**2) - 0.5
            case 3:
                return 2.5 * (x**3) - 1.5 * x
            case 4:
                return (1 / 8) * (3 - 30 * x**2 + 35 * x**4)
            case 5:
                return (1 / 8) * (15 * x - 70 * x**3 + 63 * x**5)
            case 6:
                return (1 / 16) * (231 * x**6 - 315 * x**4 + 105 * x**2 - 5)
            case _:
                raise ValueError(f"legendreP not implimented for given n={n}")
    elif deriv == 1:
        match n:
            case -1:
                return 0
            case 0:
                return 0
            case 1:
                # tmp = np.zeros(len(x))
                # tmp[:] = 1.0
                # return tmp
                return 1
            case 2:
                return 3 * x
            case 3:
                return 0.5 * (15 * (x**2) - 3)
            case 4:
                return (-60 * x + (140 * (x**3))) / 8
            case 5:
                return (15 - 210 * x**2 + 315 * x**4) / 8
            case 6:
                return (693 * x**5) / 8 - (315 * x**3) / 4 + (105 * x) / 8
            case _:
                raise ValueError(f"legendreP not implimented for given n={n}")
    elif deriv == 2:
        match n:
            case -1:
                return 0
            case 0:
                return 0
            case 1:
                return 0
            case 2:
                # tmp = np.zeros(len(x))
                # tmp[:] = 3.0
                # return tmp
                # print("about to return 3")
                return 3
            case 3:
                return 15 * x
            case 4:
                return -7.5 + 52.5 * (x**2)
            case 5:
                return (-420 * x + 1260 * x**3) / 8
            case 6:
                return (3465 * x**4) / 8 - (945 * x**2) / 4 + 105 / 8
            case _:
                raise ValueError(f"legendreP not implimented for given n={n}")
    # print("deriv=", deriv)
    # print("x=", x)
    # print("n=", n)
    raise ValueError("something went badly wrong")
    # return 0


def vecAbs(v1):
    """
    Calculate the Euclidean norm (magnitude) of a vector.

    Parameters
    ----------
    v1: numpy.ndarray
        Input vector

    Returns
    -------
    float
        Magnitude of the vector
    """
    return np.sqrt(np.dot(v1, v1))


def mandT(v1, v2):
    """
    Calculate the Mandelstam variable t.

    Parameters
    ----------
    v1: numpy.ndarray
        First four-vector
    v2: numpy.ndarray
        Second four-vector

    Returns
    -------
    float
        t=(v1-v2)^2 = v1^2 + v2^2 - 2*v1·v2
    """
    t1 = fourvecdot(v1, v1)
    t2 = fourvecdot(v2, v2)
    t3 = fourvecdot(v1, v2)
    return t1 + t2 - 2 * t3


def fourvecsquare(v1):
    """
    Calculate the Minkowski square of a four-vector.

    Parameters
    ----------
    v1: numpy.ndarray
        Input four-vector

    Returns
    -------
    float
        The Minkowski square v1·v1
    """
    return fourvecdot(v1, v1)


def fourvecdot(v1, v2):
    """
    Calculate the Minkowski dot product of two four-vectors.

    Parameters
    ----------
    v1: numpy.ndarray
        First four-vector
    v2: numpy.ndarray
        Second four-vector

    Returns
    -------
    float
        The Minkowski dot product v1·v2
    """
    return v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2] - v1[3] * v2[3]


def matDotVec(matVec, vec):
    """
    Calculate the dot product of a vector of matrices with a vector.

    Parameters
    ----------
    matVec: numpy.ndarray
        Array of matrices [M_x, M_y, M_z]
    vec: numpy.ndarray
        Vector [v_x, v_y, v_z]

    Returns
    -------
    numpy.ndarray
        Result of M_x*v_x + M_y*v_y + M_z*v_z
    """
    tmp = matVec[0]
    assert np.shape(tmp) == (2, 2)
    assert np.shape(vec) == (3,)
    return matVec[0] * vec[0] + matVec[1] * vec[1] + matVec[2] * vec[2]


def labE_to_sqrtS(E_lab):
    """
    Convert lab-frame photon energy to center-of-mass energy.

    Parameters
    ----------
    E_lab: float
        Photon energy in the lab frame

    Returns
    -------
    float
        sqrt(s) = sqrt(mN^2 + 2 * mN * E_lab), assuming the target nucleon is at rest
    """
    return np.sqrt(mN**2 + 2.0 * mN * E_lab)


if __name__ == "__main__":
    main()
