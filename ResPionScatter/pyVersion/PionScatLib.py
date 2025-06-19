# -*- coding: utf-8 -*-

"""
Python library for pion-nucleon scattering calculations.

This module provides functions to calculate scattering amplitudes,
cross-sections and various kinematic quantities for pion-nucleon interactions.

@author: alexl

In this code, q, q' is the pion momentum before/after
and p, p' is the nucleon momentum before and after
"""

import numpy as np
from parseFile import getFileData
from copy import copy
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"

mpi = 134.97
mpiPlus = 139.57
mProton = 938.272
mNeutron = 939.565
mN = (mProton + mNeutron) / 2
mpiDict = {-1: mpiPlus, 0: mpi, 1: mpiPlus}
mNuc = {-1: mNeutron, 1: mProton}

dataDict = getFileData()


def main():
    sqrtS = 1162
    isospin = 1
    piCharge = 1
    thetas = np.arange(0, 181, 30)
    print("  theta      CS       Mat CS       Ratio")
    for theta in thetas:
        x = np.cos(theta * np.pi / 180)
        # f, g = getGH(sqrtS, x, isospin, piCharge)
        DSG = getcsGH(sqrtS, x, isospin, piCharge)
        # f, g = getGH(sqrtS, x, isospin=1, piCharge=1)
        MatDSG = getCS(sqrtS, x, isospin, piCharge)
        ratio = DSG / MatDSG
        print(f"{theta:7.2f}  {DSG:9.6f}   {MatDSG:9.6f},  {ratio:9.6f}")


def getCS(sqrtS, x, isospin=1, piCharge=0):
    mat = getMat(sqrtS, x, isospin, piCharge)
    matDag = dag(mat)
    CrossSec = 10 * np.trace(np.dot(mat, matDag)).real / 4

    # TODO:find where the missing factor of 2 comes from
    # should get 2*2 spins for unpolarized target, what about the rest
    # its possbile this factor should be in getMat instead
    # But then it would have to be sqrt(2)
    fudgeFactor = 2

    return CrossSec * fudgeFactor


def getMat(sqrtS, x, isospin=1, piCharge=0):
    """
    Parameters
    ----------
    sqrtS   : float
      Center‐of‐mass energy (MeV)
    x       : float
      cos(θ)
    isospin : int
      Nucleon isospin (1=proton, −1=neutron)
    piCharge: int
      Pion charge (−1, 0, +1)
    g, h = getGH(sqrtS, x, isospin, piCharge)
    """
    m1 = mpiDict[piCharge]  # pion mass (MeV)
    m2 = mNuc[isospin]  # nucleon mass (MeV)

    # q and qp are inital and final pion 3 momentum
    _, qVec, qpVec = getKinematics(sqrtS, x, m1, m2, m1, m2)

    qHat = normalize(qVec)
    qpHat = normalize(qpVec)

    matFactor = -1j * matDotVec(sigVec, np.cross(qpHat, qHat))

    mat = iden * g + matFactor * h

    return mat


def getcsGH(sqrtS, x, isospin=1, piCharge=0):
    """
    Returns cross section in milibarns
    """
    g, h = getGH(sqrtS, x, isospin, piCharge)
    sintheta = np.sqrt(1 - x**2)
    # theta = np.arccos(x)
    # print("x=", x)
    # print("theta=", theta, "=", theta * 180 / np.pi)
    # print("sqrt(1-x**2)=", np.sqrt(1 - x**2))
    # print("np.sin(theta)=", np.sin(theta))
    # assert abs(np.sin(theta) - sintheta) < 1e-10
    DSG = abs(g) ** 2 + abs(h * sintheta) ** 2
    DSG = DSG * 10
    return DSG


def getGH(sqrtS, x, isospin=1, piCharge=0):
    """
    Build the 2×2 “bare” partial‐wave amplitude F_code(s,θ) in fm, directly
    using Höhler’s eq. (A 8.23):

      F_code(s,θ) = Σ_{I=1,3} [  ( (ℓ+1) f_{ℓ+}^{(I)} + ℓ f_{ℓ−}^{(I)} )·P_ℓ(cosθ)·𝐈
                              − i (f_{ℓ+}^{(I)} − f_{ℓ−}^{(I)})·P_ℓ'(cosθ)·(σ·n̂)  ] × ℘_I

    where ℘_I is the 2×2 isospin‐projection factor for total isospin I, and
    n̂ = (q̂′×q̂)/|q̂′×q̂| is the unit vector normal to the scattering plane.

    Parameters
    ----------
    sqrtS   : float
      Center‐of‐mass energy (MeV)
    x       : float
      cos(θ)
    isospin : int
      Nucleon isospin (1=proton, −1=neutron)
    piCharge: int
      Pion charge (−1, 0, +1)

    Returns
    -------
    g, hOverSinTheta: float
        The g and h values as defined in equatin A8.21-A8.22 in Pions And Nuclei
        By Ericson and Weise. Return h/Sin(theta) instead of just h
    """
    twoIs = np.array([1, 3])
    ells = np.arange(5, dtype=int)
    charge2Idx = {1: 0, -1: 1, 0: 2}
    isoVecDict = {1: np.array([1, 0]), -1: np.array([0, 1])}
    isoVec = isoVecDict[isospin]  # e.g. [1,0] for proton, [0,1] for neutron

    # 1) fill in masses
    m1 = mpiDict[piCharge]  # pion mass (MeV)
    m2 = mNuc[isospin]  # nucleon mass (MeV)

    # 2) build unit‐vectors q̂ and q̂′ for the spin‐flip term:
    _, qVec, _ = getKinematics(sqrtS, x, m1, m2, m1, m2)
    # 4) assemble the sum over I=1,3 and ℓ=0…4 (per A 8.23)
    gTerm = 0
    hTerm = 0

    for twoI in twoIs:
        gTermTmp = 0
        hTermTmp = 0
        # each of g/hTermtmp has to be multiplied by weight before added to actual sum
        chargeIdx = charge2Idx[piCharge]  # 0→π⁺, 1→π⁻, 2→π⁰
        projII = projI[twoI][chargeIdx]  # 2×2 isospin matrix

        weight = float(np.dot(isoVec, np.dot(projII, isoVec)))
        # print(f"For twoI={twoI}, weight=", weight)

        # 4c) loop partial waves ℓ = 0…4
        for ell in ells:
            fPlus = getF(qVec, twoI, ell, 1, sqrtS)  # f_{ℓ+}^{(I)} in fm
            fMinus = getF(qVec, twoI, ell, -1, sqrtS)  # f_{ℓ−}^{(I)} in fm

            # 4d) Legendre P_ℓ(cosθ) and its derivative P_ℓ'(cosθ):
            P_l = legP(x, ell, deriv=0)
            P_lprime = legP(x, ell, deriv=1)

            # 4e) build the two pieces from A 8.23:
            #     (ℓ+1)*fPlus + ℓ*fMinus  multiplies P_ℓ → spin‐nonflip
            gTermTmp += ((ell + 1) * fPlus + ell * fMinus) * P_l

            #     −i (fPlus − fMinus) multiplies P_ℓ'(cosθ) and (σ·n̂) → spin‐flip
            hTermTmp += (fPlus - fMinus) * P_lprime

        gTerm += gTermTmp * weight
        hTerm += hTermTmp * weight

    return gTerm, hTerm


def getF(qVec, twoI, ell, sign, sqrtS):
    """
    Calculate the partial wave amplitude. Returns f_{I,l} in units of fm

    Parameters
    ----------
    qVec : numpy.ndarray
        Momentum vector
    twoI : int
        Twice the isospin value
    ell : int
        Angular momentum
    sign : int
        Sign for j=l±1/2
    sqrtS : float
        Center-of-mass energy

    Returns
    -------
    complex
        Partial wave amplitude
    """
    # TODO: interpolate S instead of taking nearest
    sqrtSTmp = sqrtS / 2
    sqrtStmp = int(np.round(sqrtSTmp))  # round to the nearest even number
    sqrtStmp = 2 * sqrtStmp
    # print("sqrtStmp=", sqrtStmp)
    # print("2*ell+sign=", 2 * ell + sign)
    try:
        deltaRe, sr = dataDict[ell][twoI][2 * ell + sign][sqrtStmp]
    except KeyError:
        deltaRe = 0
        sr = 0

    eta = np.sqrt(1 - sr)
    qAbs = vecAbs(qVec)
    fOut1 = (eta * np.e ** (2j * deltaRe)) - 1
    fOut1 = (1 / (2j * qAbs)) * fOut1  # now in units of MeV^-1
    fOut1 = fOut1 * MeVtofm
    # deltaIm = np.log(eta) / (-2)
    # delta = deltaRe + deltaIm * 1j
    # fOut2 = np.e ** (2j * delta) - 1
    # fOut2 = (1 / (2j * qAbs)) * fOut2
    # assert fOut1 == fOut2
    return fOut1


def printMat(mat, val=None):
    """
    Print a matrix with optional label.

    Parameters
    ----------
    mat : numpy.ndarray
        Matrix to print
    val : str, optional
        Label for the matrix, by default None
    """
    assert len(np.shape(mat)) == 2
    if val is not None:
        for i in range(len(mat[0])):
            for j in range(len(mat[:, 0])):
                print(f"{val}[{i},{j}]=", f"{mat[i, j]:.6E}")
    else:
        for i in range(len(mat[0])):
            for j in range(len(mat[:, 0])):
                print(f"Mat[{i},{j}]=", f"{mat[i, j]:.6E}")

    print("")


def Wcm2Tlab(sqrtS, m1=mpiPlus, m2=mProton):
    S = sqrtS**2
    # m3 = m1
    # m4 = m2
    assert m1 < m2
    E1cm = (S + m1**2 - m2**2) / (2 * sqrtS)
    E2cm = (S - m1**2 + m2**2) / (2 * sqrtS)
    assert E2cm > E1cm
    p1cm = np.sqrt((E1cm**2) - m1**2)  # pion momentum
    p1_lab = p1cm * sqrtS / m2
    E1lab = np.sqrt(m1**2 + p1_lab**2)
    Tlab = E1lab - m1
    return Tlab


def Tlab2Wcm(Tlab, m1=mpiPlus, m2=mProton):
    # S = sqrtS**2
    E1lab = Tlab + m1
    p1_lab = np.sqrt(E1lab**2 - m1**2)
    S = (E1lab + m2) ** 2 - p1_lab**2
    sqrtS = np.sqrt(S)
    # p1cm = p1_lab * m2 / sqrtS
    return sqrtS


def getKinematics(sqrtS, x, m1, m2, m3, m4):
    """
    Calculate kinematic variables for a given reaction.

    For pion photoproduction, I am programming with the following in mind
    m1: pion mass
    m2: nucleon mass
    m3: pion mass
    m4: nucleon mass

    q,q': pion 3-momenta before and after
    p,p': nucleon 3-momenta before and after

    Parameters
    ----------
    sqrtS: float
        The square root of the Mandelstam variable S
    x: float
        cos(theta) where theta is the scattering angle
    m1: float
        Mass of particle 1 (initial pion)
    m2: float
        Mass of particle 2 (initial nucleon)
    m3: float
        Mass of particle 3 (final pion)
    m4: float
        Mass of particle 4 (final nucleon)

    Returns
    -------
    S: float Mandelstam S
    qVec: float length 3 numpy array, inital pion momentum
    qpVec: float length 3 numpy array, final pion momentum
    """
    S = sqrtS**2
    # m3=m1
    # m4=m2
    assert m1 < m2
    E1cm = (S + m1**2 - m2**2) / (2 * sqrtS)
    E2cm = (S - m1**2 + m2**2) / (2 * sqrtS)
    assert E2cm > E1cm

    E3cm = (S + m3**2 - m4**2) / (2 * sqrtS)
    E4cm = (S - m3**2 + m4**2) / (2 * sqrtS)

    # print("E3 numer",(S+m3**2 - m4**2))
    # print("E4 numer",(S-m3**2 + m4**2))
    # # assert E3cm==E1cm
    # # assert E4cm==E2cm
    # # E3cm = E1cm
    # # E4cm = E2cm
    # # kVec = np.array([0, 0, omega])
    # print("m1=",m1)
    # print("m2=",m2)
    # print("m3=",m3)
    # print("m4=",m4)
    # print("E1cm=",E1cm)
    # print("E2cm=",E1cm)
    # print("E3cm=",E1cm)
    # print("E4cm=",E4cm)

    absp = np.sqrt(E1cm**2 - m1**2)
    qVec = np.array([0, 0, 1]) * absp

    # absQ = np.sqrt(E2cm**2 - m2**2)
    # pVec = np.array([0, 0, -1]) * absQ

    absQp = np.sqrt(E3cm**2 - m3**2)
    qpVec = np.array([0, np.sqrt(1 - x**2), x]) * absQp

    # abskp = np.sqrt(E4cm**2 - m4**2)
    # kpVec = np.array([0, np.sqrt(1 - x**2), x]) * abskp

    # abspp = np.sqrt(Enucl**2 - mN**2)
    # qVec = np.array([0, np.sin(x), np.cos(x)]) * absQ

    return S, qVec, qpVec


def getfullKinematics(sqrtS, x, m1, m2, m3, m4):
    """
    Calculate full set of kinematic variables.

    Parameters
    ----------
    sqrtS: float
        The square root of the Mandelstam variable S
    x: float
        cos(theta) where theta is the scattering angle
    m1: float
        Mass of particle 1 (initial pion)
    m2: float
        Mass of particle 2 (initial nucleon)
    m3: float
        Mass of particle 3 (final pion)
    m4: float
        Mass of particle 4 (final nucleon)

    Returns
    -------
    S, qVec,qpVec,kVec, kpVec
    """
    S = sqrtS**2
    # m3=m1
    # m4=m2
    assert m1 < m2
    E1cm = (S + m1**2 - m2**2) / (2 * sqrtS)
    E2cm = (S - m1**2 + m2**2) / (2 * sqrtS)
    assert E2cm > E1cm

    E3cm = (S + m3**2 - m4**2) / (2 * sqrtS)
    E4cm = (S - m3**2 + m4**2) / (2 * sqrtS)

    absp = np.sqrt(E1cm**2 - m1**2)
    qVec = np.array([0, 0, 1]) * absp

    absk = np.sqrt(E2cm**2 - m2**2)
    kVec = np.array([0, 0, -1]) * absk

    absQp = np.sqrt(E3cm**2 - m3**2)
    qpVec = np.array([0, np.sqrt(1 - x**2), x]) * absQp

    abskp = np.sqrt(E4cm**2 - m4**2)
    kpVec = np.array([0, np.sqrt(1 - x**2), x]) * -1 * abskp

    # abspp = np.sqrt(Enucl**2 - mN**2)
    # qVec = np.array([0, np.sin(x), np.cos(x)]) * absQ

    return S, qVec, qpVec, kVec, kpVec


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
    Mtmp = copy(M)
    return np.conjugate(np.transpose(Mtmp))


def legP(x, n, deriv=0):
    """
    Calculate Legendre polynomials or their derivatives.

    Parameters
    ----------
    x: float or numpy.ndarray
        x=cos(theta)
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


def normalize(vec):
    """
    Normalize a vector to unit length.

    Parameters
    ----------
    vec : numpy.ndarray
        Vector to normalize

    Returns
    -------
    numpy.ndarray
        Normalized vector
    """
    assert len(vec) == 3
    vec = copy(vec)
    return vec / np.sqrt(np.dot(vec, vec))


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


sigx = np.array([[0, 1], [1, 0]])
sigy = np.array([[0, -1j], [1j, 0]])
sigz = np.array([[1, 0], [0, -1]])
sigVec = np.array([sigx, sigy, sigz])
iden = np.array([[1, 0], [0, 1]], dtype=complex)

# P3plus = np.array([[1, 1j / 3], [-1j / 3, 1]])
# P3minus = np.array([[1, -1j / 3], [1j / 3, 1]])
P3plus = np.array([[1, 0], [0, 1 / 3]])
P3minus = np.array([[1 / 3, 0], [0, 1]])
P3neutral = np.array([[2 / 3, 0], [0, 2 / 3]])
P3 = np.array([P3plus, P3minus, P3neutral])

P1plus = np.array([[0, 0], [0, 2 / 3]])
P1minus = np.array([[2 / 3, 0], [0, 0]])
P1neutral = np.array([[1 / 3, 0], [0, 1 / 3]])
P3 = np.array([P3plus, P3minus, P3neutral])
P1 = np.array([P1plus, P1minus, P1neutral])
projI = {1: P1, 3: P3}

for i in range(3):
    assert (P3[i] + P1[i] == np.eye(2)).all()
MeVtofm = 197.3  # 1fm =(1/197.3) MeV^-1

if __name__ == "__main__":
    main()
