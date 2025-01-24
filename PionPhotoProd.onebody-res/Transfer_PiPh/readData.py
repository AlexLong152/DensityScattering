# -*- coding: utf-8 -*-

"""
@author: alexl

TODO: implement beyond l=1 in sum using experimental data
Note for extrapolation: E_l, M_l behave like |vec{q}|^l in the threshold region
"""

# import pdb
from pdb import set_trace
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
from matplotlib import rcParams
import os
from os import listdir
from os.path import isfile, join
# from copy import copy
# import saidData as sd

DEBUG = True
trace = set_trace if DEBUG else None

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
plt.rcParams.update({"font.size": 12})

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
    calcCross()
    # check()


def calcCross():
    # MN may have to change to target mass instead of nucleon mass
    # The process gamma p to pi0 p is p12+32q results
    # The process gamma n to pi0 n is 32q-p12 results

    nucs = "pp0"
    sqrtSvals = np.array([1100, 1164, 1187, 1227])
    sqrtSvals = np.sort(sqrtSvals)
    # neg1x_vals is the value read off the graph from the Rijneveen thesis figure 5.25,
    # this is to check the ratios of the values
    neg1x_vals = {1100: 0.6, 1227: 14.9, 1164: 6.1, 1187: 10.0}
    pos1x_vals = {1100: 0.26, 1227: 17, 1164: 3.9, 1187: 7.9}
    ylabel = r"$\mathrm{d}\sigma/\mathrm{d}\Omega\;\left[\mathrm{\mu b}\, \mathrm{sr}^{-1}\right]$"
    xlabel = r"$\cos(\theta)$"

    fig, axs = plt.subplots(
        len(sqrtSvals), 1, figsize=(7, 4.5 * len(sqrtSvals)), squeeze=False
    )
    # fig, axs = plt.subplots(len(sqrtSvals), 1)

    data1 = Poles()
    data2 = SaidPoles()

    labels = {data1: "Rijneveen Thesis", data2: "SAID Poles"}
    lines = {data1: "-", data2: "--"}

    for data in [data1, data2]:
        for idx, sqrtS in enumerate(sqrtSvals):
            xs = np.arange(-1, 1, 0.05)
            crossSec = np.zeros(len(xs))

            for i, x in enumerate(xs):
                S, kVec, qVec = getKinematics(sqrtS, x)
                # print("S=", S)
                MSquare = getMSquare(x, sqrtS, qVec, kVec, data, nucs)
                crossSec[i] = (vecAbs(qVec) / vecAbs(kVec)) * MSquare
            S = sqrtS**2
            crossSec = 0.25 * crossSec * (1 / (64 * S * np.pi**2))
            # crossSec = crossSec * S * np.pi
            crossSec = crossSec * MeVtofm**2
            crossSec = crossSec / 100  # Convert to barns
            crossSec = crossSec / (10**-6)  # Convert to microbarns

            leftVal = neg1x_vals.get(sqrtS, False)
            rightVal = pos1x_vals.get(sqrtS, False)
            if leftVal and rightVal:
                print("At sqrtS=", sqrtS)
                # print(
                #     "Real Value/my value for x=-1,1=",
                #     crossSec[0] / leftVal,
                #     crossSec[0] / rightVal,
                #     "\n",
                # )

                print(
                    "Real Value/my value for x=-1,1=",
                    (float(leftVal / crossSec[0]), float(rightVal / crossSec[0])),
                    "\n",
                )
                # print("That value inverse square is:", (calcVal / crossSec[0]) ** -2, "\n")
            axs[idx, 0].plot(
                xs,
                crossSec,
                label=f"$\\sqrt{{s}}={sqrtS} \\;\\mathrm{{MeV}}$ " + labels[data],
                linestyle=lines[data],
            )
            axs[idx, 0].set_ylabel(ylabel)

            if idx == len(sqrtSvals) - 1:
                axs[idx, 0].set_xlabel(xlabel)

            axs[idx, 0].legend()

    titles = {
        "pp0": r"$\gamma p \to \pi^0 p$ Cross Section",
        "nn0": r"$\gamma n \to \pi^0 n$ Cross Section",
    }
    fig.suptitle(titles[nucs], fontsize=14)
    plt.tight_layout()
    plt.show()


def check():
    data1 = Poles()
    data2 = SaidPoles()
    target = "32q"
    ell = 1
    sign = "plus"
    # should be M_1+^(3/2) and E_1+^(3/2
    labels = {data1: "Rijneveen Thesis", data2: "SaidPoles"}
    lines = {data1: "-", data2: "--"}
    minVal = 1080
    maxVal = 1250
    fig, ax = plt.subplots(2, 1)
    for i, data in enumerate([data1, data2]):
        Mvals = data.Mpole[sign][target][ell].T
        Evals = data.Epole[sign][target][ell].T

        # Epath = data.Epole["filename"][sign][target][ell]
        # print("Epath=", Epath)
        # Mvals = np.array([xy for xy in Mvals if xy[0] < 1250 and xy[0] > 1100]).T
        # Evals = np.array([xy for xy in Evals if xy[0] < 1250 and xy[0] > 1100]).T

        Mvals = np.array(
            [xy for xy in Mvals if abs(xy[0]) > minVal and abs(xy[0]) < maxVal]
        ).T
        Evals = np.array(
            [xy for xy in Evals if abs(xy[0]) > minVal and abs(xy[0]) < maxVal]
        ).T
        # trace()

        ax[0].plot(
            Evals[0].real,
            Evals[1].real,
            label="E real " + labels[data],
            linestyle=lines[data],
        )
        ax[0].plot(
            Evals[0].real,
            Evals[1].imag,
            label="E imag " + labels[data],
            linestyle=lines[data],
        )
        ax[0].set_ylabel(f"E {ell}{sign}, {target}")
        ax[0].legend()

        ax[1].plot(
            Mvals[0].real,
            Mvals[1].real,
            label="M real " + labels[data],
            linestyle=lines[data],
        )
        ax[1].plot(
            Mvals[0].real,
            Mvals[1].imag,
            label="M imag " + labels[data],
            linestyle=lines[data],
        )
        ax[1].set_ylabel(f"M {ell}{sign}, {target}")
        ax[1].legend()

    plt.show()


def getMSquare(x, sqrtS, qVec, kVec, data, nucs):
    epsVecs = np.array([[-1, -1j, 0], [1, -1j, 0]]) / np.sqrt(2)
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
    MSquare = 0.0

    # Prefactor for the amplitude
    # prefactor = 4 * np.pi * sqrtS / mN
    prefactor = 8 * np.pi * sqrtS
    # print("prefactor=", prefactor)
    for epsVec in epsVecs:
        # Compute the amplitude Mtmp for this polarization
        Mtmp = np.zeros((2, 2), dtype=complex)
        for i, target in enumerate(targets):
            # F returns a 2x2 matrix
            Mtmp += coefs[i] * F(x, sqrtS, qVec, kVec, epsVec, data, target) * prefactor

        # Compute |Mtmp|^2 by taking Mtmp Mtmp^\dagger and then the trace
        M_conjugate_transpose = Mtmp.conjugate().T
        Mtmp_MS = np.dot(Mtmp, M_conjugate_transpose)  # This is a 2x2 matrix

        # Add the trace of Mtmp_MS to MSquare
        MSquare += np.trace(Mtmp_MS).real
        # MSquare += np.sum(Mtmp_MS).real  # both are equivalent, trace is slightly more accurate with floating point

    # MSquare now contains the sum over the two photon polarizations of |M|^2.
    return MSquare


def getKinematics(sqrtS, x):
    # theta = np.arccos(x)
    S = sqrtS**2
    omega = (S - mN**2) / (2 * sqrtS)
    kVec = np.array([0, 0, omega])

    # stu = 2 * mN**2 + mpi**2
    Epi = (S + mpi**2 - mN**2) / (2 * sqrtS)

    # Enucl = (S - mpi**2 + mN**2) / (2 * sqrtS)
    # print("Epi=", Epi)
    # print("Enucl=", Enucl)
    absQ = np.sqrt(Epi**2 - mpi**2)
    qVec = np.array([0, np.sqrt(1 - x**2), x]) * absQ

    # abspp = np.sqrt(Enucl**2 - mN**2)
    # qVec = np.array([0, np.sin(x), np.cos(x)]) * absQ

    return S, kVec, qVec


def getfullKinematics(sqrtS, x):
    # theta = np.arccos(x)
    S = sqrtS**2
    omega = (S - mN**2) / (2 * sqrtS)
    kVec = np.array([0, 0, omega])

    # stu = 2 * mN**2 + mpi**2

    Epi = (S + mpi**2 - mN**2) / (2 * sqrtS)
    Enucl = (S - mpi**2 + mN**2) / (2 * sqrtS)
    print("Epi=", Epi)
    print("Enucl=", Enucl)
    absQ = np.sqrt(Epi**2 - mpi**2)
    abspp = np.sqrt(Enucl**2 - mN**2)
    qVec = np.array([0, np.sqrt(1 - x**2), x]) * absQ
    # qVec = np.array([0, np.sin(x), np.cos(x)]) * absQ

    # omegaLab = (S - mN**2) / (2 * mN)  # lab frame

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
    print("Assertions passed: Epi==q0, Enucl==pp0, abspp==|ppVec|\n")
    return S, kVec, qVec, pVec, pVec, ppVec, k0, p0, q0, pp0, omega


def dag(M):
    return np.conjugate(np.transpose(M))


def F(x, sqrtS, qVec, kVec, epsVec, data, target):
    """
    Valid targets are
    targets = ["p12","n12", "32q"]
    """
    # in units of fm^-1
    # xs = data.Epole["plus"][target][0][0]
    F1 = getF(x, sqrtS, data, 1, target=target)
    F1Term = matDotVec(sigVec, epsVec) * F1 * 1j

    F2 = getF(x, sqrtS, data, 2, target=target)
    F2Term = matDotVec(sigVec, qVec) @ matDotVec(sigVec, np.cross(kVec, epsVec)) * F2
    F2Term = F2Term / (vecAbs(qVec) * vecAbs(kVec))

    F3 = getF(x, sqrtS, data, 3, target=target)
    F3Term = 1j * matDotVec(sigVec, kVec) * np.dot(qVec, epsVec) * F3
    F3Term = F3Term / (vecAbs(qVec) * vecAbs(kVec))

    F4Term = 0  # Contributes starting at l=2

    return F1Term + F2Term + F3Term + F4Term


def getF(x, sqrtS, data, Fi, target):
    """
    x=cos(theta)
    Fi is which F value is being use, i=1...4 in  general
    but for this its i=1,2,3
    """
    out = 0
    for ell in np.array([0, 1]):
        # ind = np.argmin(abs(sqrtSs.real - sqrtS))
        Eplus = data.Epole["plus"][target][ell][1]
        Eminus = data.Epole["minus"][target][ell][1]
        Mplus = data.Mpole["plus"][target][ell][1]
        Mminus = data.Mpole["minus"][target][ell][1]

        #############################################################
        sqrtSs = data.Epole["plus"][target][ell][0]
        if len(sqrtSs) == 0:
            Eplus = 0.0j
        else:
            ind = np.argmin(abs(sqrtSs - sqrtS))
            # print("")
            # print("For Eplus")
            # print("ind=", ind)
            # print("sqrtSs[:5], sqrtS=", sqrtSs[:5], sqrtS)
            Eplus = Eplus[ind]

        #############################################################
        sqrtSs = data.Epole["minus"][target][ell][0]
        if len(sqrtSs) == 0:
            Eminus = 0.0j
        else:
            ind = np.argmin(abs(sqrtSs - sqrtS))
            # print("")
            # print("For Eminus")
            # print("ind=", ind)
            # print("sqrtSs[:5], sqrtS=", sqrtSs[:5], sqrtS)
            Eminus = Eminus[ind]

        #############################################################
        sqrtSs = data.Mpole["plus"][target][ell][0]
        if len(sqrtSs) == 0:
            Mplus = 0.0j
        else:
            ind = np.argmin(abs(sqrtSs - sqrtS))
            # print("")
            # print("For Mplus")
            # print("ind=", ind)
            # print("sqrtSs[:5], sqrtS=", sqrtSs[:5], sqrtS)
            Mplus = Mplus[ind]

        #############################################################
        sqrtSs = data.Mpole["minus"][target][ell][0]

        if len(sqrtSs) == 0:
            Mminus = 0.0j
        else:
            ind = np.argmin(abs(sqrtSs - sqrtS))
            # print("")
            # print("For Mminus")
            # print("ind=", ind)
            # print("sqrtSs[:5], sqrtS=", sqrtSs[:5], sqrtS)
            Mminus = Mminus[ind]

        #############################################################

        # trace()
        # not val returns True for length zero lists and arrays
        """
        if Eminus == 0j:
            # assert not data.Epole["minus"][target][ell]
            assert len(data.Epole["minus"][target][ell]) == 0, "len=" + str(
                len(data.Epole["minus"][target][ell])
            )
        if Eplus == 0j:
            # assert not data.Epole["plus"][target][ell]
            assert len(data.Epole["plus"][target][ell]) == 0
        if Mplus == 0j:
            # assert not data.Mpole["plus"][target][ell]
            assert len(data.Mpole["plus"][target][ell]) == 0
        if Mminus == 0j:
            # assert not data.Mpole["minus"][target][ell]
            assert len(data.Mpole["minus"][target][ell]) == 0
        """
        #############################################################
        # a = input("waiting")
        # print(a)
        # print(100 * "#")
        # print(100 * "#")
        match Fi:
            case 1:
                out += (ell * Mplus + Eplus) * legP(x, ell + 1, deriv=1)
                out += ((ell + 1) * Mminus + Eminus) * legP(x, ell - 1, deriv=1)
            case 2:
                out += ((ell + 1) * Mplus + (ell * Mminus)) * legP(x, ell, deriv=1)
            case 3:
                out += (Eplus - Mplus) * legP(x, ell + 1, deriv=2)
                out += (Eminus + Mminus) * legP(x, ell - 1, deriv=2)

    return out


def parseSpinString(spinString):
    """
    Parse something like 'S11pE' or 'P33nM' into
      (plusMinus, ell, I, subChan='pE'/'nM' etc.)

    Returns None if something is unphysical or parse fails.

    plusMinus ∈ {'plus','minus'}
    ell ∈ {0,1,2,3,...} from S,P,D,F,...
    I = 0.5 or 1.5 (usually)
    subChan = the substring after the partial-wave block, e.g. 'pE' or 'nE'

    List of bad spinStrings:
    """
    badStrings = [
        "S11pM",
        "S11nM",
        "S31nM",
        "S31pM",
        "P11pE",
        "P11nE",
        "p31pE",
        "p31nE",
    ]
    if spinString in badStrings:
        print("You are attempting to parse an unphysical pole")

    if len(spinString) < 4:
        return None

    # partial wave label is first 3 chars, e.g. 'S11', 'P33', ...
    pw = spinString[:3]  # e.g. 'S11'
    subChan = spinString[3:]  # e.g. 'pE', 'nE', etc.

    letter = pw[0].upper()  # 'S','P','D','F',...
    try:
        twoI = int(pw[1])  # e.g. 1 => I=1/2, or 3 => I=3/2
        twoJ = int(pw[2])  # e.g. 1 => J=1/2, etc.
    except ValueError:
        return None

    # Map letter -> L
    L_map = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}
    if letter not in L_map:
        return None
    ell = L_map[letter]

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


def buildSpinString(plusMinus, ell, I, subChan):
    """
    multipole to spinstring, P33pm -> M_{1+}^{3/2}
    The inverse of parseSpinString. Given plusMinus, ell, I, and subChan ('nM','pE', etc.),
    build a single string 'P33nM' or 'S11pE' etc.
    Return None if something is unphysical.

    Parameters
    -----------
    plusMinus: str
        either "plus" or "minus"
    ell: int
        The partial wave
    I: int
        the isospin
    subChan: str
        the subchannel pM, nE for example
    """
    # invert ell-> letter
    L_map_inv = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G"}

    if ell not in L_map_inv:
        return None
    letter = L_map_inv[ell]

    # 2I
    twoI = int(round(2 * I))
    # if we only allow I=0.5 or 1.5

    # plus => J = ell+0.5 => twoJ= 2ell+1
    # minus => J= ell-0.5 => twoJ= 2ell-1
    if plusMinus == "plus":
        twoJ = 2 * ell + 1
    elif plusMinus == "minus":
        twoJ = 2 * ell - 1
    else:
        return None

    if twoJ < 1:
        return None  # e.g. 'minus' for ell=0 => unphysical

    # partial wave label e.g. 'P33'
    pw = f"{letter}{twoI}{twoJ}"

    # subChan appended
    out = pw + subChan

    badStrings = [
        "S11pM",
        "S11nM",
        "S31nM",
        "S31pM",
        "P11pE",
        "P11nE",
        "p31pE",
        "p31nE",
    ]
    if out in badStrings:
        print("You are attempting to parse an unphysical pole")
    return out


def lab_E_sqrtS(E_lab):
    """
    Convert lab-frame photon energy E_lab -> center-of-mass energy sqrt(s).
    sqrt(s) = sqrt(mN^2 + 2*mN*E_lab),
    assuming the target nucleon is at rest (proton or neutron ~ same mass).
    """
    return np.sqrt(mN**2 + 2.0 * mN * E_lab)


class Poles:
    """
    A class for loading E and M polarization data.

    After loading:
    - pol_dict[sign][target][i] will be a 2xN complex128 array: arr[0] = x-values, arr[1] = y-values.
    - pol_dict["filename"][sign][target][i] will be a list of filenames used to produce that array.
      If imaginary data is appended, its filename will be appended to this list.
    """

    def __init__(self, method="covariant", order="qTo3epsTo2", unitsFactor=None):
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
        If a file is found, store the array. If not, store an empty array.
        Also record filenames in pol_dict["filename"].
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
        Apply imaginary component data from the separate folder for plus/32q/1.
        Before applying, ensure x-values match.
        Append the imaginary data filename to pol_dict["filename"][sign][target][i].
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
        Parse a data file into a 2xN array:
        arr[0] = x-values, arr[1] = y-values*unitsFactor
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
        Assert that the x_values are spaced by 0.5 if more than one point exists.
        """
        if len(x_values) > 1:
            assert np.isclose(
                x_values[1] - x_values[0], 0.5
            ), "Unexpected spacing in x_values"


class SaidPoles:
    """
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
        Parse the entire 'said-SM22.txt' file and store data in:
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
        Reads said-SM22.txt line by line, looks for blocks like:
           PI0P S11  pE   ...
             EG       EMreal   EMimag  ...
           ...
        and then stops at 'GWU Single energy values' or next partial wave.
        We'll store all these lines as one partial-wave set.

        The lab-energy 'EG' is turned into sqrtS using lab_E_sqrtS().
        EMreal, EMimag are multiplied by 'saidFactor' to yield amplitude in MeV^-1.

        Then we place it either in self.Epole[...] or self.Mpole[...] depending on whether subChan is 'pE','nE','pM','nM', etc.
        The 'p' or 'n' letter + the isospin decides target='p12','n12', or '32q' (for I=3/2).
        The spinString 'S11pE' => parseSpinString => (plusMinus, ell, I, subChan='pE').
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

                # build spinString = 'S11pE'
                waveLabel = parts[1].strip()  # e.g. 'S11'
                subChan = parts[2].strip()  # e.g. 'pE'
                spinString = waveLabel + subChan  # => 'S11pE'

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
        Return a 2xN array [ [ sqrtS... ],
                             [ amplitude... ] ]
        for the partial wave given by e.g. 'S11pE', 'S31pE', 'S11nE', etc.
        or an empty array if none was found in the file.
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
            case _:
                print("function called with n=", n)
                raise ValueError("legendreP not implimented for given n")
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

            case _:
                print("function called with n=", n)
                raise ValueError("legendreP not implimented for given n")
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
            case _:
                print("function called with n=", n)
                raise ValueError("legendreP not implimented for given n")
    raise ValueError("something went badly wrong")
    # return 0


def legendrePTest():
    xs = np.arange(-1, 1, 0.01)
    for n in range(-1, 4):
        ys = np.array([legP(x, n, 0) for x in xs])
        plt.plot(xs, ys, label="n=" + str(n))
    plt.legend()
    plt.show()


def makeResPlots():
    data = Poles()

    _, axs = plt.subplots(3, 4)
    targets = [
        "32q",
        "p12",
        "n12",
    ]
    targTex = ["$I=3/2$", "$I=p 1/2$", "$I= n1/2$"]
    for ell in [0, 1]:
        for num, t in enumerate(targets):
            vals = data.Epole["plus"][t][ell]
            axs[num, ell].plot(vals[0], vals[1])
            if num == 0:
                axs[num, ell].set_title(f"$E_{{{ell}+}}$")
            if ell == 0:
                axs[num, ell].set_ylabel(targTex[num])

    texSign = ["+", "-"]
    for i, sign in enumerate(["plus", "minus"]):
        s = texSign[i]
        for num, t in enumerate(targets):
            vals = data.Mpole[sign][t][1]
            # print(vals)
            axs[num, 2 + i].plot(vals[0], vals[1])
            if num == 0:
                axs[num, 2 + i].set_title(f"$M_{{1{s}}}$")

    plt.show()


def vecAbs(v1):
    return np.sqrt(np.dot(v1, v1))


def mandT(v1, v2):
    """
    t=(v1-v2)^2 = v1^2 +v2^2 -2 v1 dot v2
    """
    t1 = fourvecdot(v1, v1)
    t2 = fourvecdot(v2, v2)
    t3 = fourvecdot(v1, v2)
    return t1 + t2 - 2 * t3


def fourvecsquare(v1):
    return fourvecdot(v1, v1)


def fourvecdot(v1, v2):
    return v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2] - v1[3] * v2[3]


def matDotVec(matVec, vec):
    tmp = matVec[0]
    assert np.shape(tmp) == (2, 2)
    assert np.shape(vec) == (3,)
    return matVec[0] * vec[0] + matVec[1] * vec[1] + matVec[2] * vec[2]


def invert():
    ell = 1
    target = "p12"
    sqrtSs = [1092, 1164, 1212]
    sqrtS = sqrtSs[0]
    data = Poles()

    def Eplus(x):
        out = legP(x, ell) * getF(x, sqrtS, data, 1, target)
        out += -1 * legP(x, ell + 1) * getF(x, sqrtS, data, 2, target)
        out += (
            (legP(x, ell - 1) - legP(x, ell + 1))
            * getF(x, sqrtS, data, 3, target)
            * (ell / (2 * ell + 1))
        )
        out += (
            (legP(x, ell) - legP(x, ell + 2))
            * getF(x, sqrtS, data, 4, target)
            * ((ell + 1) / (2 * ell + 3))
        )
        out = out * (1 / (2 * (ell + 1)))
        return out

    testEplus = quad(Eplus, -1, 1)[0]
    # testEplus=myInte(Eplus,-1,1)
    sqrtSreadVals = data.Epole["plus"][target][ell][0]
    ind = np.argmin(abs(sqrtSreadVals - sqrtS))

    Eplusread = data.Epole["plus"][target][ell][1]
    Eplusread = Eplusread[ind]

    print("")
    print(f"Target={target},sqrtS={sqrtS},ell={ell}")
    print("---read value of E_l+=", Eplusread)
    print("integrated value E_l+=", testEplus)


def myInte(func, a, b):
    dx = 0.05
    xs = np.arange(a, b, dx)
    ys = np.array([func(x) for x in xs])
    return dx * np.sum(ys)


def plotImag():
    data = Poles()
    _, axs = plt.subplots(2)

    xs, ys = data.Epole["plus"]["32q"][1]
    axs[0].plot(xs.real, ys.imag)

    axs[0].set_ylabel("$\\mathrm{Im}\\left[E^{3/2}_{1+}\\right] 10^{-3}/M_{\\pi+}$")

    xs, ys = data.Mpole["plus"]["32q"][1]
    axs[1].plot(xs.real, ys.imag)

    axs[1].set_ylabel("$\\mathrm{Im}\\left[M^{3/2}_{1+}\\right] 10^{-3}/M_{\\pi+}$")
    plt.show()


def labE_to_sqrtS(E_lab):
    """
    Convert lab-frame photon energy E_lab -> center-of-mass energy sqrt(s).
    sqrt(s) = sqrt(mN^2 + 2 * mN * E_lab),
    assuming the target nucleon is at rest.
    """
    return np.sqrt(mN**2 + 2.0 * mN * E_lab)


if __name__ == "__main__":
    # legendrePTest()
    # makeResPlots()
    # invert()
    # plotImag()
    main()
