# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
from matplotlib import rcParams
from os import listdir
from os.path import isfile, join

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
plt.rcParams.update({"font.size": 12})

sigx = np.array([[0, 1], [1, 0]])
sigy = np.array([[0, -1j], [1j, 0]])
sigz = np.array([[1, 0], [0, -1]])
sigVec = np.array([sigx, sigy, sigz])

MeVtofm = 197.3
mpi = 134.97
mpiPlus = 139.57
mProton = 938.272
mNeutron = 939.565
mN = (mProton + mNeutron) / 2


def main():
    """
    MN may have to change to target mass instead of nucleon mass
    The process gamma p to pi0 p is p12+32q results
    The process gamma n to pi0 n is 32q-p12 results

    """
    # target = "p12"
    targets = ["p12", "32q"]
    # for sqrtS in [1092, 1164, 1212]:
    for sqrtS in [1092]:
        S = sqrtS**2
        omega = (S - mN**2) / (2 * sqrtS)
        # omega = (S - mN**2) / (2*mN) # lab frame

        Epi = (S + mpi**2 - mN**2) / (2 * sqrtS)
        absQ = np.sqrt(Epi**2 - mpi**2)
        kVec = np.array([0, 0, omega])
        prefactor = 4 * np.pi * sqrtS / mN  # unitless

        """
        # In the center of mass frame with kVec=k zHat, p=-k zHat
        pVec = -1 * kVec  # cm frame
        k0 = omega
        p0 = np.sqrt(mN**2 + omega**2)

        ppVec= -1*qVec
        q0 = np.sqrt(mpi**2 + np.dot(qVec, qVec))
        pp0 = np.sqrt(mN**2 + np.dot(ppVec, ppVec))
        """

        epsVecs = np.array([[-1, -1j, 0], [1, -1j, 0]])
        epsVecs = epsVecs / np.sqrt(2)  # circularly polarized
        data = poles()

        # thetas = np.arange(0, np.pi, 0.02)
        # xs = np.cos(thetas)
        xs = np.arange(-1, 1, 0.02)

        # crossSec = np.zeros(len(thetas))
        crossSec = np.zeros(len(xs))

        for i, x in enumerate(xs):
            # theta = thetas[i]
            theta = np.arccos(x)
            qVec = np.array([0, np.sin(theta) * absQ, np.cos(theta) * absQ])
            M = 0
            for target in targets:
                for epsVec in epsVecs:
                    M += F(x, sqrtS, qVec, kVec, epsVec, data, target)  # has units fm
            M = prefactor * M  # prefactor is unitless

            # MSquare = np.trace(np.dot(dag(M), M))  # has units fm^2
            # units of vecs cancel, crossSec is in units fm^2
            MSquare = np.sum(M * np.conjugate(M))
            crossSec[i] += (vecAbs(qVec) / vecAbs(kVec)) * MSquare.real

        crossSec = 0.25 * crossSec * (1 / (64 * S * np.pi**2))

        microConvert = 10**-6
        crossSec = crossSec / (microConvert * 100)
        plt.plot(xs, crossSec, label=f"$\\sqrt{{s}}={sqrtS} \\;\\mathrm{{MeV}}$")

    ylabel = r"$\mathrm{d}\sigma/\mathrm{d}\Omega\;\left[\mathrm{\mu b}\, \mathrm{sr}^{-1}\right]$"
    xlabel = r"$\cos(\theta)$"
    plt.title(r"$\gamma p \to \pi^0 p$ Cross Section")
    # plt.ylim(0, 2.5 * 10**-6)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()

    """
    sqrtSs, ys=data.Epole["plus"][target][ell]
    #consistancy check

    fourp=np.array([p0,pVec[0],pVec[1],pVec[2]])
    fourpp=np.array([pp0,ppVec[0],ppVec[1],ppVec[2]])

    fourk=np.array([k0,kVec[0],kVec[1],kVec[2]])
    fourq=np.array([q0,qVec[0],qVec[1],qVec[2]])

    mandalT1= mandT(fourp,fourpp)
    mandalT2 = mpi**2 - 2 * (omega * Epi - vecAbs(kVec) * vecAbs(qVec) * x)
    diff=fourp-fourpp
    mandalT3=fourvecdot(diff,diff)
    diff=fourk-fourq
    mandalT4=fourvecdot(diff,diff)

    print("Epi1=",Epi)
    print("Epi2=",q0)
    print("mandalT1=", mandalT1)
    print("mandalT2=", mandalT2)
    print("mandalT3=", mandalT3)
    """


def dag(M):
    return np.conjugate(np.transpose(M))


def F(x, sqrtS, qVec, kVec, epsVec, data, target):
    # in units of fm^-1
    # xs = data.Epole["plus"][target][0][0]
    F1 = getF(x, sqrtS, data, 1, target=target)
    F1Term = matDotVec(sigVec, epsVec) * F1 * 1j

    F2 = getF(x, sqrtS, data, 2, target=target)
    F2Term = matDotVec(sigVec, qVec) * matDotVec(sigVec, np.cross(kVec, epsVec)) * F2
    F2Term = F2Term / (vecAbs(qVec) * vecAbs(kVec))

    F3 = getF(x, sqrtS, data, 3, target=target)
    F3Term = 1j * matDotVec(sigVec, kVec) * np.dot(qVec, epsVec) * F3
    F3Term = F3Term / (vecAbs(qVec) * vecAbs(kVec))

    return F1Term + F2Term + F3Term


def getF(x, sqrtS, data, i, target):
    out = 0
    for ell in np.array([0, 1]):
        # print("ell=",ell)
        # vals = data.Epole["plus"][target][ell]
        # xs=vals[0]
        # Eplus=vals[1]
        sqrtSs = data.Epole["plus"][target][ell][0]
        ind = np.argmin(abs(sqrtSs - sqrtS))
        # print("sqrtSs[ind]",sqrtSs[ind])
        Eplus = data.Epole["plus"][target][ell][1]
        Eminus = data.Epole["minus"][target][ell][1]
        Mplus = data.Mpole["plus"][target][ell][1]
        Mminus = data.Mpole["minus"][target][ell][1]
        Eplus = Eplus[ind]
        Eminus = Eminus[ind]
        Mplus = Mplus[ind]
        Mminus = Mminus[ind]

        match i:
            case 1:
                out += (ell * Mplus + Eplus) * legP(x, ell + 1, deriv=1)
                out += ((ell + 1) * Mminus + Eminus) * legP(x, ell - 1, deriv=1)
            case 2:
                if ell == 0:
                    out += 0
                else:
                    out += ((ell + 1) * Mplus + (ell * Mminus)) * legP(x, ell, deriv=1)
                # out += ((ell + 1) * Mplus + (ell * Mminus)) * legP(x, ell, deriv=1)

            case 3:
                if ell == 0:
                    out += 0
                else:
                    out += (Eplus - Mplus) * legP(x, ell + 1, deriv=2)
                    out += (Eminus + Mminus) * legP(x, ell - 1, deriv=2)

            # out += (Eplus - Mplus) * legP(x, ell + 1, deriv=2)
            # out += (Eminus + Mminus) * legP(x, ell - 1, deriv=2)
    return out


class poles:
    """
    A class for loading in the data all at once so the files don't have to
    be read many times.

    Creates two main properties that are of use self.Epole and self.Mpole.
    To access the data from either of these just use

    data=poles()
    tmp=data.Mpole["plus"][target][ell] #returns a 2x449 array
    tmp[0] # is x axis, sqrtS
    tmp[1] # are the y values, of the multipoles corrisponding to the x value in the same index

    returns
    """

    # def __init__(self, method="heavy_baryon", order="qTo3epsTo2"):
    def __init__(self, method="covariant", order="qTo3epsTo2"):
        self.method = method
        self.order = order
        self.folderPath = self.method + r"/data_" + self.order
        self.Epole = {}
        self.Mpole = {}
        self.dataLen = 449

        self.getFiles()

    def getFiles(self):
        """
        This looks really complex but this just parses the text files into
        complex128
        """
        unitsFactor = MeVtofm * 10**-3 / mpiPlus
        # self.unitsFactor = unitsFactor
        onlyfiles = [
            f for f in listdir(self.folderPath) if isfile(join(self.folderPath, f))
        ]
        onlyfiles = [f for f in onlyfiles if "err" not in f]
        onlyfiles = [f for f in onlyfiles if self.order in f]
        onlyfiles.sort()

        letter = "E"
        signs = ["plus", "minus"]
        targets = ["p12", "n12", "32q"]
        for sign in signs:
            self.Epole[sign] = {}
            for target in targets:
                tmp = self.Epole[sign]
                tmp[target] = {}

        for target in targets:
            for sign in signs:
                for i in [0, 1]:
                    filename = letter + str(i) + sign + target + self.order + ".dat"
                    filename = filename.replace("qq", "q")  #  handles weird edge case
                    if filename in onlyfiles:
                        fullpath = self.folderPath + r"/" + filename
                        with open(fullpath, "r") as f:
                            arr = f.read()
                        lines = arr.splitlines()
                        words = []
                        for line in lines:
                            words.extend(line.split())
                        words = np.array(words)

                        words = words.astype(np.complex128)
                        arr = np.zeros((2, len(words) // 2), dtype=np.complex128)
                        arr[0] = words[0::2]
                        arr[1] = words[1::2] * unitsFactor
                        self.Epole[sign][target][i] = arr
                        tmp = arr[0]
                        assert tmp[1] - tmp[0] == 0.5
                    else:
                        # print("fed zero for")
                        # print("filename=", filename)
                        fullpath = self.folderPath + r"/" + filename
                        # print("fullpath=", fullpath)
                        self.Epole[sign][target][i] = np.zeros(
                            (2, self.dataLen), dtype=np.complex128
                        )
                        # self.Epole[sign][target][i] = None

        folderPath2 = self.method + r"/data_im_epsTo3"
        fileName = f"im{letter}1plus32epsTo3.dat"
        fullpath = folderPath2 + r"/" + fileName

        with open(fullpath, "r") as f:
            arr = f.read()
        lines = arr.splitlines()
        words = []
        for line in lines:
            words.extend(line.split())

        words = np.array(words)
        words = words.astype(np.complex128)

        arr = np.zeros((2, len(words) // 2), dtype=np.complex128)
        arr[0] = words[0::2]
        arr[1] = words[1::2] * unitsFactor
        self.Epole["plus"]["32q"][1][1] += arr[1] * 1j

        tmp = arr[0]
        assert tmp[1] - tmp[0] == 0.5
        letter = "M"
        for sign in signs:
            self.Mpole[sign] = {}
            for target in targets:
                tmp = self.Mpole[sign]
                tmp[target] = {}

        for target in targets:
            for sign in signs:
                for i in [0, 1]:
                    filename = letter + str(i) + sign + target + self.order + ".dat"
                    filename = filename.replace("qq", "q")  #  handles weird edge case
                    if filename in onlyfiles:
                        fullpath = self.folderPath + r"/" + filename
                        with open(fullpath, "r") as f:
                            arr = f.read()
                        lines = arr.splitlines()
                        words = []
                        for line in lines:
                            words.extend(line.split())
                        words = np.array(words)

                        words = words.astype(np.complex128)
                        arr = np.zeros((2, len(words) // 2), dtype=np.complex128)
                        arr[0] = words[0::2]
                        arr[1] = words[1::2] * unitsFactor
                        self.Mpole[sign][target][i] = arr
                        tmp = arr[0]
                        assert tmp[1] - tmp[0] == 0.5
                    else:
                        # print("filename=", filename)
                        self.Mpole[sign][target][i] = np.zeros(
                            (2, self.dataLen), dtype=np.complex128
                        )

        folderPath2 = self.method + r"/data_im_epsTo3"
        fileName = f"im{letter}1plus32epsTo3.dat"
        fullpath = folderPath2 + r"/" + fileName

        with open(fullpath, "r") as f:
            arr = f.read()
        lines = arr.splitlines()
        words = []
        for line in lines:
            words.extend(line.split())

        words = np.array(words)
        words = words.astype(np.complex128)

        arr = np.zeros((2, len(words) // 2), dtype=np.complex128)
        arr[0] = words[0::2]
        arr[1] = words[1::2] * unitsFactor
        self.Mpole["plus"]["32q"][1][1] += arr[1] * 1j


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
    return 0


def legendrePTest():
    xs = np.arange(-1, 1, 0.01)
    for n in range(-1, 4):
        ys = np.array([legP(x, n, 0) for x in xs])
        plt.plot(xs, ys, label="n=" + str(n))
    plt.legend()
    plt.show()


def makeResPlots():
    data = poles()

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


def dot(v1, v2):
    if isinstance(v1, np.ndarray) and isinstance(v2, np.ndarray):
        out = np.dot(v1, v2)
    elif isinstance(v1, fourMomentum) and isinstance(v2, fourMomentum):
        threeVecA = v1.vec
        threeVecB = v2.vec
        out = (v1.fourVec[0] * v2.fourVec[0]) - np.dot(threeVecA, threeVecB)
    else:
        raise (ValueError)
    return out


def vecAbs(v1):
    return np.sqrt(np.dot(v1, v1))


def vecSqr(v1):
    return dot(v1, v1)


class fourMomentum:
    def __init__(self, vec, mass):
        self.vec = vec
        self.mass = mass
        self.fourVec = np.array(
            [np.sqrt(mass**2 + np.dot(vec, vec)), vec[0], vec[1], vec[2]]
        )

    @property
    def zeroth(self):
        return self.fourVec

    def square(self):
        return self.zeroth**2 - np.dot(self.vec, self.vec)


def mandT(v1, v2):
    """
    t=(v1-v2)^2 = v1^2 +v2^2 -2 v1 dot v2
    """
    t1 = fourvecdot(v1, v1)
    t2 = fourvecdot(v2, v2)
    t3 = fourvecdot(v1, v2)
    return t1 + t2 - 2 * t3


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
    data = poles()

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
    data = poles()
    _, axs = plt.subplots(2)

    xs, ys = data.Epole["plus"]["32q"][1]
    axs[0].plot(xs.real, ys.imag)

    axs[0].set_ylabel("$\\mathrm{Im}\\left[E^{3/2}_{1+}\\right] 10^{-3}/M_{\\pi+}$")

    xs, ys = data.Mpole["plus"]["32q"][1]
    axs[1].plot(xs.real, ys.imag)

    axs[1].set_ylabel("$\\mathrm{Im}\\left[M^{3/2}_{1+}\\right] 10^{-3}/M_{\\pi+}$")
    plt.show()


if __name__ == "__main__":
    # legendrePTest()
    # makeResPlots()
    # invert()
    # plotImag()
    main()
