# -*- coding: utf-8 -*-

"""
@author: alexl

This file makes calls to the PionPhotoLib library,
many of these functions are "junk" and left over from development
"""

import numpy as np
import PionPhotoLib as ppl
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy.integrate import quad

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
plt.rcParams.update({"font.size": 12})

MeVtofm = 197.3  # 1fm =(1/197.3) MeV^-1
# TODO: mpi masses appear to matter a lot, so add sig figs, and check for mpi+ vs mpi0
mpi = 134.97
mpiPlus = 139.57
mProton = 938.272
mNeutron = 939.565
mN = (mProton + mNeutron) / 2


def main():
    checkVsData()


def checkVsData():
    """
    https://said.phys.gwu.edu/analysis/pr_analysis.html
    under "Observables"

    Observable: DSG -> (differential cross section)
    Enter independent variable:
    Acm CosA Elab -> (select Acm, A stands for angle)
    Fixed variable Elab=160

    Display jpeg of plot
    ....[Lower]......[Upper]...Limit of Plotted Data -> 0, 1
    """
    thetaCMVSCrossSec = np.array(
        [
            [0.000, 0.1007e00],
            [10.000, 0.1027e00],
            [20.000, 0.1084e00],
            [30.000, 0.1175e00],
            [40.000, 0.1293e00],
            [50.000, 0.1432e00],
            [60.000, 0.1581e00],
            [70.000, 0.1729e00],
            [80.000, 0.1866e00],
            [90.000, 0.1982e00],
            [100.000, 0.2070e00],
            [110.000, 0.2126e00],
            [120.000, 0.2151e00],
            [130.000, 0.2151e00],
            [140.000, 0.2133e00],
            [150.000, 0.2106e00],
            [160.000, 0.2080e00],
            [170.000, 0.2061e00],
            [180.000, 0.2054e00],
        ]
    )
    nucs = "pp0"
    Elab = 160
    dataExp = {
        "data": thetaCMVSCrossSec,
        "nucs": nucs,
        "Elab": Elab,
        "source": "Updated fit SM22 W. Briscoe et al, PRC108(2023)065205",
        "mass": mProton,
    }
    dataExpCheck(dataExp)


def dataExpCheck(dataExp):
    """
    Takes a dataExp dictionary object as an argument
    """
    poleData = ppl.SaidPoles()

    # omegaLab = (S - mN**2) / (2 * mN)  # lab frame
    # omegaCM = (S - mNucl**2) / (2 * sqrtS)
    mTarget = dataExp["mass"]
    nucs = dataExp["nucs"]
    S = 2 * dataExp["Elab"] * mTarget + mTarget**2
    # omegaCM = (S - mTarget**2) / (2 * np.sqrt(S))
    # print("omegaCM=", omegaCM)
    dataExp = dataExp["data"].T
    thetas = dataExp[0]

    ccs = np.zeros(len(thetas))
    for i, theta in enumerate(thetas):
        x = np.cos(theta * np.pi / 180)
        # ccs[i] = ppl.calcCrossSection(S, x, nucs, poleData)
        ccs[i] = ppl.calcCrossSectionFromRaw(S, x, nucs, poleData)

    plt.plot(thetas, ccs, label="Calculated From Poles")
    plt.plot(dataExp[0], dataExp[1], label="Experimental")
    plt.legend()
    plt.ylim(0, 1.2 * np.max(ccs))
    plt.show()
    sqrtS = np.sqrt(S)
    print("sqrtS=", sqrtS)
    for i in range(len(thetas)):
        theta = thetas[i]
        # print("theta,cc=", f"{theta:4.2f}", ccs[i])
        print(ccs[i], "--", theta)


def test_getPoles():
    """Test the getPoles function with specific parameters that match the Fortran test."""
    data = ppl.SaidPoles()
    target = "32q"
    ell = 1
    sqrtS_test = 200.0
    Eplus, Mplus, Eminus, Mminus = ppl.getPoles(data, target, ell, sqrtS_test)

    # Print the results to match Fortran output format exactly
    print("Python getPoles results:")
    print(f"Eplus = ({Eplus.real:.8E}, {Eplus.imag:.8E})")
    print(f"Mplus = ({Mplus.real:.8E}, {Mplus.imag:.8E})")
    print(f"Eminus = ({Eminus.real:.8E}, {Eminus.imag:.8E})")
    print(f"Mminus = ({Mminus.real:.8E}, {Mminus.imag:.8E})")


def basicPlot():
    nucs = "pp0"
    # nucs = "nn0"
    sqrtS = 1100
    ylabel = r"$\mathrm{d}\sigma/\mathrm{d}\Omega\;\left[\mathrm{\mu b}\, \mathrm{sr}^{-1}\right]$"
    xlabel = r"$\cos(\theta)$"
    S = sqrtS**2
    omegaLab = (S - mN**2) / (2 * mN)  # lab frame
    print("omegaLab=", omegaLab)
    data = ppl.SaidPoles()

    thetas = np.arange(0, 181, 5)
    crossSec = np.zeros(len(thetas))

    for i, theta in enumerate(thetas):
        theta_rad = theta * np.pi / 180
        x = np.cos(theta_rad)
        S, kVec, qVec = ppl.getKinematics(sqrtS, x, nucs)
        MSquare = ppl.getMSquare(sqrtS, x, nucs, data)
        crossSec[i] = (ppl.vecAbs(qVec) / ppl.vecAbs(kVec)) * MSquare
    S = sqrtS**2
    crossSec = 0.25 * crossSec * (1 / (64 * S * np.pi**2))
    # crossSec = crossSec * S * np.pi
    crossSec = crossSec * MeVtofm**2
    crossSec = crossSec / 100  # Convert to barns
    crossSec = crossSec / (10**-6)  # Convert to microbarns

    out = np.vstack((thetas, crossSec))
    print(out.T)

    plt.plot(thetas, crossSec)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.show()


def calcCrossExample():
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

    # data1 = ppl.Poles()
    data2 = ppl.SaidPoles()

    # labels = {data1: "Rijneveen Thesis", data2: "SAID Poles"}
    # lines = {data1: "-", data2: "--"}

    labels = {data2: "SAID Poles"}
    lines = {data2: "--"}
    for data in [data2]:
        for idx, sqrtS in enumerate(sqrtSvals):
            xs = np.arange(-1, 1, 0.05)
            crossSec = np.zeros(len(xs))

            for i, x in enumerate(xs):
                S, kVec, qVec = ppl.getKinematics(sqrtS, x, nucs)
                # print("S=", S)
                MSquare = ppl.getMSquare(sqrtS, x, nucs, data)
                crossSec[i] = (ppl.vecAbs(qVec) / ppl.vecAbs(kVec)) * MSquare
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
    plt.subplots_adjust(top=0.95)
    plt.show()


def check():
    data1 = ppl.Poles()
    data2 = ppl.SaidPoles()
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


def legendrePTest():
    xs = np.arange(-1, 1, 0.01)
    for n in range(-1, 4):
        ys = np.array([ppl.legP(x, n, 0) for x in xs])
        plt.plot(xs, ys, label="n=" + str(n))
    plt.legend()
    plt.show()


def makeResPlots():
    data = ppl.Poles()

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


def invert():
    def myInte(func, a, b):
        dx = 0.05
        xs = np.arange(a, b, dx)
        ys = np.array([func(x) for x in xs])
        return dx * np.sum(ys)

    ell = 1
    target = "p12"
    sqrtSs = [1092, 1164, 1212]
    sqrtS = sqrtSs[0]
    data = ppl.Poles()

    def Eplus(x):
        out = ppl.legP(x, ell) * ppl.getF(x, sqrtS, data, 1, target)
        out += -1 * ppl.legP(x, ell + 1) * ppl.getF(x, sqrtS, data, 2, target)
        out += (
            (ppl.legP(x, ell - 1) - ppl.legP(x, ell + 1))
            * ppl.getF(x, sqrtS, data, 3, target)
            * (ell / (2 * ell + 1))
        )
        out += (
            (ppl.legP(x, ell) - ppl.legP(x, ell + 2))
            * ppl.getF(x, sqrtS, data, 4, target)
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


def plotImag():
    data = ppl.Poles()
    _, axs = plt.subplots(2)

    xs, ys = data.Epole["plus"]["32q"][1]
    axs[0].plot(xs.real, ys.imag)

    axs[0].set_ylabel("$\\mathrm{Im}\\left[E^{3/2}_{1+}\\right] 10^{-3}/M_{\\pi+}$")

    xs, ys = data.Mpole["plus"]["32q"][1]
    axs[1].plot(xs.real, ys.imag)

    axs[1].set_ylabel("$\\mathrm{Im}\\left[M^{3/2}_{1+}\\right] 10^{-3}/M_{\\pi+}$")
    plt.show()


def parse_complex(s):
    """
    Helper for converting from fortran to python
    """
    # Remove any parentheses and spaces
    # print("s=", s)
    s = s.replace(" ", "")
    s = s.strip().strip("()")
    # Split the real and imaginary parts
    parts = s.split(",")
    if len(parts) == 2:
        return complex(float(parts[0]), float(parts[1]))
    else:
        # Try Python's native complex parsing as a fallback
        try:
            return complex(s)
        except ValueError:
            print("bad s=", s)
            raise ValueError


def checkFFortran():
    # Compare the results of F function in Python and Fortran
    data = ppl.SaidPoles()

    # Parse arguments from command line
    x = float(sys.argv[1])  # cos(theta)
    sqrtS = float(sys.argv[2])  # Center of mass energy

    qVec = np.array([float(x) for x in sys.argv[3:6]])
    kVec = np.array([float(x) for x in sys.argv[6:9]])

    epsVec = np.array([parse_complex(x) for x in sys.argv[9:12]])
    # print("epsVec=", epsVec)
    # print("kVec=", kVec)
    # print("qVec=", qVec)
    # Target
    target = str(sys.argv[12]).strip()

    # Parse Fortran result matrix (2x2 complex)
    fort_result = np.zeros((2, 2), dtype=complex)
    # print("sys.argv[13:17]=",sys.argv[13:17])
    fort_result[0, 0] = parse_complex(sys.argv[13])
    fort_result[0, 1] = parse_complex(sys.argv[14])
    fort_result[1, 0] = parse_complex(sys.argv[15])
    fort_result[1, 1] = parse_complex(sys.argv[16])

    whichInd = np.zeros((2, 2), dtype=int)  # this is for error handling
    whichInd[0, 0] = 13
    whichInd[1, 0] = 14
    whichInd[0, 1] = 15
    whichInd[1, 1] = 16
    # Calculate Python result
    py_result = ppl.F(x, sqrtS, qVec, kVec, epsVec, data, target)

    # print("\nPython result matrix:")
    # print(f"[{py_result[0, 0]:.8e}, {py_result[0, 1]:.8e}]")
    # print(f"[{py_result[1, 0]:.8e}, {py_result[1, 1]:.8e}]")
    #
    # print("\nFortran result matrix:")
    # print(f"[{fort_result[0, 0]:.8e}, {fort_result[0, 1]:.8e}]")
    # print(f"[{fort_result[1, 0]:.8e}, {fort_result[1, 1]:.8e}]")

    diff = py_result - fort_result
    diff2 = diff.flatten()
    worstVal = np.max([abs(x) for x in diff2])
    if worstVal < 1e-10:
        print(f"target={target} Passed\n")
    else:
        print("max(abs(python values - fortran values))=", worstVal)
        # print("\nPython - Fortran =")
        # print(f"[{diff[0, 0]:.8e}, {diff[0, 1]:.8e}]")
        # print(f"[{diff[1, 0]:.8e}, {diff[1, 1]:.8e}]")
        # Calculate element-wise differences and ratios
        for i in range(2):
            for j in range(2):
                py_val = py_result[i, j]
                fort_val = fort_result[i, j]

                if abs(fort_val) > 1e-8:  # Avoid division by very small numbers
                    ratio = py_val / fort_val

                    if abs(abs(ratio) - 1) > 0.01:
                        print(
                            f"Difference detected: |ratio - 1| = {abs(abs(ratio) - 1):.8f}"
                            f"With py_val={py_val} and fort_val={fort_val}"
                        )
                        ind = whichInd[i, j]
                        print("sys.argv[ind]=", sys.argv[ind])

                        print(100 * "!")
                else:
                    diff = abs(py_val - fort_val)

                    if diff > 1e-10:
                        print(f"Element [{i},{j}]: abs diff = {diff:.8e}")
                        print(f"difference detected: diff = {diff:.8e}")
                        print(100 * "!")

    print(75 * "#")
    print(75 * "#")


def checkPoleReadingFortran():
    # Skip the general test and just do the direct pole calculation
    data = ppl.SaidPoles()
    sqrtS = float(sys.argv[1])

    # Ensure target string has no extra spaces for dictionary lookup
    target = str(sys.argv[2]).strip()

    ell = int(sys.argv[3])

    # Handle optional Fortran comparison values if provided

    Eplus_fort = complex(sys.argv[4])
    Mplus_fort = complex(sys.argv[5])
    Eminus_fort = complex(sys.argv[6])
    Mminus_fort = complex(sys.argv[7])

    Eplus, Mplus, Eminus, Mminus = ppl.getPoles(data, target, ell, sqrtS)
    # print("Python getPoles result:")

    # print("Eplus=",Eplus)
    # print("Eminus=",Eminus)
    # print("Mplus=",Mplus)
    # print("Mminus=",Mminus)

    mylabels = {
        "Eplus": (Eplus, Eplus_fort),
        "Mplus": (Mplus, Mplus_fort),
        "Eminus": (Eminus, Eminus_fort),
        "Mminus": (Mminus, Mminus_fort),
    }
    noIssues = True
    for key, value in mylabels.items():
        pyVal, fortVal = value
        if abs(fortVal) > 1e-10:  # Avoid division by very small numbers
            ratio = abs(pyVal / fortVal)
            if abs(ratio - 1) > 0.001:
                print(f"ISSUE HERE {key}:", ratio, 20 * "-", "ISSUE HERE")
                noIssues = False
        else:
            if abs(fortVal - pyVal) > 1e-6:
                print(
                    f" ISSUE HERE {key}:",
                    "fortVal=",
                    fortVal,
                    "pyVal=",
                    pyVal,
                    10 * "-",
                    "ISSUE HERE",
                )
                noIssues = False
    if noIssues:
        print("No Issues Found")
    else:
        print("Issue Found" + 100 * "!")

    print(" " + 50 * "#")


def checkGetFFortran():
    # Compare the results of getF function in Python and Fortran
    data = ppl.SaidPoles()

    # Parse arguments from command line
    x = float(sys.argv[1])  # cos(theta)
    sqrtS = float(sys.argv[2])  # Center of mass energy
    fi = int(sys.argv[3])  # F index (1-4)
    target = str(sys.argv[4]).strip()  # Target

    # Parse Fortran result (complex number)
    fort_result = parse_complex(sys.argv[5])
    # print("In python:")
    # print("x=",x)
    # print("sqrtS=",sqrtS)
    # print("fi=",fi)
    # print("target=",target)

    # Calculate Python result
    py_result = ppl.getF(x, sqrtS, data, fi, target)

    # Compare results

    # Calculate difference and ratio
    passing = True
    if abs(fort_result) > 1e-10:  # Avoid division by very small numbers
        ratio = py_result / fort_result
        # print(f"Ratio = {ratio:.8f}")
        if abs(abs(ratio) - 1) > 0.01:
            print("fi=", fi)
            print(
                f"|ratio - 1| = {abs(abs(ratio) - 1):.8f}  " + 30 * "-" + "Issue Here"
            )
            print(f"Python result: {py_result:.8e}")
            print(f"Fortran result: {fort_result:.8e}")
            passing = False
    else:
        diff = abs(py_result - fort_result)
        # print(f"Absolute difference = {diff:.8e}")

        if diff > 1e-10:
            print("fi=", fi)
            print(f"Difference is {diff:.8e}  " + 30 * "-" + "Issue Here")
            print(f"Python result: {py_result:.8e}")
            print(f"Fortran result: {fort_result:.8e}")
        else:
            passing = False
    if passing:
        print("Passed")



if __name__ == "__main__":
    # When called with command-line arguments, run the appropriate test
    if len(sys.argv) == 1:
        main()
    else:
        function = sys.argv[1]
        # print("function=", function)
        match function:
            case "getF":
                sys.argv.remove("getF")
                checkGetFFortran()
            case "getPoles":
                sys.argv.remove("getPoles")
                checkPoleReadingFortran()
            case "checkF":
                sys.argv.remove("checkF")
                checkFFortran()
            case _:
                sysStr = "\nsys.argv=" + str(sys.argv)
                outStr = (
                    "Did not match any case \nlen(sys.argv)="
                    + str(len(sys.argv))
                    + sysStr
                )
                raise ValueError(outStr)
