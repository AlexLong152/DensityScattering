"""
@author: alexl
"""

import numpy as np
import CrossSection as cc
from matplotlib import pyplot as plt
from matplotlib import rcParams


rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    omega2012 = 60  # in MeV
    data2012 = np.array(
        [
            [40.0, 196.0, 2.0 * np.sqrt(34)],
            [55.0, 176.0, 10.0],
            [75.0, 125.0, np.sqrt(65)],
            [90.0, 135.0, np.sqrt(65)],
            [110.0, 138.0, np.sqrt(89)],
            [125.0, 177.0, 3.0 * np.sqrt(13)],
            [145.0, 178.0, 3.0 * np.sqrt(13)],
            [159.0, 193.0, 10.0],
        ],
        dtype=float,
    ).T
    Ntotmaxs = np.array([6, 8, 10, 12, 14])

    thetas = np.array([0, 40, 55, 75, 90, 110, 125, 145, 159, 180])

    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/"
    onebody_dir = twobody_dir + r"../onebody/"
    # savefolder = twobody_dir + r"../results/"
    tmp = np.array(twobody_dir.split(r"/"))
    # title = tmp[-3]
    energy = 60
    lambdaSRG = 1.880
    lambdaCut = 550
    omegaHs = np.array([10, 12, 14, 16, 18])
    markers = ["x", ",", "o", "v", "1", "*", "D"]
    out = {}

    for i, omegaH in enumerate(omegaHs):
        for theta in thetas:
            for j, Ntot in enumerate(Ntotmaxs):
                try:
                    ccVal = cc.ccForDict(
                        onebody_dir,
                        twobody_dir,
                        energy=energy,
                        angle=theta,
                        lambdaCut=lambdaCut,
                        lambdaSRG=lambdaSRG,
                        Ntotmax=Ntot,
                        omegaH=omegaH,
                    )
                except ValueError:
                    print(
                        f"Missing file for energy={energy}, angle={theta}, lambdaCut={lambdaCut}, lambdaSRG={lambdaSRG}, Ntotmax={Ntot}, omegaH={omegaH}"
                    )
                    ccVal = None
                out[(omegaH, Ntot, theta)] = ccVal

    inputDict = {}
    for theta in thetas:
        tmp = {}
        for omegaH in omegaHs:
            for Ntot in Ntotmaxs:
                tmp[(omegaH, Ntot)] = out[(omegaH, Ntot, theta)]

        inputDict[theta] = tmp
    colors = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
    ]
    for i, theta in enumerate(thetas):
        xVal, yVal, sigma = squeezeUncer(
            inputDict[theta], omegaHs, Ntotmaxs, returnFuncs=False
        )
        sigmaTheory = 0.1 * yVal
        sigmaTotal = np.sqrt(sigmaTheory**2 + sigma**2)
        if i == 0:
            plt.errorbar(
                theta,
                yVal,
                yerr=sigmaTotal,
                fmt="o",
                capsize=5,
                linestyle=None,
                c=colors[0],
                label=r"Preliminary TDA $\mathcal{O}(e^2 \delta^3)$ Calculation",
            )

        plt.errorbar(
            theta,
            yVal,
            yerr=sigmaTotal,
            fmt="o",
            capsize=5,
            linestyle=None,
            c=colors[0],
        )

    plt.errorbar(
        data2012[0],
        data2012[1],
        yerr=data2012[2],
        linestyle=None,
        fmt="o",
        capsize=5,
        c=colors[1],
        label="L. S. Myers \\textit{et al.} Experimental Result",
    )

    title = (
        r"Preliminary ${}^6\mathrm{Li}$ Compton Scattering at $"
        + str(energy)
        + r"\;\mathrm{MeV}$"
        + "\n"
        + f"Chiral SMS N4LO+3nfN2LO $\\Lambda_{{\\mathrm{{NN}}}}={lambdaCut} \\mathrm{{MeV}}$"
        + f"$,\\;\\Lambda_{{\\mathrm{{SRG}}}}={lambdaSRG} \\mathrm{{fm}}^{{-1}}$"
    )
    axisSize = 14
    plt.xticks(np.arange(0, 181, 15))
    plt.ylim(bottom=0)
    eps = 4
    plt.xlim([0 - eps, 180 + eps])
    plt.title(title, fontsize=axisSize + 2)
    plt.ylabel(
        r"$\mathrm{d\sigma/d\Omega\; [ \mathrm{nbn}/\mathrm{sr}]}$", fontsize=axisSize
    )
    plt.xlabel(r"$\theta$", fontsize=axisSize)
    plt.tight_layout()
    plt.legend()
    plt.show()


def squeezeUncer(inputDict, omegaHs, Ntotmaxs, returnFuncs=False):
    """
    inputDict[(omegaH, Ntot)]
    #it is generated from the following
    inputDict = {}
    for theta in thetas:
        tmp = {}
        for omegaH in omegaHs:
            for Ntot in Ntotmaxs:
                tmp[(omegaH, Ntot)] = CrossSectionValue

        inputDict[theta] = tmp

    Returns
    -------
    Either:
        xIntercept, y_Intercpet, sigma, f1,f2
    or
        xIntercept, y_Intercpet, sigma
    the "xIntercept" and f1, f2 are really only useful for plotting
    """
    # start at the omega values with the best binding energies
    # omegaH1 = 14
    # omegaH2 = 16
    # ys1 = np.array([inputDict[(omegaH1, Ntot)] for Ntot in Ntotmaxs])
    # ys2 = np.array([inputDict[(omegaH2, Ntot)] for Ntot in Ntotmaxs])

    omVals = np.array([14, 16])
    Ntot1 = Ntotmaxs[-1]
    Ntot2 = Ntotmaxs[-2]
    j = 0

    om1, om2 = omVals
    while True and j < 3:
        x1a, y1a = Ntot1, inputDict[(om1, Ntot1)]
        x1b, y1b = Ntot2, inputDict[(om1, Ntot2)]

        x2a, y2a = Ntot1, inputDict[(om2, Ntot1)]
        x2b, y2b = Ntot2, inputDict[(om2, Ntot2)]
        # m1 is slope of line 1
        m1 = (y1a - y1b) / (x1a - x1b)
        # m2 is slope of line 2
        m2 = (y2a - y2b) / (x2a - x2b)
        if np.sign(m1) > 0:
            loc = np.where(omegaHs == om1)[0][0]
            om1 = omegaHs[loc + 1]

        if np.sign(m2) < 0:
            loc = np.where(omegaHs == om2)[0][0]
            om2 = omegaHs[loc - 1]

        if np.sign(m1) < 0 and np.sign(m2) > 0:
            sigma = abs(y1a - y2a) / 2
            return compute_intersection(
                x1a, y1a, x1b, y1b, x2a, y2a, x2b, y2b, returnFuncs=returnFuncs
            ) + (sigma,)
        j += 1


def compute_intersection(x1a, y1a, x1b, y1b, x2a, y2a, x2b, y2b, returnFuncs=False):
    # Compute slopes
    m1 = (y1a - y1b) / (x1a - x1b)
    m2 = (y2a - y2b) / (x2a - x2b)

    # Compute y-intercepts
    b1 = y1a - m1 * x1a
    b2 = y2a - m2 * x2a

    x_intersect = (b2 - b1) / (m1 - m2)
    y_intersect = m1 * x_intersect + b1
    if returnFuncs:

        def f1(x):
            return m1 * x + b1

        def f2(x):
            return m2 * x + b2

        return x_intersect, y_intersect, f1, f2

    else:
        return x_intersect, y_intersect


if __name__ == "__main__":
    main()
