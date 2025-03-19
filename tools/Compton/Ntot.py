# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import param_depend as pdep
from matplotlib import pyplot as plt
from matplotlib import rcParams


rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    fullPlot()


def fullPlot():
    Ntotmaxs = np.array([6, 8, 10, 12, 14])
    thetas = np.array([0, 40, 55, 75, 90, 110, 125, 145, 159, 180])
    print("len(thetas)=", len(thetas))

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
                    ccVal = pdep.ccForDict(
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
    s = 3

    # fig, axs = plt.subplots(5, 2, layout="tight", figsize=(8, 15))
    fig, axs = plt.subplots(5, 2, figsize=(8, 15), gridspec_kw={"hspace": 0})
    axs = axs.flatten()
    s = 20
    for i, ax in enumerate(axs):
        xVal, yVal, f1, f2, sigma = squeezeUncer(
            inputDict[thetas[i]], omegaHs, Ntotmaxs, returnFuncs=True
        )
        ax.scatter(xVal, yVal, c="black", marker="X")
        xsTmp = np.arange(10, xVal + 2, 0.2)
        ax.plot(xsTmp, f1(xsTmp), linestyle="--", c="black", alpha=0.4)
        ax.plot(xsTmp, f2(xsTmp), linestyle="--", c="black", alpha=0.4)
        if i % 2 == 1:
            ax.tick_params(axis="y", labelright=True, labelleft=False)
        t = thetas[i]

        thetalabel = f"$\\theta={t}^\\circ \\;$"
        ax.text(0.5, 0.02, thetalabel, transform=ax.transAxes, ha="center", va="bottom")
        ylabel = (
            r"$\mathrm{d \sigma/d\Omega}="
            + str(np.round(yVal, 1))
            + r"\pm"
            + str(np.round(sigma, 1))
            + r"\;\mathrm{n bn/sr}$"
        )
        ax.set_ylabel(ylabel)

        if ax == axs[-1] or ax == axs[-2]:
            xlabel = r"$N_{tot}$"
            ax.set_xlabel(xlabel)
        for j, omegaH in enumerate(omegaHs):
            for Ntot in Ntotmaxs:
                if Ntot == Ntotmaxs[0]:
                    ax.scatter(
                        Ntot,
                        out[(omegaH, Ntot, thetas[i])],
                        c=colors[j],
                        s=s,
                        marker=markers[j],
                        label=r"$\omega_{H}$=" + str(omegaH),
                    )
                else:
                    ax.scatter(
                        Ntot,
                        out[(omegaH, Ntot, thetas[i])],
                        c=colors[j],
                        s=s,
                        marker=markers[j],
                    )

        # ax.legend()
    axs[-1].legend()
    title = (
        r"${}^6\mathrm{Li}$ Compton Scattering"
        + "\n"
        + f"Chiral SMS N4LO+3nfN2LO $\\Lambda_{{\\mathrm{{NN}}}}={lambdaCut} \\mathrm{{MeV}}$"
        + f"$,\\;\\Lambda_{{\mathrm{{SRG}}}}={lambdaSRG} \\mathrm{{fm}}^{{-1}}$"
    )
    fig.suptitle(title)
    # fig.supylabel(r"$\mathrm{d\sigma/d\Omega\; [\mu bn]}$")
    # fig.subplots_adjust(left=0.2, right=1, bottom=0.2, top=1)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
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
