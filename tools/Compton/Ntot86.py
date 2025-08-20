# -*- coding: utf-8 -*-

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
    Ntotmaxs = np.array([6, 8, 10, 12, 14])
    thetas = np.array([40, 159])
    # print("len(thetas)=", len(thetas))

    twobody_dir = (
        "/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/2bod/86MeV/"
    )
    # onebody_dir = twobody_dir + r"../../onebody/"
    onebody_dir = (
        "/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/1bod/86MeV/"
    )
    # savefolder = twobody_dir + r"../results/"
    tmp = np.array(twobody_dir.split(r"/"))
    # title = tmp[-3]
    energy = 86
    lambdaSRGs = np.array([1.880, 2.236])
    # lambdaSRG = lambdaSRGs[0]
    lambdaCut = 500
    omegaHs = np.array([14, 16, 18])
    markers = ["x", ",", "o", "v", "1", "*", "D"]
    out = {}
    for lambdaSRG in lambdaSRGs:
        for theta in thetas:
            for i, omegaH in enumerate(omegaHs):
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
                    out[(lambdaSRG, omegaH, Ntot, theta)] = ccVal
    fig, axs = plt.subplots(len(lambdaSRGs), len(thetas), figsize=(12, 10))
    fig.suptitle(
        r"Ntot comparison for ${}^6\mathrm{Li}$ Compton Scattering at $"
        + str(energy)
        + r"\, \mathrm{MeV}$ ",
        fontsize=15,
    )
    colors = np.array(["C0", "C1", "C2", "C3", "C4", "C5"])

    # for k, ax in enumerate(axs):
    for k in range(len(thetas)):
        for a in range(len(lambdaSRGs)):
            lambdaSRG = lambdaSRGs[a]
            ax = axs[a, k]
            theta = thetas[k]
            ax.set_xlabel(r"$\mathrm{Ntot}$", fontsize=15)
            ax.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\Omega$", fontsize=15)
            ax.set_title(
                rf"\;$\theta={theta},\Lambda_{{\mathrm{{SRG}}}}={lambdaSRG}$",
                fontsize=15,
            )
            for i, omegaH in enumerate(omegaHs):
                legendPlotted = False
                for j, Ntot in enumerate(Ntotmaxs):
                    try:
                        val = out[(lambdaSRG, omegaH, Ntot, thetas[k])]
                        if not legendPlotted:
                            ax.scatter(
                                Ntot,
                                val,
                                marker=markers[i],
                                c=colors[i],
                                label="omegaH=" + str(omegaH),
                            )
                            legendPlotted = True
                        else:
                            ax.scatter(
                                Ntot,
                                val,
                                marker=markers[i],
                                c=colors[i],
                            )
                    except KeyError:
                        pass
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
