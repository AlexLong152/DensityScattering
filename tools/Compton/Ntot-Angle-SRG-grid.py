# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import CrossSection as cc
from matplotlib import pyplot as plt
from matplotlib import rcParams
import os
from nucdens import access
import pandas as pd

pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)

workdir = os.environ["HOME"] + r"/OneDrive/"
beta_webbase = "https://just-object.fz-juelich.de:9000/jikp03/densitystore-beta/"
rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
file = f"{workdir}densities_table.h5"
densdf = access.database(workdir=workdir, webbase=beta_webbase)
df = densdf.pddf
df = df.reset_index()
Li6Bind = -5332.33 * 6 / 1000


def main():
    Ntotmaxs = np.array([6, 8, 10, 12, 14])
    thetas = np.array([40, 159])
    energy = 60
    # energy = 100

    twobody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/2bod/{energy}MeV/"
    onebody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/1bod/{energy}MeV/"

    lambdaSRGs = np.array([1.880, 2.236, 3.00])
    lambdaCut = 500
    omegaHs = np.array([14, 16, 18])
    markers = ["x", ",", "o", "v", "1", "*", "D"]
    out = {}

    print("lambdaCut=", lambdaCut)
    print("energy=", energy)
    print(f"{'lambdaSRG':>10} {'omegaH':>10} {'Ntot':>8} {'theta':>8} {'ccVal':>12}")

    for theta in thetas:
        for lambdaSRG in lambdaSRGs:
            for i, omegaH in enumerate(omegaHs):
                for j, Ntot in enumerate(Ntotmaxs):
                    try:
                        ccDict = cc.ccForDict(
                            onebody_dir,
                            twobody_dir,
                            returnFull=True,
                            energy=energy,
                            angle=theta,
                            lambdaCut=lambdaCut,
                            lambdaSRG=lambdaSRG,
                            Ntotmax=Ntot,
                            omegaH=omegaH,
                        )
                        if ccDict is not None:
                            ccVal = ccDict["cc"]
                        oneHash = ccDict["hash_onebody"]
                        twoHash = ccDict["hash_twobody"]
                        #'hashname'
                        selectOne = df["hashname"] == oneHash
                        selectTwo = df["hashname"] == twoHash

                        out[(lambdaSRG, omegaH, Ntot, theta)] = ccVal
                        print(
                            f"{lambdaSRG:10} {omegaH:10} {Ntot:8} {theta:8}   {ccVal:8}"
                        )
                    except FileNotFoundError:
                        # print(
                        #     f"Missing file for energy={energy}, angle={theta}, lambdaCut={lambdaCut}, "
                        #     f"lambdaSRG={lambdaSRG}, Ntotmax={Ntot}, omegaH={omegaH}"
                        # )
                        ccVal = None
                    try:
                        oneBodyBind = df.loc[selectOne, "E[MeV]"].iloc[0]
                        twoBodyBind = df.loc[selectTwo, "E[MeV]"].iloc[0]
                        out[(lambdaSRG, omegaH, Ntot, theta, "oneBodyBind")] = (
                            oneBodyBind
                        )
                        out[(lambdaSRG, omegaH, Ntot, theta, "twoBodyBind")] = (
                            twoBodyBind
                        )
                    except UnboundLocalError:
                        oneBodyBind = 0
                        twoBodyBind = 0
                        out[(lambdaSRG, omegaH, Ntot, theta, "oneBodyBind")] = 0
                        out[(lambdaSRG, omegaH, Ntot, theta, "twoBodyBind")] = 0
                    # if oneBodyBind != 0 or twoBodyBind != 0:
                    #     print(
                    #         "lambdaSRG, omegaH,Ntot,theta,=",
                    #         lambdaSRG,
                    #         omegaH,
                    #         Ntot,
                    #         theta,
                    #     )
                    #     print(
                    #         "oneBodyBind-Li6Bind=", np.round(oneBodyBind - Li6Bind, 5)
                    #     )
                    #     print(
                    #         "twoBodyBind-Li6Bind=", np.round(twoBodyBind - Li6Bind, 5)
                    #     )
                    #     print("\n")

    # --- shared axes & figure setup ---
    fig, axs = plt.subplots(
        len(lambdaSRGs), len(thetas), sharex=True, sharey=True, figsize=(12, 10)
    )
    fig.suptitle(
        r"$N_{\mathrm{tot}}$ comparison for ${}^6\mathrm{Li}$ Compton Scattering at $"
        + str(energy)
        + r"\, \mathrm{MeV}$, "
        + r"\;$\Lambda_{\mathrm{NN}}="
        + f"{lambdaCut}"
        + r"\, \mathrm{MeV}$ ",
        fontsize=15,
        y=0.95,
    )

    colors = np.array(["C0", "C1", "C2", "C3", "C4", "C5"])
    r"""
    #This block of code plots the binding energies on top of the same plot
    for k in range(len(thetas)):
        for a in range(len(lambdaSRGs)):
            lambdaSRG = lambdaSRGs[a]
            theta = thetas[k]
            ax = axs[a, k]
            ax2 = ax.twinx()  # secondary y-axis for binding energies

            # axis labels
            ax.set_xlabel(r"$N_\mathrm{tot}$", fontsize=15)
            ax.set_ylabel(
                r"$\mathrm{d}\sigma/\mathrm{d}\Omega\; [\mathrm{nb}]$", fontsize=15
            )
            ax2.set_ylabel(r"Binding energy $E\;[\mathrm{MeV}]$", fontsize=15)

            # per-panel title
            ax.set_title(
                rf"\;$\theta={theta}^\circ,\;\Lambda_{{\mathrm{{SRG}}}}={lambdaSRG}"
                + r"\, \mathrm{fm}^{-1}$ ",
                fontsize=15,
                y=0.00,
            )

            # flags to avoid duplicate legend entries
            omega_legend_done = set()
            bind_legend_plotted = False

            for i, omegaH in enumerate(omegaHs):
                for j, Ntot in enumerate(Ntotmaxs):
                    # cross section (primary axis)
                    val = out.get((lambdaSRG, omegaH, Ntot, theta), None)
                    if val is not None:
                        lbl = (
                            r"$\omega_H=" + str(omegaH) + "$"
                            if omegaH not in omega_legend_done
                            else None
                        )
                        ax.scatter(
                            Ntot,
                            val,
                            marker=markers[i % len(markers)],
                            c=colors[i % len(colors)],
                            label=lbl,
                        )
                        omega_legend_done.add(omegaH)

                    # binding energies (secondary axis)
                    e1 = out.get((lambdaSRG, omegaH, Ntot, theta, "oneBodyBind"), None)
                    e2 = out.get((lambdaSRG, omegaH, Ntot, theta, "twoBodyBind"), None)
                    # plot as small jitter to avoid perfect overlap on same Ntot tick
                    x1 = Ntot - 0.08
                    x2 = Ntot + 0.08

                    if e1 != 0:
                        ax2.scatter(
                            x1,
                            e1,
                            marker="s",
                            s=30,
                            edgecolors="none",
                            c="k",
                            label="1-body bind" if not bind_legend_plotted else None,
                        )
                    if e2 != 0:
                        ax2.scatter(
                            x2,
                            e2,
                            marker="^",
                            s=30,
                            edgecolors="none",
                            c="0.35",
                            label="2-body bind" if not bind_legend_plotted else None,
                        )
                    if (e1 is not None) or (e2 is not None):
                        bind_legend_plotted = True

            # combine legends from both axes on the upper-right panel in the top row
            if a == 0 and k == 1:
                h1, l1 = ax.get_legend_handles_labels()
                h2, l2 = ax2.get_legend_handles_labels()
                ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=12)
    """
    for k in range(len(thetas)):
        for a in range(len(lambdaSRGs)):
            lambdaSRG = lambdaSRGs[a]
            theta = thetas[k]
            ax = axs[a, k]

            # axis labels: let outer panels show labels; others will be hidden via label_outer()
            ax.set_xlabel(r"$N_\mathrm{tot}$", fontsize=15)
            ax.set_ylabel(
                r"$\mathrm{d}\sigma/\mathrm{d}\Omega\; [\mathrm{nb}]$", fontsize=15
            )

            # per-panel title keeps theta and Lambda_SRG text
            ax.set_title(
                rf"\;$\theta={theta}^\circ,\;\Lambda_{{\mathrm{{SRG}}}}={lambdaSRG}"
                + r"\, \mathrm{fm}^{-1}$ ",
                fontsize=15,
                y=0.00,
            )

            for i, omegaH in enumerate(omegaHs):
                legendPlotted = False
                for j, Ntot in enumerate(Ntotmaxs):
                    try:
                        val = out[(lambdaSRG, omegaH, Ntot, theta)]
                        if not legendPlotted:
                            ax.scatter(
                                Ntot,
                                val,
                                marker=markers[i],
                                c=colors[i],
                                label=r"$\omega_H=" + str(omegaH) + "$",
                            )
                            legendPlotted = True
                        else:
                            ax.scatter(Ntot, val, marker=markers[i], c=colors[i])
                    except KeyError:
                        pass

    # remove spacing between subplots
    # plt.subplots_adjust(wspace=0, hspace=0)

    plt.subplots_adjust(wspace=0, hspace=0, left=0.1, right=0.95, top=0.9, bottom=0.1)
    # hide inner tick labels so things don't collide (outer axes keep theirs)
    for ax in axs.flat:
        ax.label_outer()

    # single, figure-level legend (bottom-right, inside figure bounds)
    handles, labels = axs[0, 0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.98, 0.02))

    axs[0, 1].legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),  # adjust these numbers if you want to push it out
        fontsize=12,
    )
    plt.show()


if __name__ == "__main__":
    main()
