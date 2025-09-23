# -*- coding: utf-8 -*-

"""
@title: Li6 Compton: x = lambdaSRG, vary omegaH, Ntot=14 only
@author: alexl
"""

import numpy as np
import CrossSection as cc
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    # --- configuration ---
    Ntotmax = 14  # <- ONLY Ntot=14
    thetas = np.array([40, 159])
    energy = 100
    lambdaSRGs = np.array([1.880, 2.236, 3.0])
    lambdaCut = 500
    omegaHs = np.array([14, 16, 18])  # <- multiple ω_H values

    twobody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/2bod/{energy}MeV/"
    onebody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/1bod/{energy}MeV/"

    markers = ["o", "s", "D", "v", "^", "*"]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5"]

    # --- compute / cache ---
    out = {}
    for lam in lambdaSRGs:
        for theta in thetas:
            for omegaH in omegaHs:
                try:
                    val = cc.ccForDict(
                        onebody_dir,
                        twobody_dir,
                        energy=energy,
                        angle=theta,
                        lambdaCut=lambdaCut,
                        lambdaSRG=lam,
                        Ntotmax=Ntotmax,
                        omegaH=omegaH,
                    )
                except ValueError:
                    print(
                        f"Missing file for E={energy}, θ={theta}, Λ_NN={lambdaCut}, "
                        f"Λ_SRG={lam}, Ntot={Ntotmax}, ω_H={omegaH}"
                    )
                    val = None
                out[(lam, omegaH, theta)] = val

    # --- figure setup (independent axes) ---
    fig, axs = plt.subplots(1, len(thetas), figsize=(10, 5))

    fig.suptitle(
        r"${}^6\mathrm{Li}$ Compton at $"
        + str(energy)
        + r"\,\mathrm{MeV}$; "
        + r"$\Lambda_{\mathrm{NN}}="
        + f"{lambdaCut}"
        + r"\,\mathrm{MeV}$, "
        + rf"$N_\mathrm{{tot}}={Ntotmax}$",
        fontsize=15,
        y=0.95,
    )

    for k, theta in enumerate(thetas):
        ax = axs[k] if len(thetas) > 1 else axs
        ax.set_title(rf"$\theta={theta}^\circ$", fontsize=14, y=0.02)

        for j, omegaH in enumerate(omegaHs):
            xs, ys = [], []
            for lam in lambdaSRGs:
                val = out.get((lam, omegaH, theta), None)
                if val is not None:
                    xs.append(lam)
                    ys.append(val)

            if len(xs) > 0:
                ax.scatter(
                    xs,
                    ys,
                    marker=markers[j % len(markers)],
                    linestyle="-",
                    linewidth=1.4,
                    c=colors[j % len(colors)],
                    label=rf"$\omega_H={omegaH}$",
                )

        ax.set_xlabel(r"$\Lambda_{\mathrm{SRG}}\; [\mathrm{fm}^{-1}]$")
        ax.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\Omega\;[\mathrm{nb}]$")

        # put legend only on θ = 40° subplot
        if theta == 40:
            ax.legend(loc="lower left", fontsize=11)

    plt.subplots_adjust(wspace=0.30, left=0.10, right=0.95, top=0.90, bottom=0.12)
    plt.show()


if __name__ == "__main__":
    main()
