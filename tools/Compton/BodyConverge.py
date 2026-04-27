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
rcParams["axes.labelsize"] = 12
rcParams["xtick.labelsize"] = 10
rcParams["ytick.labelsize"] = 10
rcParams["legend.fontsize"] = 10
energy = 86
Odeltaonebod = "Odelta3"
Odeltatwobod = "Odelta4"
error = 0.03


Ntot = 14
base = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/"
lambdaSRGs = np.array([1.880, 2.236])
lambdaCuts = np.array([450, 500])
omegaHs = np.array([16, 18])


def main():
    energies = [60, 75, 86, 100]
    styles = ["-", "--", "-.", ":"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    fig, ax = plt.subplots(figsize=(5.5, 3.7))

    print(f"Odelta: 1bod={Odeltaonebod}, 2bod={Odeltatwobod}")
    print(f"lambdaSRGs={lambdaSRGs}, lambdaCuts={lambdaCuts}, omegaHs={omegaHs}")

    for eng, style, color in zip(energies, styles, colors):
        onebody_dir = base + f"1bod/{eng}MeV/"
        thetas = np.sort(np.array(cc.getValuesAvailable(onebody_dir, "angle")))
        thetas = thetas[thetas != 0]

        print(f"\nEnergy = {eng} MeV, Ntotmax={Ntot}")
        print(f"{'theta':>6}  {'ratio-1':>10} {'± unc':>8}")
        print("-" * 28)

        plot_thetas = []
        plot_vals = []
        plot_uncs = []
        for theta in thetas:
            ratio, ratio_unc = bodyRatio(theta, eng)
            if ratio is None:
                print(f"{theta:6.1f}  {'missing':>10}")
            else:
                print(f"{theta:6.1f}  {ratio - 1:10.4f} {ratio_unc:8.4f}")
                plot_thetas.append(theta)
                plot_vals.append(ratio - 1)
                plot_uncs.append(ratio_unc)

        plot_thetas = np.array(plot_thetas)
        plot_vals = np.array(plot_vals)
        plot_uncs = np.array(plot_uncs)
        ax.plot(
            plot_thetas,
            plot_vals,
            ls=style,
            color=color,
            lw=1.5,
            label=rf"$\omega = {eng}$~MeV, $Q_{{\mathrm{{few}}}}={np.mean(plot_vals):.2f}$",
        )
        ax.fill_between(
            plot_thetas,
            plot_vals - plot_uncs,
            plot_vals + plot_uncs,
            color=color,
            alpha=0.15,
        )

    ax.set_xlabel(r"$\theta_\mathrm{cm}$ [deg]")
    ax.set_ylabel(
        r"$\sqrt{(\mathrm{d}\sigma/\mathrm{d}\Omega)^{1\mathrm{N}+2\mathrm{N}}"
        r"\,/\,(\mathrm{d}\sigma/\mathrm{d}\Omega)^{1\mathrm{N}}} - 1$"
    )
    ax.set_xticks(np.arange(0, 181, 30))
    ax.set_xlim([-5, 185])
    ax.set_ylim([0, np.max(plot_vals) * 1.1])
    ax.legend(framealpha=0.9, edgecolor="gray")
    fig.tight_layout()
    fig.savefig(
        "/home/alexander/OneDrive/main/compton/twobody_convergence_6Li.pdf", dpi=600
    )
    plt.show()


def bodyRatio(theta, omega):
    """
    Returns (ratio, ratio_unc) of sqrt(cross section with 1+2 body) / sqrt(cross section with 1 body only),
    averaged over lambdaSRG, lambdaCut, omegaH.
    """
    onebody_dir = base + f"1bod/{omega}MeV/"
    twobody_dir = base + f"2bod/{omega}MeV/"
    vals_one = []
    vals_full = []
    for lambdaSRG in lambdaSRGs:
        for lambdaCut in lambdaCuts:
            for omegaH in omegaHs:
                try:
                    cc_one = cc.ccForDict(
                        onebody_dir,
                        twobody_dir,
                        Odeltaonebod=Odeltaonebod,
                        Odeltatwobod=Odeltatwobod,
                        zeroTwobod=True,
                        energy=omega,
                        angle=theta,
                        lambdaCut=lambdaCut,
                        lambdaSRG=lambdaSRG,
                        Ntotmax=Ntot,
                        omegaH=omegaH,
                    )
                    cc_full = cc.ccForDict(
                        onebody_dir,
                        twobody_dir,
                        Odeltaonebod=Odeltaonebod,
                        Odeltatwobod=Odeltatwobod,
                        zeroTwobod=False,
                        energy=omega,
                        angle=theta,
                        lambdaCut=lambdaCut,
                        lambdaSRG=lambdaSRG,
                        Ntotmax=Ntot,
                        omegaH=omegaH,
                    )
                    vals_one.append(cc_one)
                    vals_full.append(cc_full)
                except FileNotFoundError:
                    pass

    if len(vals_one) == 0:
        return None, None

    sqrt_one = np.sqrt(np.mean(vals_one))
    sqrt_full = np.sqrt(np.mean(vals_full))
    unc_one = maxDev(np.sqrt(np.array(vals_one)))
    unc_full = maxDev(np.sqrt(np.array(vals_full)))

    return ratioUnc(sqrt_one, sqrt_full, unc_one, unc_full)


def maxDev(ar):
    tmp = (np.max(ar) - np.min(ar)) / 2
    mean = np.mean(ar)
    return np.sqrt((mean * error) ** 2 + tmp**2)


def ratioUnc(sqrt_one, sqrt_full, unc_one, unc_full):
    full = sqrt_full
    one = sqrt_one
    ratio = full / one
    # sigma_z^2 = sigma_x^2 df/dx^2 + sigma_y^2 df/dy^2, where z = f(x,y)=x/y
    uncerSqr = (unc_full / one) ** 2 + (full * unc_one / one**2) ** 2
    return ratio, np.sqrt(uncerSqr)


if __name__ == "__main__":
    main()
