# -*- coding: utf-8 -*-
"""Compare cross sections across regulator choices for select angles."""

import os
import numpy as np
import matplotlib.pyplot as plt
import CrossSection as cc

from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"

energies = [60, 75, 86, 100]
Odeltaonebod = "Odelta3"
Odeltatwobod = "Odelta4"
Ntotmax = 14

lambdaSRGs = np.array([1.880, 2.236])
lambdaCuts = np.array([450, 500])
omegaHs = np.array([16, 18])
thetas = np.array([40, 159], dtype=float)

saveFolder = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/results/FinalCrossSections-And-Uncertainties/"
plot_path = os.path.join(saveFolder, "Compton-6Li-cross-section-scatter.pdf")
writeOutput = True


def main():
    missing_all = []
    plot_data = {}

    for energy in energies:
        missing, energy_data = run_for_energy(energy)
        missing_all.extend(missing)
        plot_data[energy] = energy_data

    print_missing_summary(missing_all)
    make_plots(plot_data)


def run_for_energy(energy):
    outfile = os.path.join(
        saveFolder, f"Compton-6Li-{energy}MeV-distribution-output.txt"
    )
    if writeOutput and os.path.exists(outfile):
        os.remove(outfile)

    writer = make_writer(outfile)
    onebody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/1bod/{energy}MeV/"
    twobody_dir = f"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/2bod/{energy}MeV/"
    theta_values = [float(theta) for theta in thetas]

    data, missing = gather_cross_sections(
        onebody_dir, twobody_dir, theta_values, energy
    )

    header_line = (
        "lambdaSRG values: "
        + format_values(lambdaSRGs)
        + "; lambdaCut values: "
        + format_values(lambdaCuts)
        + "; omegaH values: "
        + format_values(omegaHs)
    )

    writer(f"{energy} MeV")
    writer(header_line)

    for theta in theta_values:
        writer("")
        writer(f"theta = {theta:.0f} degrees")
        writer(
            f"{'lambdaSRG':>10} {'lambdaCut':>10} {'omegaH':>8} "
            f"{'1N Cross Sec [nb]':>20} {'1N+2N Cross Sec [nb]':>24}"
        )
        for entry in data[theta]:
            writer(
                f"{entry['lambdaSRG']:10.3f} {entry['lambdaCut']:10.0f} "
                f"{entry['omegaH']:8.0f} {entry['cross_section_1N']:20.4f} "
                f"{entry['cross_section_total']:24.4f}"
            )

    print("outfile= ", outfile)
    return missing, data


def gather_cross_sections(onebody_dir, twobody_dir, theta_values, energy):
    data = {theta: [] for theta in theta_values}
    missing = []

    for theta in theta_values:
        for lambdaSRG in lambdaSRGs:
            for lambdaCut in lambdaCuts:
                for omegaH in omegaHs:
                    kwargs = dict(
                        energy=energy,
                        angle=theta,
                        lambdaCut=lambdaCut,
                        lambdaSRG=lambdaSRG,
                        Ntotmax=Ntotmax,
                        omegaH=omegaH,
                    )
                    try:
                        cross_section_1N = cc.ccForDict(
                            onebody_dir,
                            twobody_dir,
                            Odeltaonebod=Odeltaonebod,
                            Odeltatwobod=Odeltatwobod,
                            zeroTwobod=True,
                            **kwargs,
                        )
                        cross_section_total = cc.ccForDict(
                            onebody_dir,
                            twobody_dir,
                            Odeltaonebod=Odeltaonebod,
                            Odeltatwobod=Odeltatwobod,
                            zeroTwobod=False,
                            **kwargs,
                        )
                    except FileNotFoundError:
                        missing.append(
                            {
                                "energy": energy,
                                "theta": theta,
                                "lambdaSRG": lambdaSRG,
                                "lambdaCut": lambdaCut,
                                "omegaH": omegaH,
                            }
                        )
                        continue

                    data[theta].append(
                        {
                            "lambdaSRG": lambdaSRG,
                            "lambdaCut": lambdaCut,
                            "omegaH": omegaH,
                            "cross_section_1N": cross_section_1N,
                            "cross_section_total": cross_section_total,
                        }
                    )

    return data, missing


def make_plots(plot_data):
    fig, axes = plt.subplots(1, 4, figsize=(12, 8))
    axes = axes.flatten()
    colors = {40.0: "tab:blue", 159.0: "tab:orange"}

    for idx, energy in enumerate(energies):
        ax = axes[idx]
        energy_data = plot_data.get(energy, {})
        had_points = False
        for theta in sorted(energy_data.keys()):
            entries = energy_data[theta]
            if not entries:
                continue
            had_points = True
            xs = np.full(len(entries), theta, dtype=float)
            ys = [entry["cross_section_total"] for entry in entries]
            if len(entries) > 1:
                jitter = np.linspace(-0.6, 0.6, num=len(entries))
            else:
                jitter = np.zeros(1)
            label = rf"$\theta = {theta:.0f}^\circ$"
            existing_labels = set(ax.get_legend_handles_labels()[1])
            if label in existing_labels:
                label = "_nolegend_"
            ax.scatter(
                xs + jitter, ys, s=45, label=label, color=colors.get(theta, "tab:gray")
            )
        ax.set_title(rf"$E_\gamma = {energy}\,\mathrm{{MeV}}$")
        ax.set_xlabel(r"$\theta\,(^{\circ})$")
        ax.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\Omega[\mathrm{nb}]$")
        if had_points:
            ax.legend(prop={"size": 9})

    fig.suptitle(r"$^6$Li Compton Scattering Cross Sections", fontsize=16)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.show()
    plt.close(fig)
    print("Scatter plot saved to:", plot_path)


def format_values(values):
    return "[" + ", ".join(f"{val:g}" for val in values) + "]"


def make_writer(outfile):
    def _writer(*args, sep=" ", end="\n", file=None, flush=False):
        print(*args, sep=sep, end=end, file=file, flush=flush)
        if writeOutput:
            with open(outfile, "a", encoding="utf-8") as f:
                print(*args, sep=sep, end=end, file=f)

    return _writer


def print_missing_summary(missing):
    print("Missing combinations:")
    print(missing)


if __name__ == "__main__":
    main()
