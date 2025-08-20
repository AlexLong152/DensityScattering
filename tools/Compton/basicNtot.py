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
    Ntotmaxplot()


def Ntotmaxplot():
    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/"
    onebody_dir = twobody_dir + r"../onebody/"
    savefolder = twobody_dir + r"../results/"
    tmp = np.array(twobody_dir.split(r"/"))
    title = tmp[-3]

    energy = 60
    lambdaSRG = 1.880
    lambdaCut = 550
    omegaH = 16

    # Ntotmax = 14
    Ntotmaxs = np.array([6, 8, 10, 12, 14])
    thetas = np.array([0, 40, 55, 75, 90, 110, 125, 145, 159, 180])
    markers = ["x", ",", "o", "v", "1", "*", "D"]

    varStr = "Ntotmax"
    varIterable = Ntotmaxs
    total = []
    for i, var in enumerate(varIterable):
        ys = []
        for theta in thetas:
            # need to change variable manually

            ccVal = pdep.ccForDict(
                onebody_dir,
                twobody_dir,
                energy=energy,
                angle=theta,
                lambdaCut=lambdaCut,
                lambdaSRG=lambdaSRG,
                Ntotmax=var,
                omegaH=omegaH,
            )
            ys.append(ccVal)
        total.append(ys)
        labelStr = varStr + "=" + str(var)
        plt.scatter(
            thetas + Ntotmaxs[i] - np.min(Ntotmaxs),
            ys,
            label=labelStr,
            marker=markers[i],
        )
    plt.legend()
    title = "6Li Compton Scattering " + title
    title += f"\nomegaH={omegaH}, lambdaSRG={lambdaSRG}"
    plt.title(title)
    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\;\mathrm{ sr^{-1}} ]$")
    plt.show()

    array_out = pdep.array_to_table(
        np.array(total).T,
        [str(theta) for theta in thetas],
        [str(x) for x in varIterable],
    )
    out = title
    # out += f"\nNtotmax={Ntotmax}, omegaH={omegaH}"
    out += f"\nChange of {varStr} in the columns, theta values by rows\n" + str(
        array_out
    )
    print(out)
    fileName = title.replace(" ", "-")
    fileName = fileName.replace(",", "")
    fileName = fileName.replace("Scattering-", "")
    fileName = fileName + "omegaH=" + str(omegaH)
    fileName = varStr + "_vary_" + fileName
    print("fileName=", fileName)
    pdep.save_string_to_file(savefolder + r"/" + fileName + ".txt", out)


if __name__ == "__main__":
    main()
