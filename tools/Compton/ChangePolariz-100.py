# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import CrossSection as cc
import sys
import os

sys.path.insert(1, "..")


def main():
    base = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/"
    energy = 100
    for alphaDiff in [-2, 0, 2]:
        for betaDiff in [-2, 0, 2]:
            outfile = f"polarVary/6Li-{energy}MeV-Compton-alpha+={alphaDiff}-beta+={betaDiff}.txt"

            if os.path.exists(outfile):
                os.remove(outfile)

            def wprint(s):
                # Print to screen
                print(s)
                # Open in append mode: creates the file if it doesn't exist
                with open(outfile, "a") as f:  # needs outfile from global context
                    f.write(s + "\n")

            # angle = 180
            lambdaSRG = 1.880
            Ntotmax = 14
            lambdaCuts = [450, 500]
            lambdaSRGs = [1.880, 2.236]
            omegaH = 18

            thetas = np.array([1, 40, 55, 75, 90, 110, 125, 145, 159, 180])
            thetas = np.array([40, 159])
            xs = thetas
            ccs = []

            onebody_dir = base + r"1bod/ENERGY!!!MeV/"
            twobody_dir = base + r"2bod/ENERGY!!!MeV/"
            onebody_dir = onebody_dir.replace("ENERGY!!!", str(energy))
            twobody_dir = twobody_dir.replace("ENERGY!!!", str(energy))
            # print("onebody_dir=", onebody_dir)
            # print("twobody_dir=", twobody_dir)
            outDict = {}

            for lambdaSRG in lambdaSRGs:
                for lambdaCut in lambdaCuts:
                    outDict[(lambdaSRG, lambdaCut)] = []
            for x in xs:
                # kwargsfind[plotting] = x

                # print("kwargsfind=\n", kwargsfind)
                i = 0
                for lambdaSRG in lambdaSRGs:
                    for lambdaCut in lambdaCuts:
                        kwargsfind = {
                            "energy": energy,
                            "angle": x,
                            "lambdaCut": lambdaCut,
                            "lambdaSRG": lambdaSRG,
                            "Ntotmax": Ntotmax,
                            "omegaH": omegaH,
                        }
                        onebod, twobod, kwargs = cc.getMatchingFiles(
                            onebody_dir, twobody_dir, **kwargsfind
                        )

                        if onebod is not None and twobod is not None:
                            yVal = cc.crossSection(
                                onebody_dir + onebod,
                                twobody_dir + twobod,
                                deltaAlpha=alphaDiff,
                                deltaBeta=betaDiff,
                            )["cc"]
                            ccs.append(yVal)
                            outDict[(lambdaSRG, lambdaCut)].append([float(x), yVal])
                        i += 1

            for lambdaSRG in lambdaSRGs:
                for lambdaCut in lambdaCuts:
                    wprint(
                        "lambdaSRG= " + str(lambdaSRG) + " lambdaCut= " + str(lambdaCut)
                    )
                    wprint(40 * "-")
                    outList = outDict[(lambdaSRG, lambdaCut)]
                    wprint(f"{'theta':>7} | {'Cross Section':>14}")
                    for val in outList:
                        x, y = val
                        wprint(f"{x:7.2f} , {y:14.6f}")
            wprint("\n")
            wprint(40 * "-")
            wprint("AlphaOffset= " + str(alphaDiff))
            wprint("BetaOffset= " + str(betaDiff))
            wprint(40 * "-")
    wprint("cc.alphap=" + str(cc.alphap))
    wprint("cc.alphan=" + str(cc.alphan))
    wprint("cc.betap=" + str(cc.betap))
    wprint("cc.betan=" + str(cc.betan))


if __name__ == "__main__":
    main()
