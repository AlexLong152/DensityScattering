# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import CrossSection as cc
import sys
import os
from os import listdir
from os.path import isfile, join

sys.path.insert(1, "..")

# Odeltaonebod = "Odelta3"
# Odeltatwobod = "Odelta4"

outFiles = {}


def main(Odeltaonebod, Odeltatwobod, energy=86):
    base = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/"

    lambdaSRG = 1.880
    Ntotmax = 14
    lambdaCuts = [450, 500]
    lambdaSRGs = [1.880, 2.236]
    omegaH = 18
    thetas = np.array([1, 40, 55, 75, 90, 110, 125, 145, 159, 180])

    for alphaDiff in [-2, 0, 2]:
        for betaDiff in [-2, 0, 2]:
            outfile = f"polarVary/6Li-{energy}MeV-Compton-onebod={Odeltaonebod}-twobod={Odeltatwobod}-alpha+={alphaDiff}-beta+={betaDiff}.txt"

            if os.path.exists(outfile):
                os.remove(outfile)

            def wprint(s):
                # Print to screen
                print(s)
                # Open in append mode: creates the file if it doesn't exist
                with open(outfile, "a") as f:  # needs outfile from global context
                    f.write(s + "\n")

            # angle = 180
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
            for x in thetas:
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
                            "Odeltaonebod": Odeltaonebod,
                            "Odeltatwobod": Odeltatwobod,
                        }
                        onebod, twobod, kwargs = cc.getMatchingFiles(
                            onebody_dir, twobody_dir, **kwargsfind
                        )
                        assert Odeltaonebod in onebod
                        assert Odeltatwobod in twobod
                        yVal = cc.crossSection(
                            onebody_dir + onebod,
                            twobody_dir + twobod,
                            deltaAlpha=alphaDiff,
                            deltaBeta=betaDiff,
                        )["cc"]
                        ccs.append(yVal)
                        outDict[(lambdaSRG, lambdaCut)].append([float(x), yVal])
                        outFiles[
                            (x, lambdaSRG, lambdaCut, Odeltaonebod, Odeltatwobod)
                        ] = (onebod, twobod)
                        i += 1

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
            )

            data2014 = np.array(
                [
                    [40.0, 203.0, 2.0 * np.sqrt(106)],
                    [55.0, 147.0, np.sqrt(170)],
                    [75.0, 140.0, np.sqrt(193)],
                    [90.0, 159.0, np.sqrt(233)],
                    [110.0, 146.0, np.sqrt(449)],
                    [145.0, 167.0, 2.0 * np.sqrt(97)],
                    [159.0, 172.0, 15.0],
                ],
                dtype=float,
            )
            match energy:
                case 60:
                    data = data2012
                case 86:
                    data = data2014
                case _:
                    raise ValueError("Bad energy given")

            for lambdaSRG in lambdaSRGs:
                for lambdaCut in lambdaCuts:
                    wprint(
                        "lambdaSRG= "
                        + str(lambdaSRG)
                        + " lambdaCut= "
                        + str(lambdaCut)
                        + " onebody="
                        + str(Odeltaonebod)
                        + " twobody="
                        + str(Odeltatwobod)
                    )
                    wprint(40 * "-")
                    outList = outDict[(lambdaSRG, lambdaCut)]
                    wprint(f"{'theta':>7} | {'Cross Section':>14}")
                    for val in outList:
                        x, y = val
                        # wprint(f"{x:7.2f} , {y:14.6f}")
                        wprint(f"{x:7.2f} , {y:16.12f}")
            wprint("\n")
            wprint("Jerry experimental Data")
            wprint(f"{'theta':>7} | {'val':>10} | {'uncertainty':>12}")
            wprint(40 * "-")
            for theta, val, uncer in data:
                wprint(f"{theta:7.2f} | {val:10.6f} | {uncer:12.6f}")
            wprint("AlphaOffset= " + str(alphaDiff))
            wprint("BetaOffset= " + str(betaDiff))


def getData(file):
    """
    Gets lines 3-13, 16-26, 29-39 and 42-52 from a file
    """
    data = []
    with open(file, "r") as f:
        lines = f.readlines()
        for i in range(3, 13):
            data.append(lines[i].strip())
        for i in range(16, 26):
            data.append(lines[i].strip())
        for i in range(29, 39):
            data.append(lines[i].strip())
        for i in range(42, 52):
            data.append(lines[i].strip())
    return data


if __name__ == "__main__":
    # energy = 86
    # OdeltaOneBod = "Odelta2"
    ones = ["Odelta0", "Odelta3", "Odelta3"]
    twos = ["Odelta0", "Odelta2", "Odelta4"]
    for energy in [60, 86]:
        for one, two in zip(ones, twos):
            main(one, two, energy)

    # Rest of this code checks for duplicates in the output files
    folder = r"./polarVary/"
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    paths = []
    for f in files:
        paths.append(folder + f)

    # path1 and path2 are output files, loop over all the output files
    for path1 in paths:
        for path2 in paths:
            if path1 != path2:
                data = getData(path1)
                data2 = getData(path2)
                comp = [a == b for a, b in zip(data, data2)]
                if sum(comp) > 7:
                    if "Odelta4" in path1 and "Odelta2" in path2:
                        # WLOG pick these settings
                        theta = 1
                        lambdaSRG = 450
                        lambdaCut = 1.880
                        files1 = outFiles[
                            (theta, lambdaSRG, lambdaCut, "Odelta3", "Odelta2")
                        ]
                        files2 = outFiles[
                            (theta, lambdaSRG, lambdaCut, "Odelta3", "Odelta4")
                        ]
                        for i, f in enumerate([files1, files2]):
                            onebod, twobod = f
                            print(f"file {i + 1}:")
                            print("onebod=\n", onebod)
                            print("twobod=\n", twobod)
                        print(50 * "-")
                        # print("files1=", files1)
                        # print("files2=", files2)
                        # print("data=", data)
                        # print("data2=", data2)
                        # print("path1=", path1)
                        # print("path2=", path2)
                        print(f"kdiff3 {path1} {path2}")
                        print(50 * "-")
                        print("\n")
                        print("\n")
                        assert False
