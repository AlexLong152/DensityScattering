# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
from os import listdir, system
from os.path import isfile, join
from copy import copy
from pathlib import Path

# below is the .pyinput.dat file
numDiagrams = 1

basefile = ".pyinput.dat"  # system auto creates this file
folder = r"/home/alexander/OneDrive/densities-3He/2Ndensities/132MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/"

if folder[-1] != r"/":
    folder += r"/"

if outputfolder[-1] != r"/":
    outputfolder += r"/"


def main():
    remove = input("remove old .pyinput files? will run 'rm .pyinput-*' [y/n]  ")
    if remove.lower() == "y":
        system("rm .pyinput-*")

    remove = input("remove small output files currently in directory? [y/n] ")
    if remove.lower() == "y":
        removeSmallOut()

    writepyInput()
    writeParCommands()
    # tmp = listdir(folder)
    # f = tmp[0]
    # print("f[-3:]=", f[-3:])
    onlyfiles = [
        f
        for f in listdir(folder)
        if isfile(join(folder, f)) and (f[-3:] == ".h5" or f[-3:] == ".gz")
    ]
    # print(onlyfiles)
    word = "output-"
    onlyfiles = [f for f in onlyfiles if f[: len(word)] != word]

    runcommand = []
    runcommand.append(r"./.pCommand.sh ")
    print("runcommand=", runcommand)
    # removeSmallOut()
    j = 0
    for i, f in enumerate(onlyfiles):
        path = folder + f  # path to density
        output = outputfolder + "output-" + f[:-2] + "dat"
        # outputfolder = folder
        # outputfolder = r"./output-for-dropbox/"

        output = outputfolder + generate2BodOutputName(f, Odelta=2, j12max=2)
        # print("output=", output)
        if not isfile(output):
            j += 1
            runfile = basefile[:-4] + "-" + str(i) + ".dat"  # input file for fortran
            runcommand.append(
                r'"./run.twobodyvia2Ndensity.PionPhotoProdThresh ' + runfile + r'"'
            )

            system("cp " + basefile + " " + runfile)
            omega = str(getomega(f))
            theta = str(gettheta(f))

            with open(runfile) as fout:
                s = fout.read()
                # print(s)

            with open(runfile, "w") as fout:
                s = s.replace("XXX", omega)
                s = s.replace("YYY", theta)
                # s = s.replace("OUTPUT", r"'" + output + r"'")
                s = s.replace("OUTPUT", output)
                s = s.replace("INPUT", r"'" + path + r"'")
                fout.write(s)

            if j % 5 == 0:
                runcommand.append(r"; ./.pCommand.sh ")

    finalCommand = " ".join(np.array(runcommand))
    print("finalCommand=", finalCommand)
    yn = input("Run this command? [y/n]").lower()
    yn = True if yn == "y" else False
    if yn:
        system(finalCommand)


parallelCommands = """
#!/bin/bash

for cmd in "$@"; do {
  echo "Process \"$cmd\" started";
  $cmd & pid=$!
  PID_LIST+=" $pid";
} done

trap "kill $PID_LIST" SIGINT

echo "Parallel processes have started";

wait $PID_LIST

echo
echo "All processes have completed";
"""[1:]

pyinputText = """
XXX XXX 10                   omegaLow, omegaHigh, omegaStep
YYY YYY 15                             thetaLow, thetaHigh, thetaStep
OUTPUT
INPUT
cm_ymmetry_verbos                     frame, symmetry, verbosity of STDOUT
Odelta4_j12max=2_numDiagrams=1 		    Calctype, maximal total ang mom in (12) subsystem, and number of diagrams
14 2 		    		    NP12A, NP12B
1.1 5.0 15.0 			    P12A, P12B, P12C
2 50     			    AngularType12,(Nordth12 OR Nanggrid12),Nthbins12,Nordphi12,Nphibins12
COMMENTS:_v2.0
# in density filename, code replaces XXX and YYY automatically by energy and angle.
"""[1:]


def writeParCommands():
    filename = ".pCommand.sh"
    system("rm " + filename)
    system("touch " + filename)
    with open(filename, "w") as file:
        file.write(parallelCommands)

    system("chmod +x " + filename)


def writepyInput():
    system("rm " + basefile)
    system("touch " + basefile)
    with open(basefile, "w") as file:
        file.write(pyinputText)


def removeSmallOut():
    """
    Deletes all the small files in the folder, so if theres an issue
    with the code running, you can just use this to auto-delete and run again
    """

    onlyfiles = [
        f for f in listdir(folder) if isfile(join(folder, f)) and f[-3:] == "dat"
    ]

    word = "output-"
    onlyfiles = [f for f in onlyfiles if f[: len(word)] == word]
    for f in onlyfiles:
        size = Path(folder + f).stat().st_size
        if size <= 12000:
            command = "rm " + folder + f
            # print(command)
            system(command)


def generate2BodOutputName(densName, Odelta=2, j12max=2):
    # Identify the position of ".denshash="
    insertion_point = densName.find(".denshash=")
    if insertion_point == -1:
        # If ".denshash=" is not found, return the original string or raise an error
        print("densName=", densName)
        raise ValueError("'.denshash=' not found in the input string.")

    # Build the insertion string
    insertion_str = f".Odelta{Odelta}-j12max={j12max}-"

    # Insert the new substring before ".denshash="
    outString = densName[:insertion_point] + insertion_str + densName[insertion_point:]

    last_pos = outString.rfind(".")

    assert last_pos != -1  # char found in output
    # Keep everything up to and including `last_pos` and then add your replacement
    outString = outString[: last_pos + 1] + "dat"
    return outString


def getomega(filename):
    # Find the positions of '.' after the nucleus and 'MeV'
    # Format: ...<NUC>.<ENERGY>MeV-<ANGLE>deg...
    # Step 1: Find the dot after the nucleus name
    dot_index = filename.find(".")
    if dot_index == -1:
        return ""

    # Step 2: Find 'MeV'
    mev_index = filename.find("MeV", dot_index)
    if mev_index == -1:
        return ""

    # Extract the substring between '.' and 'MeV'
    omega_str = filename[dot_index + 1 : mev_index]
    return int(omega_str)


def gettheta(filename):
    # Find 'MeV-' and 'deg'
    # Format: ...<ENERGY>MeV-<ANGLE>deg...
    mev_hyphen_index = filename.find("MeV-")
    if mev_hyphen_index == -1:
        return ""

    deg_index = filename.find("deg", mev_hyphen_index)
    if deg_index == -1:
        return ""

    # Extract substring between 'MeV-' and 'deg'
    theta_str = filename[mev_hyphen_index + 4 : deg_index]  # 'MeV-' is 4 chars long
    return int(theta_str)


if __name__ == "__main__":
    main()
