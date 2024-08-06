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

basefile = ".pyinput.dat"  # system auto creates this file
folder = r"/home/alex/OneDrive/DENSITIES/twobodyvia2Ndensity/XiangXiang-files/Experimental_6Li_Aug6/N4LO-400MeV-part1/"

if folder[-1] != r"/":
    folder += r"/"


def main():
    remove = input(
        "remove old .pyinput files? will run 'rm .pyinput-*' [y/n]  ")
    if remove.lower() == "y":
        system('rm .pyinput-*')

    remove = input(
        "remove small output files currently in directory? [y/n] ")
    if remove.lower() == "y":
        removeSmallOut()

    writepyInput()
    writeParCommands()

    onlyfiles = [f for f in listdir(folder) if isfile(
        join(folder, f)) and (f[-3:] == ".h5" or f[-3:] == "*.gz")]
    word = "output-"
    onlyfiles = [f for f in onlyfiles if f[:len(word)] != word]

    # print(onlyfiles)
    runcommand = []
    runcommand.append(r"./.pCommand.sh ")

    # removeSmallOut()
    for i, f in enumerate(onlyfiles):

        path = folder+f  # path to density
        output = folder+"output-"+f[:-2]+"dat"
        if not isfile(output):
            runfile = basefile[:-4]+"-"+str(i)+".dat"  # input file for fortran
            runcommand.append(r'"./run.twobodyvia2Ndensity '+runfile+r'"')

            system("cp "+basefile + " "+runfile)
            omega = str(getomega(path))
            theta = str(gettheta(path))

            with open(runfile) as fout:
                s = fout.read()
                # print(s)

            with open(runfile, "w") as fout:
                s = s.replace("XXX", omega)
                s = s.replace("YYY", theta)
                s = s.replace("OUTPUT", output)
                s = s.replace("INPUT", r"'"+path+r"'")
                fout.write(s)

            if len(runcommand)%11==0:
                runcommand.append(r"; ./.pCommand.sh ")

    finalCommand = " ".join(np.array(runcommand))
    print("finalCommand=", finalCommand)

    # runcommand = np.array(runcommand)
    # finalcommand = "(trap 'kill 0' SIGINT;make clean; make; SIGINT;"
    # finalcommand += runcommand[0]
    # for i in range(1, len(runcommand)):
    # c = runcommand[i]
    # finalcommand += " & "+c
    # finalcommand += ")"
    # print(finalcommand)
    system(finalCommand)


def getomega(file):
    # works for the form: omega=5.00E+01-theta=3.00E+01.h5
    file = copy(file)
    contains_omega = file.split("omega=")[1]
    omega = contains_omega.split("-")[0]
    return int(float(omega))


def gettheta(file):
    # works for the form: omega=5.00E+01-theta=3.00E+01.h5
    file = copy(file)
    contains_theta = file.split("theta=")[1]
    theta = contains_theta.split(".h5")[0]
    return int(float(theta))


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
Odelta2_j12max=2 		    Calctype, maximal total ang mom in (12) subsystem
14 2 		    		    NP12A, NP12B
1.1 5.0 15.0 			    P12A, P12B, P12C
2 50     			    AngularType12,(Nordth12 OR Nanggrid12),Nthbins12,Nordphi12,Nphibins12
COMMENTS:_v2.0
# in density filename, code replaces XXX and YYY automatically by energy and angle.
"""[1:]


def writeParCommands():
    filename = ".pCommand.sh"
    system("rm "+filename)
    system("touch "+filename)
    with open(filename, 'w') as file:
        file.write(parallelCommands)

    system("chmod +x "+filename)


def writepyInput():
    system("rm "+basefile)
    system("touch "+basefile)
    with open(basefile, 'w') as file:
        file.write(pyinputText)


# print("pyinputText[0]=", pyinputText[0])
# print("pyinputText=", pyinputText)
# print("parallelCommands[0]=", parallelCommands[0])
# print("parallelCommands=", parallelCommands)

def removeSmallOut():
    """
    Deletes all the small files in the folder, so if theres an issue
    with the code running, you can just use this to auto-delete and run again
    """

    onlyfiles = [f for f in listdir(folder) if isfile(
        join(folder, f)) and f[-3:] == "dat"]

    word = "output-"
    onlyfiles = [f for f in onlyfiles if f[:len(word)] == word]
    for f in onlyfiles:
        size = Path(folder+f).stat().st_size
        if size <= 12000:
            command = "rm "+folder+f
            # print(command)
            system(command)


if __name__ == "__main__":
    main()
    # removeSmallOut()
