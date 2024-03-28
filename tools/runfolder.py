# -*- coding: utf-8 -*-

"""
@author: alexl
"""
import numpy as np
from os import listdir, system
from os.path import isfile, join
from copy import copy

# below is the .pyinput.dat file
"""
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
"""
basefile = ".pyinput.dat"
folder = r"/home/alex/OneDrive/DENSITIES/twobodyvia2Ndensity/XiangXiang-files"

if folder[-1] != r"/":
    folder += r"/"


def main():
    onlyfiles = [f for f in listdir(folder) if isfile(
        join(folder, f)) and (f[-3:] == ".h5" or f[-3:] == "*.gz")]
    # onlyfiles = [onlyfiles[0], onlyfiles[1]]

    runcommand = []
    for i, f in enumerate(onlyfiles):
        path = folder+f  # path to density
        output = folder+"output-"+f
        runfile = basefile[:-4]+"-"+str(i)+".dat"  # input file for fortran
        runcommand.append("./run.twobodyvia2Ndensity "+runfile)

        # print("cp "+basefile + " "+runfile)
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

    runcommand = np.array(runcommand)
    finalcommand = "(trap 'kill 0' SIGINT;make clean; make; SIGINT;"
    finalcommand += runcommand[0]
    for i in range(1, len(runcommand)):
        c = runcommand[i]
        finalcommand += " & "+c
    finalcommand += ")"
    print(finalcommand)
    system(finalcommand)


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


if __name__ == "__main__":
    main()
