# -*- coding: utf-8 -*-

"""
@author: alexl
"""
import numpy as np
from os import listdir, system
from os.path import isfile, join
from copy import copy
from pathlib import Path
import time
"""
Automatically runs all densities in a folder provided.
TODO: add a taskNumber so that multiple instances can run at the same time.

when running this file automatically creates a bash script to run things in
parrallel, along with input files for the code
"""

basefile = ".pyinput.dat"  # system auto creates this file
folder=r'/home/alexander/OneDrive/DENSITIES/onebodyvia1Ndensity/chiralsmsN4LO+3nfN2LO-lambda400-SRG/onebody/'


if folder[-1] != r"/":
    folder += r"/"

runStandard=True
varyA=True

def main():


    system(r'rm .pyinput-*')
    writepyInput()

    onlyfiles = [f for f in listdir(folder) if isfile(
        join(folder, f)) and (f[-3:] == "dat") ]

#   print("onlyfiles=",onlyfiles)


    word = "output-"
    onlyfiles = [f for f in onlyfiles if f[:len(word)] != word]
    runcommand = []

    # removeSmallOut()
    for i, f in enumerate(onlyfiles):

        path = folder+f  # path to density
        output = folder+"output-ORDER-"+f[:-3]+"dat"
        if not isfile(output):
            if runStandard:
                runfile = basefile[:-4]+"-"+str(i)+".dat"  # input file for fortran

                runcommand.append(r'./run.onebodyvia1Ndensity '+runfile+';')

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

            if varyA:
                for val in ["p","n"]:
                    for j in range(1,7):
                        output = folder+"output-VaryA"+str(j)+val+f[:-3]+"dat"

                        runfile = basefile[:-4]+"-"+str(i)+"VaryA"+str(j)+val+".dat"  # input file for fortran

                        runcommand.append(r'./run.onebodyvia1Ndensity '+runfile+';')

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
                            s = s.replace("Odelta3", "VaryA"+str(j)+val)
                            s = s.replace("INPUT", r"'"+path+r"'")
                            fout.write(s)

#   Its possible for the command to become too long for the program to handle
#   so break it up into chuncks of size step
    runcommand=np.array(runcommand)
    step=50
    for i in range(step,len(runcommand), step):
        finalCommand = " ".join(runcommand[i-step:i])
        system(finalCommand)
    system(f"cd {folder}")
    system(r"mkdir outputs;mv output-* outputs/;tar -zcvf outputs-lambda400-Odelta2.tar.gz outputs/")



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
    theta = contains_theta.split(".dat")[0]
    return int(float(theta))


pyinputText = """
XXX XXX 10                           omegaLow, omegaHigh, omegaStep
YYY YYY 15                             thetaLow, thetaHigh, thetaStep
OUTPUT
INPUT
cm_symmetry_verbos                    frame, symmetry, verbosity of STDOUT
Odelta2				                                    Calctype -- VaryAXN varies NucelonN's amplitude AX (N=p,n; X=1-6)
16			             number of points Nx in Feynman parameter integration
COMMENTS:_v2.0
# in output & density filenames, code replaces ORDER, XXX and YYY automatically by order, energy and angle.

"""[1:]


def writepyInput():
    system("rm "+basefile)
    system("touch "+basefile)
    with open(basefile, 'w') as file:
        file.write(pyinputText)

# def removeSmallOut():
#   """
#   Deletes all the small files in the folder, so if theres an issue
#   with the code running, you can just use this to auto-delete and run again
#   """

#   onlyfiles = [f for f in listdir(folder) if isfile(
#       join(folder, f)) and f[-3:] == "dat"]

#   word = "output-"
#   onlyfiles = [f for f in onlyfiles if f[:len(word)] == word]
#   for f in onlyfiles:
#       size = Path(folder+f).stat().st_size
#       if size <= 12000:
#           command = "rm "+folder+f
#           # print(command)
#           system(command)


if __name__ == "__main__":
    main()
