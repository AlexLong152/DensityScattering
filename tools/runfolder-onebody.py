# -*- coding: utf-8 -*-

"""
@author: alexl
"""

from os import listdir, system
from os.path import isfile, join
import os
from copy import copy
import readDensity as rd

"""
Automatically runs all densities in a folder provided.
TODO: add a taskNumber so that multiple instances can run at the same time.

when running this file automatically creates a bash script to run things in
parrallel, along with input files for the code
"""

basefile = ".pyinput.dat"  # system auto creates this file
folder = (
    r"/home/alexander/OneDrive/densities-6Li/1Ndensities/60MeV/lambda400/lambdaSRG1.88/"
)

if folder[-1] != r"/":
    folder += r"/"

runStandard = True  # runs Odelta3
runVaryA = True  # runs the polarizability calculation

Odelta = "Odelta3"
Odelta = "Odelta2"


def main():
    system(r"rm .pyinput-*")
    writepyInput()
    print("folder=", folder)
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]

    runcommand = []
    out = r"/output/"
    outfolder = folder + out
    # outfolder = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda400/onebody/"
    outfolder = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/newfiles/lambda400/"

    outfolder = outfolder.replace(r"//", r"/")
    outfiles = []

    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    # removeSmallOut()

    outfolder = outfolder.replace(r"//", r"/")
    for i, f in enumerate(onlyfiles):
        path = folder + f  # path to density
        if runStandard:
            output = outfolder + getOutname(f, Odelta)

            if not isfile(output):
                # if True:
                runfile = (
                    basefile[:-4] + "-" + str(i) + ".dat"
                )  # input file for fortran

                command = r"./run.onebodyvia1Ndensity " + runfile + ";"
                runcommand.append(command)

                outfiles.append([output, command])
                system("cp " + basefile + " " + runfile)
                omega = str(getomega(path))
                theta = str(gettheta(path))

                with open(runfile) as fout:
                    s = fout.read()
                    # print(s)

                with open(runfile, "w") as fout:
                    s = s.replace("XXX", omega)
                    s = s.replace("YYY", theta)
                    s = s.replace("OUTPUT", output)
                    s = s.replace("OdeltaVALUE", Odelta)
                    s = s.replace("INPUT", r"'" + path + r"'")
                    fout.write(s)

        if runVaryA:
            nucNumber = getNucNumber(f)
            for val in ["p", "n"]:
                for j in range(1, nucNumber + 1):  # 7 should be numbers of nucelons+1
                    # output = folder + "output-VaryA" + str(j) + val + f[:-3] + "dat"
                    varyAString = "VaryA" + str(j) + val
                    output = outfolder + getOutname(f, varyAString)
                    # print("output=", output)
                    # test = input("blash")
                    if not isfile(output):
                        runfile = (
                            basefile[:-4]
                            + "-"
                            + str(i)
                            + "VaryA"
                            + str(j)
                            + val
                            + ".dat"
                        )  # input file for fortran
                        command = r"./run.onebodyvia1Ndensity " + runfile + ";"
                        runcommand.append(command)
                        outfiles.append([output, command])

                        system("cp " + basefile + " " + runfile)
                        omega = str(getomega(path))
                        theta = str(gettheta(path))

                        with open(runfile) as fout:
                            s = fout.read()
                            # print(s)

                        with open(runfile, "w") as fout:
                            s = s.replace("XXX", omega)
                            s = s.replace("YYY", theta)
                            s = s.replace("OUTPUT", output)
                            s = s.replace("OdeltaVALUE", "VaryA" + str(j) + val)
                            s = s.replace("INPUT", r"'" + path + r"'")
                            fout.write(s)

    #   Its possible for the command to become too long for the program to handle
    #   so break it up into chunks of size step
    print(runcommand)
    step = 50
    for i in range(step, len(runcommand), step):
        finalCommand = " ".join(runcommand[i - step : i])
        system(finalCommand)
    # system(f"cd {folder}")
    # system(
    #     r"mkdir outputs;mv output-* outputs/;tar -zcvf outputs-lambda400-Odelta2.tar.gz outputs/"
    # )

    print("Running these files again")
    for file, command in outfiles:
        if not isfile(file):
            print("The following file was not created")
            print(f"{file}\n")
            os.system(command)


def getOutname(inName, varyAStr):
    """
    varyAstr can be either actually varyA or Odelta3
    """
    tmpName = copy(inName)
    spliter = ".denshash"
    # print("tmpName=", tmpName)
    prefix, suffix = tmpName.split(spliter)
    output = prefix + "." + varyAStr + "-" + spliter + suffix

    last_pos = output.rfind(".")

    assert last_pos != -1  # char found in output
    # Keep everything up to and including `last_pos` and then add your replacement
    output = output[: last_pos + 1] + "dat"
    # print("output=", output)
    # a = input("test")
    # print(a)
    return output


def getomega(filename):
    # Find the positions of '.' after the nucleus and 'MeV'
    # Format: ...<NUC>.<ENERGY>MeV-<ANGLE>deg...
    # Step 1: Find the dot after the nucleus name
    filename = copy(filename)
    filename = filename.split(r"/")[-1]
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


pyinputText = """
XXX XXX 10                           omegaLow, omegaHigh, omegaStep
YYY YYY 15                             thetaLow, thetaHigh, thetaStep
OUTPUT
INPUT
cm_symmetry_verbos                    frame, symmetry, verbosity of STDOUT
OdeltaVALUE                                  Calctype -- VaryAXN varies NucelonN's amplitude AX (N=p,n; X=1-6)
16			             number of points Nx in Feynman parameter integration
COMMENTS:_v2.0
# in output & density filenames, code replaces ORDER, XXX and YYY automatically by order, energy and angle.

"""[1:]


def writepyInput():
    system("rm " + basefile)
    system("touch " + basefile)
    with open(basefile, "w") as file:
        file.write(pyinputText)


def getNucNumber(f):
    tmp = copy(f)
    return int(tmp[8])  # works for nucNumber<10


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


def checkGoodOutput(filename, expectedLength=0):
    try:
        vals = rd.getQuantNums(filename)
    except rd.BadDataError:
        return False

    if expectedLength != 0:
        if len(vals) != expectedLength:
            return False
        else:
            return True
    else:
        if (vals == 0).all():
            return False
        else:
            return True


if __name__ == "__main__":
    main()
