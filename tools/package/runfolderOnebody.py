# -*- coding: utf-8 -*-

"""
@author: alexl
"""

from os import listdir, system
from os.path import isfile, join
import os
from copy import copy
import readDensity as rd
import runfolderTwobody as rtb

"""
Automatically runs all densities in a folder provided.
when running this file automatically creates a bash script to run things in
parrallel, along with input files for the code

Really there should be three seperate values for running the 
varyA, Odelta3, and Odelta2 scripts but theres not need to fix it for now 
"""


basefile = ".pyinput.dat"  # system auto creates this file
Odelta = -1


def main(folder, outfolder):
    if folder[-1] != r"/":
        folder += r"/"

    if outfolder[-1] != r"/":
        outfolder += r"/"
    assert Odelta > -1

    OdeltaStr = f"Odelta{Odelta}"
    system(r"rm .pyinput-*")
    writepyInput()
    print("folder=", folder)
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]

    runcommand = []

    outfiles = []

    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    # removeSmallOut()

    # outfolder = outfolder.replace(r"//", r"/")
    for i, f in enumerate(onlyfiles):
        path = folder + f  # path to density
        output = outfolder + getOutname(f, OdeltaStr)
        # if not isfile(output):
        runfile = basefile[:-4] + "-" + str(i) + ".dat"  # input file for fortran

        command = r"./run.onebodyvia1Ndensity " + runfile + ";"
        runcommand.append(command)

        outfiles.append([output, command])
        system("cp " + basefile + " " + runfile)
        theta = str(gettheta(path))

        with open(runfile) as fout:
            s = fout.read()
            # print(s)
        omega = rtb.getomega(f)
        with open(runfile, "w") as fout:
            s = s.replace("XXX", str(omega))
            s = s.replace("YYY", theta)
            s = s.replace("OUTPUT", output)
            s = s.replace("OdeltaVALUE", OdeltaStr)
            s = s.replace("INPUT", r"'" + path + r"'")
            fout.write(s)

    step = 50
    commandTmp = copy(runcommand)
    while True:
        runcommand = commandTmp[:step]
        del commandTmp[:step]
        finalCommand = " ".join(runcommand)
        system(finalCommand)
        if len(commandTmp) == 0:
            break

    badfiles = [file for file in outfiles if not (isfile)]
    while len(badfiles) > 0:
        print("Running these files again")
        print(badfiles)
        for file, command in outfiles:
            if not isfile(file):
                # print("The following file was not created")
                # print(f"{file}\n")
                os.system(command)

        badfiles = [file for file in outfiles if not (isfile)]
        again = input("Run inputs again? [y/n] ")
        again = again.lower() == "y"
        if not again:
            break


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
cm_symmetry_verbos_extQnumlimit=3                    frame, symmetry, verbosity of STDOUT
OdeltaVALUE                                  Calctype -- VaryAXN varies NucelonN's amplitude AX (N=p,n; X=1-6)
32			             number of points Nx in Feynman parameter integration
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
