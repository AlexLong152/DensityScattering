from os import listdir
from os.path import isfile, join
import os
from copy import copy

folder = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda400/"
"""
rewrite a new "newname" function for every use case
in all likelyhood you'll have to run this a few times unless youre a god and can do it perfectly 
every time
"""

if folder[-1] != r"/":
    folder += r"/"


def newname(f):
    # rewrite me every time
    f = f.replace("onebodyonebody", "onebody")
    return f


# def replaceString(f):
#     return f.replace("output", "6Li")


def endsIn(string, substring):
    return string.endswith(substring)


def skipFile(f):
    # if skipFile returns true then the file isn't renamed
    # return False
    return endsIn(f, ".gz")


def main():
    onlyfiles = [f for f in listdir(folder) if (isfile(join(folder, f)))]
    moveFolder = input(
        "Move renamed filed into a new folder when done? [y/n]: "
    ).lower()
    moveFolder = True if moveFolder == "y" else False

    if moveFolder:
        newfolder = input("Name new folder without slashes: ")
        newfolder = r"/" + newfolder + r"/"

    else:
        newfolder = ""

    commands = []
    for f in onlyfiles:
        if not skipFile(f):
            sourceFile = folder + r"/" + f
            # newfile = (
            #     folder + newname(f)
            #     if moveFolder is False
            #     else
            # )
            newfile = folder + newfolder + newname(f)
            command = f"mv {sourceFile} {newfile}"
            command = command.replace(r"//", r"/")
            print(command)
            tmp = copy(command)
            for x in tmp.split():
                print(x)
            print("")
            commands.append(command)

    Q = input("Run these commands? [y/n]: ")
    if Q.lower() == "y":
        os.system("mkdir " + folder + newfolder)
        for c in commands:
            os.system(str(c))


if __name__ == "__main__":
    main()
