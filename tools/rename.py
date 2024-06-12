from os import listdir
from os.path import isfile, join
import os
# adds prefix to files of given filetype
prefix = "varsub-"
filetype = ".f"
onlyfiles = [f for f in listdir(".") if (isfile(join(".", f)) and
                                         f[-len(filetype):] == ".f")]
commands = []
for f in onlyfiles:
    if f[:len(prefix)] != prefix:  # to avoid accidently running this twice
        command = "mv " + f + " varsub-" + f
        commands.append(command)
        print(command)

Q = input("Run these commands? [y/n]: ")
if Q.lower() == "y":
    for c in commands:
        os.system(c)
