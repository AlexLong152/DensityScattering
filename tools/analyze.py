# -*- coding: utf-8 -*-

"""
@author: alexl
"""
import numpy as np
# from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join
from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'

mypath = "."
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and
             f[-3:] == "dat"]


def main():
    breakStr = 20*"#"+"\n"
    breakStr = breakStr+breakStr

    for f in onlyfiles:
        print(breakStr)
        out = output(f, [0, 10])
        print(out)


class output:
    """
    A class for reading in data from an output file

    whichnums: ndarray,optional
        which indicies from the array to use, for example in
        pionphotoproduction whichnums=[0,10], remember indicies start at 0

    myfilter: function, optional
        a function that is applied to self.quantNums after whichnums is applied

        self.myfilter=myfilter
        self.quantNums=myfilter(self.quantNums)

        This could for example, just take the real part of the function.
        def filter(nums):
            return np.array([x.real for x in nums])

    """

    def __init__(self, filename, whichnums=None, myfilter=None):
        self.quantNums = getQuantNums(filename)
        # This block gets the indicies that are passed in whichnums
        if not isinstance(whichnums, type(None)):
            # tmp = np.empty(len(whichnums), dtype=np.complex128)
            tmp = []
            for i in range(len(self.quantNums)):
                if i in whichnums:
                    tmp.append(self.quantNums[i])
        self.quantNums = np.array(tmp)
        self.myfilter = myfilter
        if not isinstance(self.myfilter, type(None)):
            self.quantNums = myfilter(self.quantNums)

        self.omega = getOmega(filename)
        self.theta = getOmega(filename)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Called when you print the class"""
        outStr = "omega="+str(self.omega)+"\n"
        outStr += "theta="+str(self.theta)+"\n"
        outStr += "Quantnums:\n----------\n"+str(self.quantNums)
        return outStr


def getOffset(filename):
    part = filename.partition("-72ang")[0]
    part = part.partition("offset=")[2]
    # part=float(part)
    # print(part)
    return float(part)


def getQuantNums(filename):
    """
    Reads in quantum numbers as real numbers from filename
    """
    with open(filename, "r") as f:
        contents = f.read()
        lines = np.array(contents.splitlines())[2:]
        lines = "".join(lines)
        lines = lines.split("=")
        lines = lines[0]
        lines = lines.replace("(", "")
        lines = lines.split(")")
        lines = np.array([x.strip() for x in lines if x.strip() != ''])
        vals = np.zeros(len(lines), dtype=np.complex128)
        for i in range(len(vals)):
            tmp = lines[i].split(",")
            vals[i] = float(tmp[0])+1j*float(tmp[1])
    return vals


def getOmega(filename):
    with open(filename, "r") as f:
        contents = f.read()
        line = np.array(contents.splitlines())[0]
        omega, matches, theta = line.partition("thetacm =   ")
        _, __, omega = omega.partition(" cm omega =  ")
        omega = omega.strip()

    return float(omega)


def getTheta(filename):
    with open(filename, "r") as f:
        contents = f.read()
        line = np.array(contents.splitlines())[0]
        omega, matches, theta = line.partition("thetacm =   ")
        theta = theta.strip()
    return float(theta)


if __name__ == "__main__":
    main()
