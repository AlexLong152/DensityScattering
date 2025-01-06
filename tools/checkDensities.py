import readDensity as rd
import numpy as np
from os.path import isfile, join
from os import listdir

omegaVals = np.array([60])
angles = np.array([40, 55, 75, 90, 110, 125, 145, 159])
Ntotmaxs = np.array([6, 8, 10, 12, 14])
omegaHs = np.array([10, 12, 14, 16, 18, 20, 22, 24])

densFolder = (
    r"/home/alexander/OneDrive/densities-6Li/1Ndensities/60MeV/lambda400/lambdaSRG1.88/"
)

onlyfiles = [
    f for f in listdir(densFolder) if isfile(join(densFolder, f)) and f[-3:] == ".h5"
]


def main():
    densDict = {}

    for f in onlyfiles:
        f = densFolder + f
        omegaVal = rd.getOmega(f)
        thetaVal = rd.getTheta(f)
        NtotmaxVal = rd.getNtotmax(f)
        omegaHVal = rd.getOmegaH(f)
        densDict[(omegaVal, thetaVal, NtotmaxVal, omegaHVal)] = f

    badKeys = []
    for omega in omegaVals:
        for theta in angles:
            for Nto in Ntotmaxs:
                for omegaH in omegaHs:
                    key = (omega, theta, Nto, omegaH)
                    value = densDict.get(key, False)
                    if isinstance(value, str):
                        value = True
                    if not value:
                        badKeys.append(key)
    output(badKeys)

    print(checkIfExists(densDict, 60, 159, 10, 20))


def output(badKeys):
    print(f"Input folder is:\n{densFolder}")
    num = len(badKeys)
    print(f"There are {num} densities missing")
    print(
        "Missing densities are listed in the following order:\n omega, theta,Ntotmax,omegaH"
    )
    for key in badKeys:
        print(tuple(int(x) for x in key))


def checkIfExists(densDict, omega, theta, Nto, omegaH):
    key = (omega, theta, Nto, omegaH)
    return densDict.get(key, False)


if __name__ == "__main__":
    main()
