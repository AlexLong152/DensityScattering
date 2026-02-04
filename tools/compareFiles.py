import os
import numpy as np
from nucdens import access
import pandas as pd
from readDensity import getHash

from os import listdir
from os.path import isfile, join

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_colwidth", None)
pd.set_option(
    "display.width", None
)  # Auto-detect terminal width, or use a number like 200, 300, etc.

colsel = [
    "E[MeV]",
    "addtime",
    "kind",
    "N",
    "Z",
    "Jtot",
    "E[MeV]",
    "omega",
    "theta",
    "orderNN",
    "LambdaNN",
    "tnforder",
    "hashname",
]


def getPaths(folder):
    if folder[-1] != r"/":
        folder += r"/"

    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    paths = []
    for f in files:
        paths.append(folder + f)
    return paths


def main(df):
    # hash1 = "f7af76348422787e794c858a6ce3e3b46ff04d8bd0d837cc7a86739a9a038efc"
    # hash2 = "778ff350a0f12ba9915d142b5b543534ed59e9a9bb367d55fecee55c70f1b45d"
    # hash3 = "9956771e795fda5a40e93012e7f087d89c7aeb88adb330ce8d179a36bf673661"

    folder = "/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/"
    folder = r"/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/"
    folder = r"/home/alexander/OneDrive/densities-3H/1Ndensities/132MeV/60deg/"
    files = getPaths(folder)

    hashes = [getHash(f) for f in files]
    select = df["hashname"].isin(hashes)
    print(df[select][colsel].to_string(col_space=2))


beta_webbase = "https://just-object.fz-juelich.de:9000/jikp03/densitystore-beta/"
workdir = os.environ["HOME"] + r"/OneDrive/"
densdb = access.database(workdir=workdir, webbase=beta_webbase)

workdir = os.environ["HOME"] + r"/OneDrive/"
beta_webbase = "https://just-object.fz-juelich.de:9000/jikp03/densitystore-beta/"

file = f"{workdir}densities_table.h5"
try:
    os.remove(file)
except FileNotFoundError:
    pass
densdf = access.database(workdir=workdir, webbase=beta_webbase)
df = densdf.pddf

if __name__ == "__main__":
    print("\n")
    main(df)
