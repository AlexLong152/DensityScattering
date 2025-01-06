import os

import pandas as pd
from nucdens import access
import numpy as np

# from datetime import datetime as dt
from requests.exceptions import HTTPError

pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)
energyEps = 0.001
workdir = os.environ["HOME"] + r"/OneDrive/"
beta_webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore-beta/"
webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore/"

filename = "_KIND_body-_NUC_._ENERGY_MeV-_ANGLE_deg.dens-_MODENN__ORDERNN__TNFORDER_-lambda_LAMBDANN_-lambdasrg_LAMBDASRGNN_setNtotmax_NTOMAXVAL_omegaH_OMEGAHVAL_.denshash=_HASHNAME_.v2.0.dat"


def main():
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    df = densdf.pddf
    df = df.reset_index()

    select = (
        (df["N"] == 3)
        & (df["Z"] == 3)
        # & (df["LambdaNN"] == 450)
        & (df["lambdaSRGNN"] == 2.236)
        & (df["kind"] == "one")
    )

    df = df[select]
    ufn = "uniquefilename"
    names = df[ufn].tolist()
    dlNames = []
    # df = df.reset_index()
    # df = df.iloc[0:1]
    # print(df)
    i = 0
    for _, row in df.iterrows():
        kind = row["kind"]
        kindStr = f"-{kind}body"
        dlName = (
            getName(row["Z"], row["N"])
            + kindStr
            + "-"
            + row[ufn]
            + "-hashname="
            + row["hashname"]
            + ".h5"
        )
        dlNames.append(dlName)
        print(f"{i}. ", dlName)
        i += 1
    length = len(dlNames)

    yn = input(f"Download these {length} files?").lower()
    yn = True if yn == "y" else False
    if not yn:
        return
    slash = r"/"

    i = 0
    errorout = []
    for _, row in df.iterrows():
        nucName = "densities-" + getName(row["Z"], row["N"])
        bodyFolder = "2Ndensities" if "twobody" in names else "1Ndensities"
        energy = str(int(np.round(row["omega"])))
        folder = workdir
        folder += slash + nucName
        folder += slash + bodyFolder
        folder += slash + energy + "MeV"
        folder = folder.replace(r"//", r"/")

        if not os.path.isdir(folder):
            print("Directory named\n" + folder + "\nDoes not exist, creating one now")
            os.makedirs(folder)

        fullPath = folder + slash + dlNames[i]
        folder = folder.replace(r"//", r"/")

        try:
            label = {"hashname": row["hashname"], "kind": row["kind"]}
            hashname, _ = densdf.get_file(**label)

            assert os.path.isfile(hashname)
            os.replace(hashname, fullPath)

            print("Downloaded file to:\n" + fullPath)
            print("")
        except HTTPError as err:
            name = dlNames[i]
            hashname = row["hashname"]
            errorout.append([name, err, hashname])

            print(50 * "#")
            print(f"HTTP error was raised for {name} ")
            print(f"Error is {err}")
            print(50 * "#")
            print("")

        i += 1

    if len(errorout) > 0:
        print("Errors raised when attempting to download the following files")
        for i in range(len(errorout)):
            name, err, hashname = errorout[i]
            print("Error created for")
            print(f"name={name}")
            print(f"hash={hashname}")
            print("Gave error", err)
            print("")
            print("")


def getName(Z, N):
    if Z == 3 and N == 3:
        return "6Li"
    if Z == 2 and N == 2:
        return "4He"

    raise ValueError("Impliment me")


def getNuc(name):
    if name == "3H":
        Z = 1
        N = 2
    elif name == "3He":
        Z = 2
        N = 1
    elif name == "4He":
        Z = 2
        N = 2
    elif name == "6Li":
        Z = 3
        N = 3
    else:
        raise ValueError("Wrong name entered")
    return Z, N


def tableCheck(df, hash=None):
    """
    Prints each column name and its corresponding value for the row
    where 'hashname' starts with '405470b'.

    Parameters:
    df (pandas.DataFrame): The DataFrame to process.
    """
    if hash is None:
        hash = "405470b"

    filtered_df = df[df["hashname"].str.startswith(hash, na=False)]

    if filtered_df.empty:
        print(f"No rows found with 'hashname' starting with {hash}")
        return

    # If multiple rows match, we can choose to handle all or just the first one.
    # Here, we'll print the values for each matching row.
    for index, row in filtered_df.iterrows():
        print(f"\nValues for row with index {index}:")
        for col in df.columns:
            print(f"{col}:", f"{row[col]}\n")


nameList = [
    "kind",
    "lambdaSRGNN",
    "omega",
    "theta",
    "LambdaNN",
    "tnforder",
    "OmegaHO",
    "Nmax",
    "MODENN",
    "orderNN",
    "hashname",
]
if __name__ == "__main__":
    main()
