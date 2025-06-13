import os
import pandas as pd
from nucdens import access
import numpy as np
from requests.exceptions import HTTPError

pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)

workdir = os.environ["HOME"] + r"/OneDrive/"
beta_webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore-beta/"
beta_webbase = "https://just-object.fz-juelich.de:9000/jikp03/densitystore-beta/"

# webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore/"

filename = "_KIND_body-_NUC_._ENERGY_MeV-_ANGLE_deg.dens-_MODENN__ORDERNN__TNFORDER_-lambda_LAMBDANN_-lambdaSRG_LAMBDASRGNN_setNtotmax_NTOMAXVAL_omegaH_OMEGAHVAL_.denshash=_HASHNAME_.v2.0.h5"


def main():
    # User interface: the `select` variable filters the DataFrame
    # Adjust this selection as needed
    # os.system(f"rm {workdir}//densities_table.h5")
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    df = densdf.pddf
    df = df.reset_index()
    # tableCheck(df)
    ts1 = pd.Timestamp("2025-04-26 1:00:00")

    select = (
      (df["N"] == 3)
    & (df["LambdaNN"]>10)
    & (df["Z"] == 3)
    & (df["addtime"] < ts1)
    & (df["omega"]<60.1)
    & (df["omega"]>59.9)
    & (df["theta"]==40)
    & (df["LambdaNN"]==400)
    & (df["Nmax"]==10.0)
    & (df["OmegaHO"]==16.0)
    #& (df["lambdaSRGNN"]==1.880)
)
    df = df[select]

    dlNames = []
    i = 0
    for _, row in df.iterrows():
        dlName = buildFilename(row)
        dlNames.append(dlName)
        print(f"{i}. {dlName}")
        i += 1
    length = len(dlNames)

    yn = input(f"Download these {length} files? (y/n): ").lower()
    yn = True if yn == "y" else False
    if not yn:
        return

    slash = r"/"
    i = 0
    errorout = []

    nuc = getName(row["Z"], row["N"])
    if nuc == "4He" or nuc == "6Li":
        yn = input(
            "Place these files into seperate folders based on lambda/lambdaSRG? (y/n): "
        ).lower()
        splitFolder = True if yn == "y" else False
    else:
        splitFolder = False

    for _, row in df.iterrows():
        nucName = "densities-" + getName(row["Z"], row["N"])
        bodyFolder = "2Ndensities" if row["kind"] == "two" else "1Ndensities"
        energy = str(int(np.round(row["omega"])))
        srg = str(row["lambdaSRGNN"])

        folder = workdir
        folder += slash + nucName
        folder += slash + bodyFolder
        folder += slash + energy + "MeV"
        if splitFolder:
            folder += slash + "lambda" + str(int(row["LambdaNN"]))
            folder += slash + "lambdaSRG" + srg

        folder = folder.replace(r"//", r"/")

        if not os.path.isdir(folder):
            print("Directory named\n" + folder + "\nDoes not exist, creating one now")
            os.makedirs(folder)

        fullPath = folder + slash + dlNames[i]
        folder = folder.replace(r"//", r"/")
        if not os.path.exists(fullPath):
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
        else:
            print(f"File already exists with name:\n{fullPath}\n")

        i += 1

    if len(errorout) > 0:
        print("Errors raised when attempting to download the following files")
        for name, err, hashname in errorout:
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
    if Z == 2 and N == 1:
        return "3He"
    if Z == 1 and N == 2:
        return "3H"

    raise ValueError(f"Implement me for Z={Z}, N={N}")


def buildFilename(row):
    # Extract and prepare values
    kind = row["kind"]  # 'one' or 'two'
    # Map _KIND_ to create 'onebody' or 'twobody'
    # The filename template uses "_KIND_body" so replacing _KIND_ with `kind` gives "onebody" or "twobody"
    # For example: if kind='one', _KIND_body -> 'onebody'

    Nuc = getName(row["Z"], row["N"])

    # Energy and angle
    # Rounding to int as in original code
    Energy = str(int(np.round(row["omega"])))
    Energy = Energy.zfill(3)
    Angle = str(int(np.round(row["theta"])))
    Angle = Angle.zfill(3)

    # MODENN mapping
    # If MODENN is 'chsms', map to 'chiralSMS', else use as is
    MODENN = row["MODENN"] if "MODENN" in row and not pd.isna(row["MODENN"]) else ""
    if MODENN == "chsms":
        MODENN = "chiralsms"

    # orderNN mapping
    orderNN = row["orderNN"] if "orderNN" in row and not pd.isna(row["orderNN"]) else ""
    if orderNN == "n4lo+":
        orderNN = "N4LO+"

    # tnforder mapping
    tnforder = (
        row["tnforder"] if "tnforder" in row and not pd.isna(row["tnforder"]) else ""
    )
    tnforder = tnforder.replace("-combine", "")
    if tnforder == "n2lo":
        tnforder = "3nfN2LO"

    # LambdaNN and lambdaSRGNN
    LambdaNN = ""
    if "LambdaNN" in row and not pd.isna(row["LambdaNN"]):
        # Convert to int if needed
        LambdaNN = str(int(float(row["LambdaNN"])))

    lambdaSRGNN = ""
    if "lambdaSRGNN" in row and not pd.isna(row["lambdaSRGNN"]):
        lambdaSRGNN = row["lambdaSRGNN"]
        lambdaSRGNN = f"{lambdaSRGNN:.3f}"

    if "Nmax" in row and not pd.isna(row["Nmax"]):
        NTOMAXVAL = int(row["Nmax"])
        NTOMAXVAL = f"{NTOMAXVAL:02d}"
    else:
        NTOMAXVAL = ""

    if "OmegaHO" in row and not pd.isna(row["OmegaHO"]):
        OMEGAHVAL = int(row["OmegaHO"])
        OMEGAHVAL = f"{OMEGAHVAL:02d}"
    else:
        OMEGAHVAL = ""

    # Hashname
    full_hash = row["hashname"]
    # filename template adds ".v2.0.dat" at the end anyway, so we don't need .h5

    # Perform the replacements in the filename template
    outString = filename
    outString = outString.replace("_KIND_", kind)
    outString = outString.replace("_NUC_", Nuc)
    outString = outString.replace("_ENERGY_", Energy)
    outString = outString.replace("_ANGLE_", Angle)
    outString = outString.replace("_MODENN_", MODENN)
    outString = outString.replace("_ORDERNN_", orderNN)
    outString = outString.replace("_TNFORDER_", tnforder)
    outString = outString.replace("_LAMBDANN_", LambdaNN)
    outString = outString.replace("_LAMBDASRGNN_", lambdaSRGNN)
    outString = outString.replace("_NTOMAXVAL_", NTOMAXVAL)
    outString = outString.replace("_OMEGAHVAL_", OMEGAHVAL)
    outString = outString.replace("_HASHNAME_", full_hash)

    return outString


def tableCheck(df, hash=None):
    """
    Prints each column name and its corresponding value for the row
    where 'hashname' starts with '405470b'.

    Parameters:
    df (pandas.DataFrame): The DataFrame to process.
    """
    if hash is None:
        hash = "a35b77934a5f5e5625e12c610795f23deab89c9e6572c435239cb4de1f71dcbf"
    """
    Prints each column name and its corresponding value for the row
    where 'hashname' starts with '405470b'.

    Parameters:
    df (pandas.DataFrame): The DataFrame to process.
    """

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
    val = filtered_df["lambdaSRGNN"]
    # print("val=", val)


if __name__ == "__main__":
    main()
