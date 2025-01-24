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
webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore/"

filename = "_KIND_body-_NUC_._ENERGY_MeV-_ANGLE_deg.dens-_MODENN__ORDERNN__TNFORDER_-lambda_LAMBDANN_-lambdaSRG_LAMBDASRGNN_setNtotmax_NTOMAXVAL_omegaH_OMEGAHVAL_.denshash=_HASHNAME_.v2.0.h5"


def main():
    # User interface: the `select` variable filters the DataFrame
    # Adjust this selection as needed
    # os.system(f"rm {workdir}//densities_table.h5")
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    # densdf = access.database(workdir=workdir, webbase=webbase)
    df = densdf.pddf
    df = df.reset_index()
    tableCheck(df)
    eps = 0.0001
    omega = 60
    srg = "lambdaSRGNN"
    srgVal = 1.880
    # srgVal = 2.236

    LambdaNNs = np.array([400, 550])
    thetas = np.array([40, 55, 75, 90, 110, 125, 145, 159])
    Nmaxs = np.array([6, 8, 10, 12, 14])
    omegaHOs = np.array([10, 12, 14, 16, 18, 20, 22, 24])
    kinds = np.array(["one", "two"])
    j = 0
    i = 0
    for kind in kinds:
        for theta in thetas:
            for LambdaNN in LambdaNNs:
                for Nmax in Nmaxs:
                    for omegaHO in omegaHOs:
                        select = (
                            (df["N"] == 3)
                            & (df["Z"] == 3)
                            & (df["omega"] > omega - eps)
                            & (df["omega"] < omega + eps)
                            & (df[srg] > srgVal - eps)
                            & (df[srg] < srgVal + eps)
                            & (df["LambdaNN"] == LambdaNN)
                            & (df["kind"] == kind)
                            & (df["OmegaHO"] == omegaHO)
                            & (df["Nmax"] == Nmax)
                            & (df["theta"] == theta)
                        )

                        dfTmp = df[select]
                        if len(dfTmp) == 0:
                            print(
                                "kind='",
                                kind,
                                "', lambdaSRGNN=",
                                srgVal,
                                ", theta=",
                                theta,
                                ", LambdaNN=",
                                LambdaNN,
                                ", Nmax=",
                                Nmax,
                                ", OmegaHO=",
                                omegaHO,
                                sep="",
                            )
                            j += 1
                        else:
                            i += 1

    print(j, "missing files")
    print(i, "found files")


def getName(Z, N):
    if Z == 3 and N == 3:
        return "6Li"
    if Z == 2 and N == 2:
        return "4He"
    if Z == 2 and N == 1:
        return "3He"

    raise ValueError("Implement me")


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
