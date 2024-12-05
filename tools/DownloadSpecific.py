import os

import pandas as pd
from nucdens import access
import numpy as np
from datetime import datetime as dt

pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)
energyEps = 0.001
workdir = os.environ["HOME"] + r"/OneDrive/"
beta_webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore-beta/"
webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore/"
"""
Everything needs to be defined at the top in terms of the columns in the database
Each entry in the dict corrisponds to a column in the database
All entries need to support being list or constants, besides energy which must
be a constant
A 'None value means not to filter by this attribute

IMPORTANT: if something needs to be multiple values it needs to be in 
a list, not a numpy array or something else
"""


def main():
    kind = "two"
    # srgback = 2.0
    name = "6Li"
    omega = 60
    theta = [30]

    lambdaSRGNN = 1.88
    LambdaNN = None
    tnforder = None
    ntotmax = None  # actually named nmax in the database
    OmegaHO = None
    MODENN = None  # chsms for example
    orderNN = None  # n4lo+ for example

    """
    dont touch anything below here
    """
    Z, N = getNuc(name)

    dens = {}
    dens["kind"] = kind
    # dens["name"] = name # name actually isn't in the database
    # The name is not a column in the database, instead the name is implied by Z and N

    if kind == "two":
        dens["lambdaSRGNN"] = lambdaSRGNN
        dens["omega"] = omega
        dens["theta"] = theta
        dens["LambdaNN"] = LambdaNN
        dens["tnforder"] = tnforder
        dens["OmegaHO"] = OmegaHO
        dens["Nmax"] = ntotmax
        dens["Z"] = Z
        dens["N"] = N
        dens["MODENN"] = MODENN
        dens["orderNN"] = orderNN

    getDensities(dens, name)


def getDensities(dens, nucName):
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    df = densdf.pddf
    tableCheck(df)  #  prints off all columns and values for one row for easier coding
    length = len(densdf.pddf)
    printState(dens, nucName)

    # selection = np.array([True for x in range(length)])  # init selection
    for key, value in dens.items():
        # print("key,value=", key, value)
        if key not in df.columns:
            raise ValueError(f"Key value {key} not in columns")
        if isinstance(value, type(None)):
            pass
        else:
            # if not isinstance(value, list):
            #     val = [value]
            # else:
            #     val = value
            val = value if isinstance(value, list) else [value]
            if key != "omega":
                df = df.query(f"{key} in @val")
            else:
                if len(val) > 1:
                    raise ValueError("can only handle one energy at a time")
                for en in val:
                    df = df[df["omega"].between(en - energyEps, en + energyEps)]

    labels = []
    names = []
    for i in range(len(df)):
        labels.append(densdf.pddf.to_dict("records")[i])

    ignoreList = ["Z", "N"]
    addList = ["hashname"]

    nameList = [k for (k, _) in dens.items() if k not in ignoreList]

    for val in addList:
        nameList.append(val)

    for i in range(len(df)):
        name = ""
        for n in nameList:
            tmp = df.iloc[i]  # Correct way to access a row
            tmp2 = tmp[n]
            # print("tmp2=", tmp2)
            if isinstance(tmp2, float):
                tmp2 = np.round(tmp2, 2)
            name += f"{n}={tmp2}-"  # Use formatted strings for readability

        name = nucName + "-" + name[:-1] + ".h5"
        name = name.replace("kind=two", "twobody")
        name = name.replace("kind=one", "onebody")
        names.append(name)

    for i, n in enumerate(names):
        print(str(i) + ".", n)
    num = len(names)
    print(f"To download these {num} files enter y or to stop enter n")
    yn = input(
        "Or download specific files enter their numbers on the seperated by a space [y/n]: "
    ).lower()
    if yn != "n":
        # newfolder = input(
        #     "Move new densities to a new subfolder with date-time name?[y/n]: "
        # ).lower()

        # newfolder = True if newfolder == "y" else False
        if yn != "y":
            indexVals = [int(num) for num in yn.split()]
        else:
            indexVals = [i for i in range(len(labels))]

        for i, label in enumerate(labels):
            if i in indexVals:
                if "twobody" in names[i]:
                    subfolder = "2Ndensities"
                else:
                    subfolder = "1Ndensities"

                subfolder += r"/" + str(int(np.round(dens["omega"]))) + "MeV"

                path = workdir + "densities-" + nucName + r"/" + subfolder + r"/"
                fullPath = path + name
                if not os.path.exists(path):
                    print(
                        "Directory named\n"
                        + path
                        + "\nDoes not exist, creating one now"
                    )
                    os.makedirs(path)

                hashname, uniquename = densdf.get_file(**label)
                assert os.path.isfile(hashname)
                os.replace(hashname, fullPath)

                print("Downloaded file to: " + fullPath)


def printState(dens, nucName):
    print("")
    print(f"Searching for {nucName} Densities Matching")
    print(45 * "-")
    items = list(dens.items())
    for i, (k, v) in enumerate(items):
        comma = "," if i < len(items) - 1 else ""
        if v is None:
            vRep = "Any"
        else:
            vRep = v if isinstance(v, str) else repr(v)
        print(f"{k} = {vRep}{comma}")
    print("")


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


onebodOut = "onebody-6Li.060MeV-045deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG2.236setNtotmax14omegaH24.VaryA6p-.denshash=XXXXXXXXXXXXXXXXXXXXXXXXXXXX.v2.0.dat"

twobodOut = "twobody-6Li.060MeV-159deg.dens-chiralsmsN4LO+3nfN2LO-lambda550-lambdaSRG1.880setNtotmax14omegaH24.Odelta2-j12max=2-.denshash=XXXXXXXXXXXXXXXXXXXXXXXXXXXX.v2.0.dat"


if __name__ == "__main__":
    main()
