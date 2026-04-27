"""
binding6Li.py

Checks the calculated binding energies of 6Li (N=3, Z=3) densities
from the nucdens database and compares them to the experimental
binding energy of 31.995 MeV.

Prints |E_calc| - |E_exp| for each density along with key parameters.
"""

import os
import pandas as pd
from nucdens import access

# -- Display settings for pandas so the full table is visible --
pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)
pd.set_option("display.width", None)

# -- Experimental binding energy of 6Li in MeV (negative, as a bound state) --
E_EXP = -31.995

# -- Working directory and database URL (same as DownloadDens.py) --
workdir = os.environ["HOME"] + r"/OneDrive/"
beta_webbase = "https://just-object.fz-juelich.de:9000/jikp03/densitystore-beta/"


def main():
    # Step 1: Load the density database into a pandas DataFrame
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    df = densdf.pddf

    # Step 2: Select only 6Li densities (N=3, Z=3)
    select = (df["N"] == 3) & (df["Z"] == 3)
    df = df[select]

    # Step 3: Drop rows that do not have a calculated binding energy.
    # Entries with E[MeV] == 0.0 or NaN are treated as missing.
    df = df[df["E[MeV]"].notna() & (df["E[MeV]"] != 0.0)]

    # Step 4: Compute E_calc - E_exp for each density
    df = df.copy()
    df["E_calc - E_exp [MeV]"] = df["E[MeV]"] - E_EXP

    # Step 5: Filter to the specific parameter values of interest
    df = df[
        (df["Nmax"].isin([14]))
        & (df["OmegaHO"].isin([16, 18]))
        & (df["lambdaSRGNN"].isin([1.880, 2.236, 3.0]))
        & (df["LambdaNN"].isin([400, 450, 500, 550]))
    ]

    # Step 6: Deduplicate by structural parameters (binding energy is the
    # same for all files sharing these values, regardless of omega/theta/kind).
    group_cols = [
        "Nmax",
        "OmegaHO",
        "lambdaSRGNN",
        "LambdaNN",
        "E[MeV]",
        "E_calc - E_exp [MeV]",
    ]
    result = df.groupby(group_cols, dropna=False).size().reset_index(name="_drop")
    result = result.drop(columns=["_drop"])

    # Step 7: Split into two groups and print lambda=400,550 first, then 450,500
    result = result.sort_values(by=["Nmax", "OmegaHO", "lambdaSRGNN", "LambdaNN"])
    result = result.reset_index(drop=True)

    group1 = result[result["LambdaNN"].isin([400, 550])]
    group2 = result[result["LambdaNN"].isin([450, 500])]

    # Step 8: Print the results
    print(f"Experimental binding energy of 6Li: {E_EXP} MeV\n")
    print("LambdaNN = 400, 550:")
    print(group1.to_string(index=False))
    print()
    print("LambdaNN = 450, 500:")
    print(group2.to_string(index=False))


if __name__ == "__main__":
    main()
