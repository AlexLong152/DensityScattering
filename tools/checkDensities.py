import numpy as np
from nucdens import access
import pandas as pd
import os

pd.set_option("display.max_rows", 1000)
omegas = np.array([60])
angles = np.array([40, 55, 75, 90, 110, 125, 145, 159])
# Ntotmaxs = np.array([6, 8, 10, 12, 14])
Ntotmaxs = np.array([14])
# omegaHs = np.array([10, 12, 14, 16, 18, 20, 22, 24])
omegaHs = np.array([14])
kinds = ["one", "two"]
Z = 3
N = 3
nucDict = {
    "6Li": [3, 3],
    "3He": [2, 1],
    "3H": [1, 2],
    "4He": [2, 2],
}  # in order of z value, n value
# srgVal = 1.88
srgVal = 2.236
lambdaNNs = np.array([400, 450, 550])
epsSRG = 0.001
eps = 0.5

srg = "lambdaSRGNN"
workdir = os.environ["HOME"] + r"/OneDrive/"
# workdir = "."  # or whatever your system needs
beta_webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore-beta/"


def main():
    os.remove(workdir + "densities_table.h5")
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    df = densdf.pddf
    missing = []
    # for nuc in ["3H", "3He", "4He", "6Li"]:
    for nuc in ["6Li"]:
        print(f"Begining {nuc} requirements")
        print(75 * "-")
        for omega in omegas:
            for lambdaNN in lambdaNNs:
                for theta in angles:
                    for Ntotmax in Ntotmaxs:
                        for omegaH in omegaHs:
                            for kind in kinds:
                                select = (
                                    (df["N"] == N)
                                    & (df["Z"] == Z)
                                    & (df["omega"] > omega - eps)
                                    & (df["omega"] < omega + eps)
                                    & (df[srg] > srgVal - epsSRG)
                                    & (df[srg] < srgVal + epsSRG)
                                    & (df["theta"] > theta - eps)
                                    & (df["theta"] < theta + eps)
                                    & (df["LambdaNN"] == lambdaNN)
                                    & (df["Nmax"] == Ntotmax)
                                    & (df["OmegaHO"] == omegaH)
                                    & (df["kind"] == kind)
                                )
                                dfTmp = df[select]
                                if dfTmp.empty:
                                    missing.append(
                                        [
                                            omega,
                                            theta,
                                            srgVal,
                                            lambdaNN,
                                            Ntotmax,
                                            omegaH,
                                            kind,
                                        ]
                                    )
    outputlabels = ["omega", "theta", "srgVal", "LambdaNN", "Ntotmax", "omegaH", "kind"]
    if len(missing) > 0:
        missing = np.array(missing)
        # print("missing=", missing)
        # print("np.shape(missing)=", np.shape(missing))
        out = array_to_table(missing, None, outputlabels)
        print(out)
        save_string_to_file("missing.txt", out)
    else:
        print("No densities missing")


def array_to_table(array, row_labels, col_labels):
    """
    Converts a MxN numpy array into a table with labeled rows and columns.

    Parameters:
        array (numpy.ndarray): A 2D numpy array of shape (M, N).
        row_labels (list): A list of labels for the rows (length M).
        col_labels (list): A list of labels for the columns (length N).

    Returns:
        pandas.DataFrame: A DataFrame representing the table with the given labels.
    """
    # Check if the array is 2-dimensional.
    if array.ndim != 2:
        raise ValueError("Input array must be 2-dimensional")

    # Check if the number of labels matches the dimensions of the array.
    # if len(row_labels) != array.shape[0]:
    #     raise ValueError(
    #         "The number of row labels must match the number of rows in the array"
    #     )
    if len(col_labels) != array.shape[1]:
        raise ValueError(
            "The number of column labels must match the number of columns in the array"
        )

    # Create a DataFrame with the specified row and column labels.
    df = pd.DataFrame(array, index=row_labels, columns=col_labels)
    return df.to_string(index=False)
    # print(checkIfExists(densDict, 60, 159, 10, 20))


def save_string_to_file(file_path: str, content: str):
    """
    Saves the provided string content to a text file at the specified file path.

    Parameters:
        file_path (str): The path (including filename) where the text file will be saved.
        content (str): The string content to write into the file.
    """
    with open(file_path, "w", encoding="utf-8") as file:
        file.write(content)


# def getComptonKinematics(theta_out, m_i, m_f, M_nucl, E_i_cm):
#     E_probe=m_f*(m_f+2*m_nucl)/(2*m_f


def getAbsk(S, m1, M):
    """
    m1 is probe mass
    M is target (nucleus) mass
    S is mandalstam S
    """
    kSquare = (M**4 + (m1**2 - S) ** 2 - 2 * (M**2) * (m1**2 + S)) / (4 * S)
    return np.sqrt(kSquare)


if __name__ == "__main__":
    main()
