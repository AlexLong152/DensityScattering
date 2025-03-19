# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import pandas as pd
from CrossSection import ccForDict
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"


def main():
    twobody_dir = r"/home/alexander/Dropbox/COMPTON-RESULTS-FROM-DENSITIES/results-6Li/chiralsmsN4LO+3nfN2LO-lambda550/twobody/"
    onebody_dir = twobody_dir + r"../onebody/"
    savefolder = twobody_dir + r"../results/"
    tmp = np.array(twobody_dir.split(r"/"))
    title = tmp[-3]

    energy = 60
    lambdaSRG = 1.880
    lambdaCut = 550
    omegaH = 14
    Ntotmax = 14
    thetas = np.array([0, 40, 55, 75, 90, 110, 125, 145, 159, 180])
    markers = ["x", ",", "o"]

    outfile = title + "\n"
    outfile += "      -2,    0    2" + "\n"
    total = []
    Deltas = [-1, 0, 1]
    # Deltas = [[1, -1], [0, 0], [2, -2]]#uncomment this to plot difference instead
    isIterable = {True: r"$\alpha+\beta=", False: r"$\alpha-\beta="}
    for i, delta in enumerate(Deltas):
        ys = []
        for theta in thetas:
            ccVal = ccForDict(
                onebody_dir,
                twobody_dir,
                delta=delta,
                energy=energy,
                angle=theta,
                lambdaCut=lambdaCut,
                lambdaSRG=lambdaSRG,
                Ntotmax=Ntotmax,
                omegaH=omegaH,
            )
            ys.append(ccVal)
        total.append(ys)
        alphaBetaStr = isIterable[not isinstance(Deltas[0], list)]
        if isinstance(delta, int):
            plt.scatter(
                thetas,
                ys,
                label=alphaBetaStr + str(2 * delta) + "$",
                marker=markers[i],
            )
        else:
            plt.scatter(
                thetas,
                ys,
                label=alphaBetaStr + str(np.sum([abs(x) for x in delta])) + "$",
                marker=markers[i],
            )
    plt.legend()
    title = (
        "6Li Compton Scattering " + r"varying " + alphaBetaStr + "$\n" + title + "\n"
    )
    title += f"Ntotmax={Ntotmax}, omegaH={omegaH}, lambdaSRG={lambdaSRG}"
    plt.title(title)
    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$d \sigma/ d \Omega\;\; [\mu \mathrm{b}\;\mathrm{ sr^{-1}} ]$")
    plt.show()

    array_out = array_to_table(
        np.array(total).T,
        [str(theta) for theta in thetas],
        [str(x) for x in Deltas],
    )
    out = title
    out += "\nChange of delta in the columns, theta values by rows\n" + str(array_out)
    print(out)
    fileName = title.replace(" ", "-")
    fileName = fileName.replace(",", "")
    fileName = fileName.replace("\n", "")
    fileName = fileName.replace("$", "")
    fileName = fileName.replace("\\", "")
    fileName = fileName.replace("beta=", "beta")
    fileName = fileName.replace("Scattering-", "")
    print("fileName=", fileName)
    save_string_to_file(savefolder + r"/" + fileName + ".txt", out)


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
    return df


def save_string_to_file(file_path: str, content: str):
    """
    Saves the provided string content to a text file at the specified file path.

    Parameters:
        file_path (str): The path (including filename) where the text file will be saved.
        content (str): The string content to write into the file.
    """
    with open(file_path, "w", encoding="utf-8") as file:
        file.write(content)


if __name__ == "__main__":
    main()
