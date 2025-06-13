# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import re
from sys import argv


def getFileData():
    """
    Parameters: None

    Returns
    -----------
    dataDict:
    """
    file = "said-pi.txt"

    # Map letters to indices
    letter_map = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}

    # The main nested dictionary
    dataDict = {}

    # Track current section
    current_letter_index = None
    current_2I = None
    current_2L = None

    with open(file, "r") as f:
        for line in f:
            # Detect headers like "S11", "D13", "F15", "G21"
            header_match = re.search(r"\b([SPDFG])(\d)(\d)\b", line)
            #                        ----^  add ‘P’
            if header_match:
                letter, n1, n2 = header_match.groups()
                current_letter_index = letter_map[letter]
                current_2I = int(n1)
                current_2L = int(n2)
                # Initialize nested dicts
                dataDict.setdefault(current_letter_index, {}).setdefault(
                    current_2I, {}
                ).setdefault(current_2L, {})
                continue  # skips the rest of the for loop

            # Parse data rows
            parts = line.split()
            if current_letter_index is not None and len(parts) >= 4:
                try:
                    wcm = float(parts[0])  # WCM column
                    del_val = float(parts[1]) * np.pi / 180  # Del column
                    sr_val = float(parts[3])  # Sr column
                except ValueError:
                    # print("something went wrong")
                    # print("parts=",parts)
                    continue

                # Assign into the nested structure
                dataDict[current_letter_index][current_2I][current_2L][wcm] = np.array(
                    [del_val, sr_val]
                )
    return dataDict


if __name__ == "__main__":
    dataDict = getFileData()
    if len(argv) == 1:
        print("should be: 1.61, 0")
        print(dataDict[0][1][1][1080.00])
        print("should be .02, 0")
        print(dataDict[3][3][7][1296])

        print("should be .18, 0")
        print(dataDict[2][3][7][1236])
    else:
        assert len(argv) == 5
        arg1 = int(argv[1])
        arg2 = int(argv[2])
        arg3 = int(argv[3])
        arg4 = float(argv[4])
        try:
            v1 = dataDict[arg1][arg2][arg3][arg4][0]
            v2 = dataDict[arg1][arg2][arg3][arg4][1]
            print(
                f"Python result:    del={v1:9.5e} sr={v2:9.5e}",
            )
        except KeyError:
            print("Entry not found\n")
