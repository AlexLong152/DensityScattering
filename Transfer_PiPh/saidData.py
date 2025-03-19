#!/usr/bin/env python3
import re
import numpy as np
from collections import defaultdict


filename = "said-SM22.txt"  # Replace with the actual data filename


def main():
    data_dict = load_pi0p_data(filename)

    # Show some of the structure
    for spin in sorted(data_dict.keys()):
        for targ in sorted(data_dict[spin].keys()):
            print(f"\n=== {spin}, {targ} ===")
            E = data_dict[spin][targ]["energy"]
            cen = data_dict[spin][targ]["center"]
            err = data_dict[spin][targ]["err"]

            print(f"E array shape:   {E.shape}")
            print(f"center shape:    {cen.shape}")
            print(f"err shape:       {err.shape}")

            # Print first few lines as a sample
            for i in range(min(3, len(E))):
                print(
                    f"  E={E[i]:.2f}  center=({cen[i].real:.3f} + {cen[i].imag:.3f}j)"
                    f"  err=({err[i].real:.3f} + {err[i].imag:.3f}j)"
                )


def parse_partial_wave_block(lines):
    """
    Given a list of lines corresponding to one partial-wave block (e.g. from
    after "PI0P S11 pE..." until the next partial-wave header or end),
    parse the multi-energy lines (not in "GWU" block) and the single-energy
    lines (in "GWU Single energy values" block). Return a dict:

      blockDict[energy_float] = {
         "center": complex(...) or None,
         "err":    complex(...) or None
      }

    This function doesn't care about spinString or targetString. It just
    extracts data from lines of either type. Missing lines => no 'center' or 'err'.
    """

    multi_line_regex = re.compile(
        r"^\s*(?P<EG>\S+)\s+(?P<EMr>\S+)\s+(?P<EMi>\S+)\s+"
        r"(?P<BORN>\S+)\s+(?P<Phs>\S+)\s+(?P<Xsect>\S+)"
    )

    single_line_regex = re.compile(
        r"^\s*(?P<EG>\S+)\s+"
        r"(?P<rCenter>[\-\+\d\.]+)\(\s*(?P<rErr>[\-\+\d\.]+)\)\s+"
        r"(?P<iCenter>[\-\+\d\.]+)\(\s*(?P<iErr>[\-\+\d\.]+)\)"
    )

    blockDict = {}
    in_gwu_block = False

    for line in lines:
        line_stripped = line.strip()

        # Check if "GWU Single energy values" has begun
        if "GWU" in line_stripped and "Single energy values" in line_stripped:
            in_gwu_block = True
        else:
            if not in_gwu_block:
                # Attempt multi-energy parse
                m_mult = multi_line_regex.match(line_stripped)
                if m_mult:
                    try:
                        EG_val = float(m_mult.group("EG"))
                        EMr_val = float(m_mult.group("EMr"))
                        EMi_val = float(m_mult.group("EMi"))
                    except ValueError:
                        EG_val = None
                    if EG_val is not None:
                        if EG_val not in blockDict:
                            blockDict[EG_val] = {"center": None, "err": None}
                        blockDict[EG_val]["center"] = complex(EMr_val, EMi_val)

            else:
                # We are in the GWU single-energy block
                m_single = single_line_regex.match(line_stripped)
                if m_single:
                    try:
                        EG_val = float(m_single.group("EG"))
                        rC_val = float(m_single.group("rCenter"))
                        rE_val = float(m_single.group("rErr"))
                        iC_val = float(m_single.group("iCenter"))
                        iE_val = float(m_single.group("iErr"))
                    except ValueError:
                        EG_val = None
                    if EG_val is not None:
                        if EG_val not in blockDict:
                            blockDict[EG_val] = {"center": None, "err": None}
                        # Overwrite/supplement center
                        blockDict[EG_val]["center"] = complex(rC_val, iC_val)
                        blockDict[EG_val]["err"] = complex(rE_val, iE_val)

    return blockDict


def load_pi0p_data(filename):
    """
    Parse the entire pi0p data file into a nested dictionary:

      myDict[spinString][targetString]["energy"] -> 1D array of energies
      myDict[spinString][targetString]["center"] -> 1D complex array
      myDict[spinString][targetString]["err"]    -> 1D complex array

    For each partial-wave header like "PI0P S11 pE ...", we accumulate lines
    until the next partial-wave header or end, pass them to parse_partial_wave_block,
    and store the results in temp[spinString][targetString][energy].
    Finally, unify energies across targets for each spin, inserting 0 if missing.
    """

    # Regex that identifies the partial-wave header line, e.g.
    # "PI0P S11  pE        1/21/25"
    header_regex = re.compile(r"^\s*PI0P\s+([SDPFG]\d{1,2})\s+([pn][EM])")

    # We will read the file lines first
    with open(filename, "r") as f:
        all_lines = f.read().splitlines()

    # We'll store them in a structure temp[spin][targ][energy] = { center, err }
    temp = defaultdict(lambda: defaultdict(dict))

    # We'll collect lines for each partial-wave block
    current_spin = None
    current_targ = None
    current_block = []
    blocks = []

    # We'll iterate over lines, detect new headers, and store blocks
    for line in all_lines:
        line_stripped = line.strip()
        m_head = header_regex.search(line_stripped)
        if m_head:
            # We found a new partial-wave header
            # 1) if we already had a block in progress, parse it
            if current_block and current_spin and current_targ:
                blockDict = parse_partial_wave_block(current_block)
                # Merge blockDict into temp
                for E_val, subdict in blockDict.items():
                    if E_val not in temp[current_spin][current_targ]:
                        temp[current_spin][current_targ][E_val] = {
                            "center": None,
                            "err": None,
                        }
                    # Possibly overwrite center
                    if subdict["center"] is not None:
                        temp[current_spin][current_targ][E_val]["center"] = subdict[
                            "center"
                        ]
                    # Possibly overwrite err
                    if subdict["err"] is not None:
                        temp[current_spin][current_targ][E_val]["err"] = subdict["err"]

            # 2) start a new block
            current_spin = m_head.group(1)
            current_targ = m_head.group(2)
            current_block = []
        else:
            # not a new header line => belongs to the current block
            current_block.append(line)

    # after the loop ends, parse the final block if any
    if current_block and current_spin and current_targ:
        blockDict = parse_partial_wave_block(current_block)
        for E_val, subdict in blockDict.items():
            if E_val not in temp[current_spin][current_targ]:
                temp[current_spin][current_targ][E_val] = {"center": None, "err": None}
            if subdict["center"] is not None:
                temp[current_spin][current_targ][E_val]["center"] = subdict["center"]
            if subdict["err"] is not None:
                temp[current_spin][current_targ][E_val]["err"] = subdict["err"]

    # Now unify energies across each spin => fill zeros for missing data
    final_dict = {}
    for spin, targ_dict in temp.items():
        # gather union of energies for all targets in this spin
        union_energies = set()
        for targ, e_dict in targ_dict.items():
            union_energies.update(e_dict.keys())
        union_energies = sorted(union_energies)

        final_dict[spin] = {}
        for targ, e_dict in targ_dict.items():
            center_list = []
            err_list = []
            for E_val in union_energies:
                if E_val not in e_dict:
                    # missing => set 0
                    center_list.append(0 + 0j)
                    err_list.append(0 + 0j)
                else:
                    cval = e_dict[E_val]["center"]
                    eval_ = e_dict[E_val]["err"]
                    if cval is None:
                        cval = 0 + 0j
                    if eval_ is None:
                        eval_ = 0 + 0j
                    center_list.append(cval)
                    err_list.append(eval_)

            # convert to arrays
            E_arr = np.array(union_energies, dtype=float)
            c_arr = np.array(center_list, dtype=complex)
            e_arr = np.array(err_list, dtype=complex)

            final_dict[spin][targ] = {
                "energy": E_arr,
                "center": c_arr,
                "err": e_arr,
            }

    return final_dict


if __name__ == "__main__":
    main()
