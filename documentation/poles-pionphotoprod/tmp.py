def parseSpinString(spinString):
    """
    Parse a string like 'P33nM', 'S11pE', etc. The first 3 chars define the partial-wave:
       e.g.  P(2I)(2J) => 'P33'
    The remainder (e.g. 'nM') defines the sub-channel.

    Return a tuple:
        (plusMinus, ell, I, poleLabel, subChan)
    or None if unphysical or parse fails.

    plusMinus  in {'plus','minus'}
    ell        in {0,1,2,3,...}   # from S,P,D,F,G
    I          = 1/2 or 3/2 typically
    poleLabel  = 'S11','P33', etc. exactly 3 chars
    subChan    = e.g. 'nM' or 'pE'
    """

    if len(spinString) < 4:
        return None

    # 1) extract the partial-wave substring (first 3 chars)
    poleLabel = spinString[:3]  # e.g. 'P33'
    subChan = spinString[3:]  # e.g. 'nM', 'pE', etc.

    # parse the partial-wave
    letter = poleLabel[0].upper()  # e.g. 'P'
    try:
        twoI = int(poleLabel[1])  # e.g. 3 => I=1.5
        twoJ = int(poleLabel[2])  # e.g. 3 => J=1.5
    except ValueError:
        return None

    # Map letter->ell
    L_map = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}
    if letter not in L_map:
        return None
    ell = L_map[letter]

    I = 0.5 * twoI
    J = 0.5 * twoJ

    # plus/minus from J-L
    delta = J - ell
    if abs(delta - 0.5) < 1e-9:
        plusMinus = "plus"
    elif abs(delta + 0.5) < 1e-9:
        plusMinus = "minus"
    else:
        return None  # unphysical

    # We could check if I in {0.5, 1.5} if needed, but maybe we allow any half-integer

    # If we want to skip unphysical, do:
    if not (abs(I - 0.5) < 1e-9 or abs(I - 1.5) < 1e-9):
        # e.g. if I=1.0 => not standard
        return None

    # that's it. We just return
    return (plusMinus, ell, I, poleLabel, subChan)


def buildSpinString(plusMinus, ell, I, subChan):
    """
    The inverse. Given plusMinus, ell, I, and subChan ('nM','pE', etc.),
    build a single string 'P33nM' or 'S11pE' etc.
    Return None if something is unphysical.
    """
    # invert ell-> letter
    L_map_inv = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G"}

    if ell not in L_map_inv:
        return None
    letter = L_map_inv[ell]

    # 2I
    twoI = int(round(2 * I))
    # if we only allow I=0.5 or 1.5
    if twoI not in (1, 3):
        return None

    # plus => J = ell+0.5 => twoJ= 2ell+1
    # minus => J= ell-0.5 => twoJ= 2ell-1
    if plusMinus == "plus":
        twoJ = 2 * ell + 1
    elif plusMinus == "minus":
        twoJ = 2 * ell - 1
    else:
        return None

    if twoJ < 1:
        return None  # e.g. 'minus' for ell=0 => unphysical

    # partial wave label e.g. 'P33'
    pw = f"{letter}{twoI}{twoJ}"

    # subChan appended
    return pw + subChan


# ---------------------
# Example usage / test:

plusmin = {"plus": "+", "minus": "-"}
if __name__ == "__main__":
    samples = ["S11pE", "P33nM", "D15nM", "S31pM", "X99foo", "P11pE", "S12pE"]
    # partWave = ["S11", "S31"]
    # for l, let in enumerate(["P","D","F"]):
    #     for I in [1,3]:
    #         for j in l
    for s in samples:
        out = parseSpinString(s)
        if out is None:
            print(f"{s} -> parseSpinString => None (unphysical or parse fail)")
        else:
            (pm, L, I, poleLbl, subCh) = out
            sign = plusmin[pm]
            tmp = str(L) + sign
            print(f"{s} -> parse ->{subCh}_({tmp})^{I}", end="; ")
            # now rebuild
            rebuild = buildSpinString(pm, L, I, subCh)
            print(f"build => {rebuild}")
