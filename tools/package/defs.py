import numpy as np
from typing import Dict

_Sigma: Dict[int, np.ndarray] = {
    1: np.array([[0, 1], [1, 0]], dtype=complex),  # τ^1 = σx
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),  # τ^2 = σy
    # 2: np.array(
    # [[0, 1j], [-1j, 0]], dtype=complex
    # ),  # τ^2 = σy, but its swapped bc of row/coluumn vs column row in fortran
    3: np.array([[1, 0], [0, -1]], dtype=complex),  # τ^3 = σz
}
_TAU = _Sigma

physicalOpers = {
    -1: (1 / np.sqrt(2)) * (_TAU[1] - 1j * _TAU[2]),
    1: (1 / np.sqrt(2)) * (_TAU[1] + 1j * _TAU[2]),
    0: _TAU[3],
}

pauliVec = np.array([_Sigma[1], _Sigma[2], _Sigma[3]])
spinHalfVec = np.array([_Sigma[1], _Sigma[2], _Sigma[3]]) / 2

Sx = (1 / np.sqrt(2)) * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=float)

Sy = (1 / np.sqrt(2)) * np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]], dtype=complex)

Sz = np.array([[1, 0, 0], [0, 0, 0], [0, -0, -1]], dtype=float)


_Spin1: Dict[int, np.ndarray] = {1: Sx, 2: Sy, 3: Sz}
spin1Vec = np.array([Sx, Sy, Sz])


def getSpinOper(TwoSpin):
    """
    Use `TwoSpin` instead of just spin, so that it can be an integer
    """
    match TwoSpin:
        case 0:
            # Just return array full of "1", since it doesn't make sense
            # to have a spin-0 operator, but we still want to read the outputs
            return np.array([1.0, 1.0, 1.0], dtype=complex)
        case 1:
            return spinHalfVec
        case 2:
            return spin1Vec
        case _:
            raise ValueError(f"TwoSpin={TwoSpin} value not supported")


def commute(A, B):
    """
    Calculates [A,B]=AB-BA
    """
    return np.dot(A, B) - np.dot(B, A)


def testCommutes(spinArr):
    """
    Tests:
    [Sx,Sy]=iSz
    [Sy,Sz]=iSx
    [Sz,Sx]=iSy
    """
    passed = True
    Spinx, Spiny, Spinz = spinArr
    t1 = commute(Spinx, Spiny) - 1j * Spinz
    t2 = commute(Spiny, Spinz) - 1j * Spinx
    t3 = commute(Spinz, Spinx) - 1j * Spiny
    strs = ["t1", "t2", "t3"]
    for i, tArray in enumerate([t1, t2, t3]):
        tmp = abs(tArray)
        # print("tmp=\n", tmp)
        # print("tArray=", tArray)
        if np.max(tmp) > 1e-10:
            strTmp = strs[i]
            print(f"Issue for {strTmp}")
            print("tArray=\n", tArray)
            passed = False

    if passed:
        print("Test Passed\n")
    else:
        print("Test Failed\n")


if __name__ == "__main__":
    print("Testing pauli vec")
    testCommutes(spinHalfVec)
    print(75 * "%")
    print("Testing spin 1 vec")
    testCommutes(spin1Vec)
