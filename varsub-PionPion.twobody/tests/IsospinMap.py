import numpy as np
from typing import Callable, Tuple, Dict


def main():
    pass


def allowed_mt(t):
    if t == 0:
        return [0]
    elif t == 1:
        return [-1, 0, 1]
    else:
        raise ValueError("t must be 0 or 1")


# ============================================================
# Two spin/isospin-1/2 particles: coupled <-> uncoupled helpers
# ============================================================

# Pauli matrices (τ^a); use for isospin or spin as needed
_TAU: Dict[int, np.ndarray] = {
    1: np.array([[0, 1], [1, 0]], dtype=complex),  # τ^1 = σx
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),  # τ^2 = σy
    3: np.array([[1, 0], [0, -1]], dtype=complex),  # τ^3 = σz
}

# Unitary map between coupled and uncoupled bases.
# Uncoupled/product basis order: |++>, |+->, |-+>, |-->
# Coupled basis order: (1,+1), (1,0), (1,-1), (0,0)
# Coefficents given by clepsh-gordans, not directly from pion field definitions
_inv_sqrt2 = 1.0 / np.sqrt(2.0)
_U = np.array(
    [
        [1, 0, 0, 0],  # |++>
        [0, _inv_sqrt2, 0, _inv_sqrt2],  # |+->
        [0, _inv_sqrt2, 0, -_inv_sqrt2],  # |-+>
        [0, 0, 1, 0],  # |-->
    ],
    dtype=complex,
)

_U = np.array(
    [
        [1, 0, 0, 0],  # |++>
        [0, _inv_sqrt2, 0, 1j * _inv_sqrt2],  # |+->
        [0, _inv_sqrt2, 0, -1j * _inv_sqrt2],  # |-+>
        [0, 0, 1, 0],  # |-->
    ],
    dtype=complex,
)

_U = np.array(
    [
        [1, 0, 0, 0],  # |++>
        [0, _inv_sqrt2, 0, 1j * _inv_sqrt2],  # |+->
        [0, _inv_sqrt2, 0, -1j * _inv_sqrt2],  # |-+>
        [0, 0, 1, 0],  # |-->
    ],
    dtype=complex,
)


def default_mapping(t: int, mt: int, tp: int, mtp: int) -> Tuple[int, int]:
    """
    Map coupled quantum numbers (t, mt), (t', mt') to matrix indices
    in the coupled-basis ordering:
        (1,+1), (1,0), (1,-1), (0,0)

    mt are in units of 1:
      - if t=1: mt ∈ {+1, 0, -1}
      - if t=0: mt = 0
    Returns (row_index, col_index) = (bra_index, ket_index).
    """
    idx = {(1, 1): 0, (1, 0): 1, (1, -1): 2, (0, 0): 3}
    if (t, mt) not in idx or (tp, mtp) not in idx:
        raise ValueError(
            f"Passed t,mt={t},{mt} and tp, mtp={tp},{mtp} "
            + "Invalid (t,mt). Use t=1 with mt in {+1,0,-1} or t=0 with mt=0."
        )
    return idx[(tp, mtp)], idx[(t, mt)]


def GetOperFromCharge(charge):
    if charge == -1:
        mat = (1 / np.sqrt(2)) * (_TAU[1] + 1j * _TAU[2])
    elif charge == 1:
        mat = (1 / np.sqrt(2)) * (_TAU[1] + 1j * _TAU[2])
    else:
        mat = _TAU[3]
    return mat


def combineOpers(oper1: np.ndarray, oper2: np.ndarray) -> np.ndarray:
    """
    Combine two single-particle 2x2 operators into a two-body operator in the COUPLED basis.

    Steps:
      1) Build uncoupled/product-space operator: oper1 ⊗ oper2  (4x4)
      2) Rotate to coupled basis: U† (oper1 ⊗ oper2) U         (4x4)

    Returns
    -------
    (4,4) complex ndarray
        The combined operator in the coupled basis.
    """
    oper1 = np.asarray(oper1, dtype=complex)
    oper2 = np.asarray(oper2, dtype=complex)
    if oper1.shape != (2, 2) or oper2.shape != (2, 2):
        raise ValueError(
            f"oper1 and oper2 must be 2x2. Got {oper1.shape} and {oper2.shape}."
        )

    total_cpl = np.kron(oper1, oper2)
    # total_cpl = _U.conj().T @ total_cpl @ _U
    return total_cpl


def PionPionBC(
    t: int,
    mt: int,
    tp: int,
    mtp: int,
    extQnum: int,
    mapping: Callable[[int, int, int, int], Tuple[int, int]] = default_mapping,
) -> complex:
    r"""
    Matrix element in coupled basis of the isotensor-like structure (a=b=extQnum):
        Σ_{i=1}^3 (τ_1^i τ_2^i)  - 2 (τ_1^a τ_2^a)
    with a = extQnum.

    Factor of 2 comes from a=b

    Implemented term-by-term:
      add  τ^1⊗τ^1
      add  τ^2⊗τ^2
      add  τ^3⊗τ^3
      add -2 τ^a⊗τ^a

    Returns
    -------
    complex
        <t' mt' | O | t mt>
    """
    if extQnum not in (1, 2, 3):
        raise ValueError("extQnum must be 1, 2, or 3.")

    total_cpl = np.zeros((4, 4), dtype=complex)

    # Σ_{i=1}^3 τ_1^i ⊗ τ_2^i
    for i in (1, 2, 3):
        total_cpl += combineOpers(_TAU[i], _TAU[i])

    # -2 τ_1^a ⊗ τ_2^a
    total_cpl += (-2.0) * combineOpers(_TAU[extQnum], _TAU[extQnum])

    irow, icol = mapping(t, mt, tp, mtp)
    return total_cpl[irow, icol]


def piPhotoOper(
    t: int,
    mt: int,
    tp: int,
    mtp: int,
    mapping: Callable[[int, int, int, int], Tuple[int, int]] = default_mapping,
) -> complex:
    r"""
    Matrix element in coupled basis of:
        (τ_1 · τ_2 - τ_1^z τ_2^z)
      = τ_1^x τ_2^x + τ_1^y τ_2^y
      = τ_1^1 τ_2^1 + τ_1^2 τ_2^2

    Implemented term-by-term:
      add τ^1⊗τ^1
      add τ^2⊗τ^2

    Returns
    -------
    complex
        <t' mt' | O | t mt>
    """
    total_cpl = np.zeros((4, 4), dtype=complex)

    for i in (1, 2):  # x,y only
        total_cpl += combineOpers(_TAU[i], _TAU[i])

    irow, icol = mapping(t, mt, tp, mtp)
    return total_cpl[irow, icol]


if __name__ == "__main__":
    main()
