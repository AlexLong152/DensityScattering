import numpy as np
from typing import Callable, Tuple, Dict

"""
|1,1>=|t1=1/2,mt2=1/2>|t2=1/2,mt2=1/2> :index 0
|1,0>=(1/sqrt(2))*(  |t=1/2,mt=1/2>|t2=1/2,mt2=-1/2>  +  |t=1/2,mt=-1/2>|t2=1/2,mt2=1/2>) :index 1
|0,0>=(1/sqrt(2))*(  |t=1/2,mt=1/2>|t2=1/2,mt2=-1/2>  -  |t=1/2,mt=-1/2>|t2=1/2,mt2=1/2>) :index 2
|1,1>=|t1=1/2,mt2=-1/2>|t2=1/2,mt2=-1/2> :index 3
"""


def main():
    extQnum = 1

    # pipiDecoupled(twom1, twom2, twom1p, twom2p, extQnum)

    print("Testing DecoupledFromCoupled vs PionPionBC...")
    print("=" * 60)

    # Test all allowed quantum number combinations
    test_states = [
        (1, 1),  # |1, 1>
        (1, 0),  # |1, 0>
        (1, -1),  # |1, -1>
        (0, 0),  # |0, 0>
    ]

    max_diff = 1e-5
    all_passed = True

    for t, mt in test_states:
        for tp, mtp in test_states:
            # Calculate using coupled basis method
            coupled_result = PionPionBC(t, mt, tp, mtp, extQnum).real

            # Calculate using decoupled basis method
            decoupled_result = DecoupledFromCoupled(t, mt, tp, mtp, extQnum).real

            # Compare results
            diff = abs(coupled_result - decoupled_result)
            max_diff = max(max_diff, diff)

            passed = diff < 1e-10
            all_passed = all_passed and passed
            status = "Pass" if passed else "Fail"
            if not passed:
                print(
                    f"{f'{status} <{tp},{mtp}|O|{t},{mt}>: ':<20}"
                    f"Coupled={coupled_result:+.3f} !- "
                    f"{decoupled_result:+.2f} = Coupled "
                )
    print("=" * 60)
    print(f"Maximum difference: {max_diff:.2e}")
    if all_passed:
        print("All tests PASSED!")
    else:
        print("Some tests FAILED!")

    return all_passed


def DecoupledFromCoupled(t, mt, tp, mtp, extQnum):
    """
    This function is for checking PionPionBC is working correctly.

    Given t, mt, tp, mtp, this function calculates the matrix element using
    pipiDecoupled(twom1, twom2, twom1p, twom2p, extQnum)
    which is returned and then can be checked against PionPionBC
    """

    # Clebsch-Gordan coefficients for two spin-1/2 particles
    # mt values are in units of 1 (not 1/2), so we need to convert
    # Individual particle m values: ±1/2 → twom values: ±1
    invSqrt = 1 / np.sqrt(2)

    def coupled_to_uncoupled(T, MT):
        """
        Returns array of (coef0, coef1, coef2, coef3) for |T,MT> state
        Basis: [|++>, |+->, |-+>, |-->] or [(1,1), (1,-1), (-1,1), (-1,-1)]
        """
        if T == 1:
            if MT == 1:
                # |1,1> = |++>
                return np.array([1, 0, 0, 0], dtype=float)
            elif MT == 0:
                # |1,0> = (1/√2)(|+-> + |-+>)
                return invSqrt * np.array([0, 1, 1, 0], dtype=float)
            elif MT == -1:
                # |1,-1> = |-->
                return np.array([0, 0, 0, 1], dtype=float)
        elif T == 0:
            if MT == 0:
                # |0,0> = (1/√2)(|-+> - |+->)
                return invSqrt * np.array([0, -1, 1, 0], dtype=float)
        raise ValueError(f"Invalid quantum numbers: T={T}, MT={MT}")

    # Expand bra and ket states
    bra_expansion = np.conj(coupled_to_uncoupled(tp, mtp))
    ket_expansion = coupled_to_uncoupled(t, mt)

    # Calculate matrix element as sum over all combinations
    # fill this in
    # Map index to quantum numbers (Kronecker basis): 0->(1,1), 1->(1,-1), 2->(-1,1), 3->(-1,-1)
    twom_states = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    i = 0
    out = 0
    for i, (twom1p, twom2p) in enumerate(twom_states):  # bra side
        ci = bra_expansion[i]
        if ci == 0:
            continue
        for j, (twom1, twom2) in enumerate(twom_states):  # ket side
            cj = ket_expansion[j]
            if cj == 0:
                continue
            me_unc = pipiDecoupled(twom1, twom2, twom1p, twom2p, extQnum)
            out += ci * cj * me_unc
    return out


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

# matNeg1 = (1 / np.sqrt(2)) * (_TAU[1] + 1j * _TAU[2])
# mat1 = (1 / np.sqrt(2)) * (_TAU[1] + 1j * _TAU[2])
# mat0 = _TAU[3]

physicalOpers = {
    -1: (1 / np.sqrt(2)) * (_TAU[1] - 1j * _TAU[2]),
    1: (1 / np.sqrt(2)) * (_TAU[1] + 1j * _TAU[2]),
    0: _TAU[3],
}
# Unitary map between coupled and uncoupled bases.
# Uncoupled/product basis order: |++>, |+->, |-+>, |-->
# Coupled basis order: (1,+1), (1,0), (1,-1), (0,0)


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
    idx = {
        (1, 1): 0,
        (1, 0): 1,
        (1, -1): 2,
        (0, 0): 3,
    }  # in fortran we have to increase index by one bc indicies start at 1

    return idx[(tp, mtp)], idx[(t, mt)]


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

    assert np.shape(oper1) == (2, 2)
    assert np.shape(oper2) == (2, 2)

    total_cpl = np.kron(oper1, oper2)
    return total_cpl


def combinePhysical(charge):
    """
    Helper function, given integer charge=-1,0,1
    returns τ_1^i τ_2^i as a 4x4 matrix
    """
    return combineOpers(physicalOpers[charge], physicalOpers[charge])


def PionPionBC(t, mt, tp, mtp, extQnum, mapping=default_mapping) -> complex:
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
    # charge = extQnum - 2
    # total_cpl += (-2.0) * combineOpers(physicalOpers[charge], physicalOpers[charge])
    total_cpl += (-2.0) * combineOpers(_TAU[extQnum], _TAU[extQnum])

    irow, icol = mapping(t, mt, tp, mtp)
    return total_cpl[irow, icol]


def pipiDecoupled(twom1, twom2, twom1p, twom2p, extQnum):
    """
    Calculates the matrix elements seperately for a cross check
    isospins are |t1,mt1>, |t2,mt2> but we drop the t.
    t1=t2=1/2 so we don't inlcude that. in this function

    the values of m are +-1, with -1 corresponding to index 0, 1 corresponding to index 1 in a length 2 vector

    multiplying m values by two so everything can be an integer
    """

    def getInd(twom):
        # maps -1 to 0, 1 to 1
        return (twom + 1) // 2

    out = 0

    ind1 = getInd(twom1)
    ind2 = getInd(twom2)
    ind1p = getInd(twom1p)
    ind2p = getInd(twom2p)
    for i in range(1, 4):
        tau1 = _TAU[i]
        tau2 = _TAU[i]
        out += tau1[ind1p, ind1] * tau2[ind2p, ind2]
        # out += tau1[ind1, ind1p] * tau2[ind2, ind2p]
    # charge = extQnum - 2
    # tau1a = physicalOpers[charge]
    # tau2a = physicalOpers[charge]
    tau1a = _TAU[extQnum]
    tau2a = _TAU[extQnum]
    out += -2.0 * tau1a[ind1p, ind1] * tau2a[ind2p, ind2]
    # out += -2.0 * tau1a[ind1, ind1p] * tau2a[ind2, ind2p]
    return out


def piPhotoOper(t, mt, tp, mtp, mapping=default_mapping) -> complex:
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
