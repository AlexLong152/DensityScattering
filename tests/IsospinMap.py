import numpy as np
from typing import Tuple, Dict
from defs import _TAU, physicalOpers

"""
|1,1>=|t1=1/2,mt2=1/2>|t2=1/2,mt2=1/2> :index 0
|1,0>=(1/sqrt(2))*(  |t=1/2,mt=1/2>|t2=1/2,mt2=-1/2>  +  |t=1/2,mt=-1/2>|t2=1/2,mt2=1/2>) :index 1
|1,-1>=|t1=1/2,mt2=-1/2>|t2=1/2,mt2=-1/2> :index 3
|0,0>=(1/sqrt(2))*(  |t=1/2,mt=1/2>|t2=1/2,mt2=-1/2>  -  |t=1/2,mt=-1/2>|t2=1/2,mt2=1/2>) :index 2
"""


def main():
    test1()
    test2()
    testGetInd()
    testPiPhoto()
    testPiPiA()


def testPiPiA():
    test_states = [
        (1, 1),  # |1, 1>
        (1, 0),  # |1, 0>
        (1, -1),  # |1, -1>
        (0, 0),  # |0, 0>
    ]

    all_passed = True
    tolerance = 1e-5
    # for extQnum in (1, 2, 3):
    extQnum = 3
    for t, mt in test_states:
        for tp, mtp in test_states:
            # Calculate using coupled basis method
            coupled_result = PionPionA(t, mt, tp, mtp, extQnum)

            # Calculate using decoupled basis method
            closed_result = 2 * ClosedPiPhoto(t, mt, tp, mtp, extQnum)

            # Compare results
            diff = abs(coupled_result - closed_result)
            passed = diff < tolerance
            status = "Pass" if passed else "Fail"
            if not passed:
                if all_passed:
                    print("Test 5 Failed: closed form PiPi vs coupled ")
                print(
                    f"{f'{status} <{tp},{mtp}|O|{t},{mt}>: ':<20}"
                    f"Coupled={coupled_result:+.2f} != "
                    f"{closed_result:+.2f} = closed form result "
                )

            all_passed = all_passed and passed
    if all_passed:
        print("Test 5 passed: closed form PiPi vs coupled ")


def testPiPhoto():
    test_states = [
        (1, 1),  # |1, 1>
        (1, 0),  # |1, 0>
        (1, -1),  # |1, -1>
        (0, 0),  # |0, 0>
    ]

    all_passed = True
    tolerance = 1e-5
    # for extQnum in (1, 2, 3):
    extQnum = 3
    for t, mt in test_states:
        for tp, mtp in test_states:
            # Calculate using coupled basis method
            coupled_result = piPhotoOper(t, mt, tp, mtp, extQnum)

            # Calculate using decoupled basis method
            closed_result = ClosedPiPhoto(t, mt, tp, mtp, extQnum)

            # Compare results
            diff = abs(coupled_result - closed_result)
            passed = diff < tolerance
            status = "Pass" if passed else "Fail"
            if not passed:
                if all_passed:
                    print("Test 4 Failed: closed form vs coupled pion photo")
                print(
                    f"{f'{status} <{tp},{mtp}|O|{t},{mt}>: ':<20}"
                    f"Coupled={coupled_result:+.2f} != "
                    f"{closed_result:+.2f} = closed form result "
                )

            all_passed = all_passed and passed
    if all_passed:
        print("Test 4 Passed: closed form vs coupled pion photo")


def test1():
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

    all_passed = True
    tolerance = 1e-5
    for extQnum in (1, 2, 3):
        for t, mt in test_states:
            for tp, mtp in test_states:
                # Calculate using coupled basis method
                coupled_result = PionPionBC(t, mt, tp, mtp, extQnum)

                # Calculate using decoupled basis method
                decoupled_result = DecoupledFromCoupled(t, mt, tp, mtp, extQnum)

                # Compare results
                diff = abs(coupled_result - decoupled_result)
                passed = diff < tolerance
                all_passed = all_passed and passed
                status = "Pass" if passed else "Fail"
                if not passed:
                    print(
                        f"{f'{status} <{tp},{mtp}|O|{t},{mt}>: ':<20}"
                        f"Coupled={coupled_result:+.2f} != "
                        f"{decoupled_result:+.2f} = Coupled "
                    )
    if all_passed:
        print("Test 1 passed: check explicit decoupled vs coupled result")
    else:
        print("In Test 1 Some tests failed")


def test2():
    """
    Tests explicit index, and sanity check on kronecker product
    """

    def op_matrix(extQnum):
        Op = np.zeros((4, 4), dtype=complex)
        for i in (1, 2, 3):
            Op += np.kron(_TAU[i], _TAU[i])
        Op += -2 * np.kron(_TAU[extQnum], _TAU[extQnum])
        return Op

    for extQnum in range(1, 4):
        twom_states = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
        Op = op_matrix(extQnum)

        for bra, (m1p, m2p) in enumerate(twom_states):
            for ket, (m1, m2) in enumerate(twom_states):
                val = pipiDecoupled(m1, m2, m1p, m2p, extQnum)

                assert np.allclose(val, Op[bra, ket])

    def tau2_me(mprime, m):
        return _TAU[2][getInd(mprime), getInd(m)]

    indAnswers = [[(1, -1), -1j], [(-1, 1), 1j]]
    passes = []
    for indAnswer in indAnswers:
        inds, ans = indAnswer
        a, b = inds
        diff = tau2_me(a, b) - ans
        test = abs(diff) < 1e-5
        passes.append(test)
        if not test:
            if all(passes[:-1]):  # first time failing
                print("Test 2 Failed: explicit index checking and kronocker product")
            print(f"tau2_me({a},{b})=", tau2_me(a, b), "!=", ans)

    if all(passes):
        print("Test 2 Passed: explicit index checking and kronocker product")
    else:
        print("\n")


def testGetInd():
    """
    Tests for index consistency
    """
    twom_states = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    all_pass = True

    for extQnum in (1, 2, 3):
        Op = np.zeros((4, 4), dtype=complex)
        for i in (1, 2, 3):
            Op += np.kron(_TAU[i], _TAU[i])
        Op += -2 * np.kron(_TAU[extQnum], _TAU[extQnum])

        for bra, (m1p, m2p) in enumerate(twom_states):
            for ket, (m1, m2) in enumerate(twom_states):
                val = pipiDecoupled(m1, m2, m1p, m2p, extQnum)
                if not np.allclose(val, Op[bra, ket]):
                    if all_pass:
                        print("In test_getInd")
                    print(f"getInd broken at <{m1p},{m2p}|O|{m1},{m2}>")
                    all_pass = False
    if all_pass:
        print("Test 3 Passed: getInd, for index check")
    else:
        print("Test Failed: getInd -- for index check")


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


# Unitary map between coupled and uncoupled bases.
# Uncoupled/product basis order: |++>, |+->, |-+>, |-->
# Coupled basis order: (1,+1), (1,0), (1,-1), (0,0)


def mapping(t: int, mt: int, tp: int, mtp: int) -> Tuple[int, int]:
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


def combineTaus(i):
    return combineOpers(_TAU[i], _TAU[i])


def combineOpers(oper1: np.ndarray, oper2: np.ndarray) -> np.ndarray:
    """
    Combine two spin 1/2 operators, and gives the result in the spin 1 + 0 basis

    |1,1>=|t1=1/2,mt2=1/2>|t2=1/2,mt2=1/2> :index 0
    |1,0>=(1/sqrt(2))*(  |t=1/2,mt=1/2>|t2=1/2,mt2=-1/2>  +  |t=1/2,mt=-1/2>|t2=1/2,mt2=1/2>) :index 1
    |1,-1>=|t1=1/2,mt2=-1/2>|t2=1/2,mt2=-1/2> :index 3
    |0,0>=(1/sqrt(2))*(  |t=1/2,mt=1/2>|t2=1/2,mt2=-1/2>  -  |t=1/2,mt=-1/2>|t2=1/2,mt2=1/2>) :index 2
    """
    invSqrt = 1 / np.sqrt(2)
    U = np.column_stack(
        [
            np.array([1, 0, 0, 0], dtype=complex),  # |1, 1>
            invSqrt * np.array([0, 1, 1, 0], dtype=complex),  # |1, 0>
            np.array([0, 0, 0, 1], dtype=complex),  # |1,-1>
            invSqrt * np.array([0, -1, 1, 0], dtype=complex),  # |0, 0>
        ]
    )
    oper1 = np.asarray(oper1, dtype=complex)
    oper2 = np.asarray(oper2, dtype=complex)

    O_prod = np.kron(oper1, oper2)  # product basis
    O_cpl = U.conj().T @ O_prod @ U  # coupled basis
    return O_cpl


def combinePhysical(extQnum):
    """
    Helper function, given integer charge=-1,0,1
    returns τ_1^i τ_2^i as a 4x4 matrix
    """
    charge = extQnum - 2
    return combineOpers(physicalOpers[charge], physicalOpers[charge])


def PionPionA(t, mt, tp, mtp, extQnum, mapping=mapping) -> complex:
    r"""
    Matrix element in coupled basis of the isotensor-like structure (a=b=extQnum):
        Σ_{i=1}^3 2(τ_1^i τ_2^i)  - 2 (τ_1^a τ_2^a)
    with a = extQnum.

    Factor of 2 comes from a=b

    Implemented term-by-term:
      add  2τ^1⊗τ^1
      add  2τ^2⊗τ^2
      add  2τ^3⊗τ^3
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
    # ExtCharge = extQnum - 2
    for i in (1, 2, 3):
        total_cpl += 2 * combineOpers(_TAU[i], _TAU[i])
        # total_cpl += combinePhysical(i-2)

    # -2 τ_1^a ⊗ τ_2^a
    # total_cpl += (-2.0) * combineOpers(
    #     physicalOpers[ExtCharge], physicalOpers[ExtCharge]
    # )
    total_cpl += (-2.0) * combineOpers(_TAU[extQnum], _TAU[extQnum])

    irow, icol = mapping(t, mt, tp, mtp)
    return total_cpl[irow, icol]


def PionPionBC(t, mt, tp, mtp, extQnum, mapping=mapping) -> complex:
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
    # ExtCharge = extQnum - 2
    for i in (1, 2, 3):
        total_cpl += combineOpers(_TAU[i], _TAU[i])
        # total_cpl += combinePhysical(i-2)

    # -2 τ_1^a ⊗ τ_2^a
    # total_cpl += (-2.0) * combineOpers(
    #     physicalOpers[ExtCharge], physicalOpers[ExtCharge]
    # )
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


def piPhotoOper(t, mt, tp, mtp, extQnum) -> complex:
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
    # total_cpl2 = np.zeros((4, 4), dtype=complex)

    for i in (1, 2, 3):  # x,y only
        total_cpl += combineTaus(i)
    # for i in (-1, 0, 1):
    #     total_cpl2 += combinePhysical(i)
    # print(np.allclose(total_cpl, total_cpl2))
    total_cpl += -combineTaus(extQnum)

    irow, icol = mapping(t, mt, tp, mtp)
    return total_cpl[irow, icol]


def ClosedPiPhoto(t, mt, tp, mtp, extQnum):
    assert extQnum == 3
    val = 2 * (-1) ** (t + 1)
    return val * delta(t, tp) * delta(mt, mtp) * delta(mtp, 0)


def delta(a, b):
    return int(a == b)


def getInd(twom):
    # maps 1 to 0, -1 to 1
    return (1 - twom) // 2


if __name__ == "__main__":
    assert _TAU[2][0, 1] == -1j
    assert _TAU[2][1, 0] == +1j
    main()
