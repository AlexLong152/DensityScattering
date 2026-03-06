# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np


Mp = 938.27231
Mn = 939.56563
mpiCharge = 139.57
mpi0 = 134.97
mpi = (mpiCharge + mpi0) / 2

MAve = (Mp + Mn) / 2
M6Li = 6.0151228 * 931.4943
M3He = 2808.4
M3H = 3.01604928 * 931.494
M4He = 3727.4

Ms = np.array([M3H, M3He, M4He, M6Li])
names = ["3H ", "3He", "4He", "6Li"]
Zs = np.array([1, 2, 2, 3])
Ns = np.array([2, 1, 2, 3])

# Central values and uncertainties (all in 10^-3 mpi^-1)
aPlusTilde0 = 1.0
daPlusTilde0 = 1.0
aPlusDelta = -3.35
daPlusDelta = 0.28
aMinusBase = 86.5
daMinusBase = 1.2
aMinusDelta = 1.39
daMinusDelta = 1.33

# Combined central values
aPlus = aPlusTilde0 + aPlusDelta  # -2.35
aMinus = aMinusBase + aMinusDelta  # 87.89

# Combined uncertainties (added in quadrature)
daPlus = np.sqrt(daPlusTilde0**2 + daPlusDelta**2)  # ~1.04
daMinus = np.sqrt(daMinusBase**2 + daMinusDelta**2)  # ~1.79

# Grid of (aPlus, aMinus) at corners of the uncertainty box
aPlusVals = [aPlus - daPlus, aPlus, aPlus + daPlus]
aMinusVals = [aMinus - daMinus, aMinus, aMinus + daMinus]

charge = np.array([-1, 0, 1])
chargeLabels = ["pi-", "pi0", "pi+"]
mpiArr = np.array([mpiCharge, mpi0, mpiCharge])  # pi-, pi0, pi+

if __name__ == "__main__":
    print(f"a+(eff) = {aPlus:.2f} +/- {daPlus:.2f}  [10^-3 mpi^-1]")
    print(f"a-(eff) = {aMinus:.2f} +/- {daMinus:.2f}  [10^-3 mpi^-1]")
    print()

    for i in range(len(names)):
        name = names[i]
        Z = Zs[i]
        N = Ns[i]
        A = Z + N
        print(f"--- {name} (Z={Z}, N={N}, A={A}) ---")
        for ic, Qpi in enumerate(charge):
            label = chargeLabels[ic]
            mpi_ic = mpiArr[ic]
            numer = 1 + mpi_ic / MAve
            denom = 1 + mpi_ic / (A * MAve)
            kinFactor = numer / denom
            vals = []
            for ap in aPlusVals:
                for am in aMinusVals:
                    val = kinFactor * (A * ap - Qpi * (Z - N) * am)
                    vals.append(val)

            central = kinFactor * (A * aPlus - Qpi * (Z - N) * aMinus)
            lo = min(vals)
            hi = max(vals)
            unc = (hi - lo) / 2.0
            Two_T3 = Z - N
            print(
                f"  {label}:  a(1N) = {central:8.2f}   (+/- {unc:.2f}),  T_3={Two_T3}/2"
            )
        print()
