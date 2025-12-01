import numpy as np

# — replace this with the actual name of your script/module —
from PionScatLib import getKinematics, getF, mpiDict, mNuc

# kinematics setup
sqrtS = 1162  # same √s you’re using
theta = 30 * np.pi / 180  # 30°
x = np.cos(theta)

# π⁻ p channel
piCharge = -1
isospin = 1

# fill in masses exactly as in your main code
m1 = mpiDict[piCharge]
m2 = mNuc[isospin]

# build the CM momentum vector for the initial pion
_, qVec, _ = getKinematics(sqrtS, x, m1, m2, m1, m2)

# partial-wave of interest: ℓ = 1 (P-wave) in the I = 3/2 channel
twoI = 3  # means I = 3/2
ell = 1

# f⁺ is the J = ℓ + 1/2 = 3/2 piece → “P33”
f_P33 = getF(qVec, twoI, ell, +1, sqrtS)

# f⁻ is the J = ℓ – 1/2 = 1/2 piece → “P31”
f_P31 = getF(qVec, twoI, ell, -1, sqrtS)

print("π⁻ p @ 30° (I=3/2):")
print(f"  P33 (f⁺) = {f_P33}")
print(f"  P31 (f⁻) = {f_P31}")
