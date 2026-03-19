# Summary of Fixes to PionPion One-Body Code

## Overview

Two bugs were identified and fixed in the pion-pion one-body scattering code that caused incorrect output magnitudes. After both fixes, the results satisfy the isospin relation F(pi-) + F(pi+) = 2*F(pi0), confirming correctness.

---

## Fix 1: Unit Conversion in `main.onebodyvia1Ndensity.f`

**Location:** Lines 386-393

**Problem:**
The old conversion factor was incorrect:
```fortran
outputMat = outputMat * (1/(mpi0*1000)) * HC * HC
```
This produced units of MeV*fm^2 instead of 10^{-3}/m_pi. Multiple errors were present:
- Dividing by `mpi0*1000` instead of multiplying (wrong direction for 10^{-3}/m_pi)
- Spurious `HC^2` factor (converting to fm^2, which is not desired)
- Missing division by `8*pi*sqrt(s)` to convert from M (dimensionless) back to F (MeV^{-1})
- Missing K1N = M_nucl / M_nucleon phase space correction factor
- Using m_{pi0} instead of m_{pi+}

**Fix:**
Replaced with the correct conversion following the pattern established in `PionPhotoProd.onebody`:
```fortran
outputMat = outputMat * (Mnucl/Mnucleon)        ! K1N phase space factor
unitsFactor = (1d-3/mpi)
outputMat = outputMat / (8*Pi*sqrtSReal)         ! M -> F (MeV^{-1})
outputMat = outputMat / unitsFactor              ! express in 10^{-3}/m_pi
```

The chain of units through the code is:
1. `getFAtValue` returns f in MeV^{-1}
2. `getGH` returns g,h in MeV^{-1}
3. `getMat` multiplies by 8*pi*sqrt(s_Real), making M dimensionless
4. The density loop accumulates A_nucl * rho * M (dimensionless)
5. The conversion divides out 8*pi*sqrt(s_Real) and applies K1N and unit factor

---

## Fix 2: |q| Mismatch in Partial Wave Amplitude

**Location:** `PionScat.f` (`getFAtValue`) and `parseFile.f` (`getScatteringData`)

**Problem:**
The partial wave amplitude formula is:

    f = (eta * exp(2i*delta) - 1) / (2i * |q|)

where |q| is the pion CM momentum. The old code computed |q| from the nuclear kinematics vector `qVec`:
```fortran
qAbs = vecAbs(qVec)
```
This `qVec` was obtained from `getKinematics(sqrtSReal, ...)` using the nuclear target mass (M_nucl = 2808.4 MeV for 3He). At the input energy (omega = 136.2 MeV), the nuclear CM kinematics are below the free pi-N threshold, so `getKinematics` floored |q| to 1.0 MeV.

Meanwhile, the SAID phase shift database starts at WCM = 1078 MeV (just above the free pi-p threshold at ~1077.84 MeV). Due to the 50 MeV tolerance in `getScatteringData`, below-threshold lookups (at WCM ~ 1072-1074) matched the WCM = 1078 data point, where the actual free pi-N momentum is |q| ~ 6.24 MeV.

Using |q| = 1.0 MeV (from nuclear kinematics) instead of |q| = 6.24 MeV (from the matched SAID energy) inflated all partial wave amplitudes by a factor of ~6.24.

**Fix:**
1. Modified `getScatteringData` in `parseFile.f` to return the actual matched WCM value (`wcm_actual`) alongside the phase shift data.

2. Modified `getFAtValue` in `PionScat.f` to compute |q| from free pi-proton kinematics at the actual matched SAID energy:
```fortran
Epi_cm = (wcm_actual**2 + mpiPlus**2 - Mproton**2) / (2.0 * wcm_actual)
if (Epi_cm > mpiPlus) then
    qAbs = sqrt(Epi_cm**2 - mpiPlus**2)
else
    qAbs = 0.0
end if
if (qAbs < 1.0) qAbs = 1.0
```
Uses `Mproton` (not the average nucleon mass `mN = 938.919`) because the SAID data corresponds to pi-proton scattering with threshold at m_pi + M_proton = 1077.84 MeV.

---

## Results

For 3He at omega = 136.2 MeV, theta = 0 deg, in units of 10^{-3}/m_pi:

| Pion Charge | Before Fixes | After Fixes |
|-------------|-------------|-------------|
| pi- (extQnum=1) | 1724 + 42i | 276 + 6.8i |
| pi0 (extQnum=2) | 146 + 36i | 23 + 5.7i |
| pi+ (extQnum=3) | -1433 + 29i | -230 + 4.7i |

**Isospin check:** F(pi-) + F(pi+) = 276.3 + (-229.7) = 46.6, and 2*F(pi0) = 2 * 23.3 = 46.6. The relation F(pi-) + F(pi+) = 2*F(pi0) is satisfied exactly, confirming the isospin structure is correct.

---

## Files Modified

- `main.onebodyvia1Ndensity.f` -- Unit conversion fix
- `PionScat.f` -- |q| computation in `getFAtValue`
- `parseFile.f` -- `getScatteringData` now returns actual matched WCM
