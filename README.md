# DensityScattering

This project implements the Transition Density Amplitude (TDA) method for
calculating scattering observables off light nuclei (³H, ³He, ⁴He, ⁶Li).
It comprises the work of Alexander P. Long's PhD thesis at The George Washington University.

The TDA approach factorises the scattering amplitude into an irreducible
few-body kernel (encoding probe-nucleon interactions) and a transition
density amplitude (encapsulating nuclear structure). The densities are
computed once per nucleus and then combined with the appropriate kernel for
any elastic reaction, providing a modular framework within chiral effective
field theory (χEFT).

Processes implemented:
- **Compton scattering**: differential cross sections at various photon energies
- **Threshold neutral pion photoproduction**: S-wave multipole amplitudes E₀₊ and L₀₊ (one-body and two-body)
- **Elastic pion-nucleus scattering**: scattering lengths from SAID partial-wave amplitudes and leading two-body diagrams

## Quickstart

For someone who already has HDF5 + h5z-zfp + MPI installed (see
[INSTALL.md](INSTALL.md) otherwise):

```bash
pip install -r requirements.txt && pip install -e tools/package
make -C common-densities/density-modules && make -C PionPion.onebody
python tools/DownloadDens.py        # edit the `select` filter inside first
# Then point input.dat's `densityFileName` line at the downloaded file:
cd PionPion.onebody && ./run.sh
```

The rest of this README explains each step in detail.

## Dependencies

- Fortran compiler with MPI support (`mpif90`, e.g. via OpenMPI or MPICH)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with parallel and Fortran support
- [h5z-zfp](https://github.com/LLNL/H5Z-ZFP) HDF5 filter plugin (the density files distributed for this project are ZFP-compressed and will not load without it)
- GNU Make
- Python 3: not required to run the Fortran code itself, but required
  for the scripting and bulk processing (density download, run
  automation, and analysis)

See [INSTALL.md](INSTALL.md) for source-build instructions for HDF5
and h5z-zfp.

The kernel Makefiles read `HDF5_DIR` from a single `config.mk` at the
repo root (default: `/usr/local/hdf5/openmpi`). Edit that file once, or
override per-shell with `HDF5_DIR=/path/to/hdf5 make`. To locate your
HDF5 prefix:

```bash
pkg-config --variable=prefix hdf5-openmpi-fortran
```

## Building

Each process folder has its own `Makefile`. The shared density modules must
be compiled first. For the non-varsub, one-body systems:

```bash
cd common-densities/density-modules
make

# Then build whichever process you need, e.g.:
cd ../../PionPion.onebody
make
```

The two-body code (`varsub-*.twobody/`) additionally depends on the shared
two-body modules:

```bash
cd common-densities/varsub-density-modules
make

cd ../../varsub-PionPion.twobody
make
```

## Running a calculation

Once a kernel is built, three pieces are needed to run it: a density
file, the Python helpers in `tools/package/` on the path, and an
`input.dat` describing the kinematics.

### 1. Install the Python dependencies and `tools/package` helpers

Install the third-party Python packages used by the download and
analysis scripts (`numpy`, `pandas`, `matplotlib`, `nucdens`, `h5py`,
`requests`):

```bash
pip install -r requirements.txt
```

The kernel runners and analysis helpers (`CrossSection`, `readDensity`,
`runfolderOnebody`, `runfolderTwobody`, `defs`) live in `tools/package/`
and are imported by name. Install the package (`myphysics-local`) once
in editable mode so edits to the source take effect immediately:

```bash
pip install -e tools/package
```

### 2. Download densities

Density files are not bundled with the repository. Fetch them from the
Jülich density store (this uses the
[`nucdens`](https://pypi.org/project/nucdens/) package installed above):

```bash
python tools/DownloadDens.py
```

The script writes to `~/OneDrive/`; edit the `workdir` and the `select`
filter at the top of the file to choose nucleus, energies, angles, and
chiral order.

> **Note**: the Compton helpers in `tools/Compton/` ultimately call into
> `tools/package/{ComptonCross,CrossSection}.py`, both of which hardcode
> the target mass to `M6Li = 6.0151228874 * 931.49432` MeV (the ⁶Li
> atomic mass). Edit that constant if running Compton scattering on a
> different nucleus.

### 3. Run with `run.sh`

Each kernel folder contains a `run.sh` wrapper. With an `input.dat` in
the same directory:

```bash
cd PionPion.onebody
./run.sh
```

`run.sh` rebuilds the kernel (`make clean && make`) and runs the
executable on `input.dat`, writing all `stdout` to the output file
named in `input.dat`.

## Project Structure

| Directory | Description |
|---|---|
| `common-densities/` | Shared density modules and utilities used by all processes |
| `onebodyvia1Ndensity/` | Generic one-body amplitude code (process-independent) |
| `twobodyvia2Ndensity/` | Generic two-body amplitude code (process-independent), without singularity, not maintained |
| **Kernel Examples** | |
| `PionPhotoProd.onebody/` | One-body pion photoproduction (above threshold), using experimentally determined amplitudes |
| `PionPhotoProdThresh.onebody/` | One-body threshold pion photoproduction, using Feynman diagrams |
| `PionPion.onebody/` | One-body elastic pion scattering (above threshold), using experimentally determined amplitudes |
| `varsub-PionPhotoProdThresh.twobody/` | Two-body threshold pion photoproduction, using Feynman diagrams |
| `varsub-PionPion.twobody/` | Two-body elastic pion scattering, using Feynman diagrams |
| **Notes Etc** | |
| `documentation/` | Derivations and notes on kernels, kinematics, and code structure (see [Documentation](#documentation) below) |
| `tools/` | Python analysis and plotting scripts |
| `tests/` | Verification tests |

## Usage

Each process folder contains an `input.dat` file showing the expected input
format. **Note**: the bundled `input.dat` examples have a `densityFileName`
line pointing at the author's local density paths (e.g.
`/home/alexander/OneDrive/densities-.../...h5`); replace it with the path
to your own download from `tools/DownloadDens.py` before running.

A typical run:

```bash
cd PionPion.onebody
make
./run.onebodyvia1Ndensity input.dat
```

Some `run` scripts have a process name
```bash
cd varsub-PionPhotoProdThresh.twobody
make
run.twobodyvia2Ndensity.PionPhotoProdThresh input.dat
```
These `run` scripts are wrappers which write all output sent to `stdout` to an output file specified in the input file.

The python package in `tools` is useful for running many calculations in parallel. 
The file `runfolder.py` present in many kernel folders can be 
used to automatically run all densities in a given folder with a specified format for easy processing.


Density files (not included in this repository) can be downloaded with `tools/DownloadDens.py`

## Input file format

Each kernel folder contains an `input.dat` consumed by the corresponding
`run.*` executable. The file is plain ASCII; values are read positionally,
one record per line. Trailing text on each line is treated as a comment
and is intended for human readers.

### Common lines (all kernels)

| Line | Fields | Meaning |
|---|---|---|
| 1 | `omegaLow omegaHigh omegaStep` | Beam energy scan in MeV (low, high, step) |
| 2 | `thetaLow thetaHigh thetaStep` | Scattering angle scan in degrees (low, high, step) |
| 3 | `outfile` | Output filename. Placeholders `ORDER`, `XXX`, `YYY` are replaced at runtime by calculation order, energy, and angle |
| 4 | `densityFileName` | Absolute path to the HDF5 density file. The target nucleus (`2H`, `3H`, `3He`, `4He`, `6Li`) is detected from this filename |
| 5 | calcstring 1 | Frame, symmetry, verbosity, and `extQnumlimit` (see below) |
| 6 | calcstring 2 | `Calctype` and (for two-body) `j12max`, `numDiagrams` (see below) |

After the calcstrings the format diverges between the one-body and
two-body kernels.

### One-body kernels (`*.onebody/`)

| Line | Fields | Meaning |
|---|---|---|
| 7 | `Nx` | Number of quadrature points for the Feynman-parameter integral in the single-nucleon amplitude |
| 8 | `COMMENTS:_<tag>` | Free-form descriptor appended to the output filename |
| 9+ | `# ...` | Optional comment lines, ignored by the reader |

Example (from `PionPion.onebody/input.dat`):

```
167.664 167.664 10                           omegaLow, omegaHigh, omegaStep
73.091 73.091 15                             thetaLow, thetaHigh, thetaStep
output.dat
'/path/to/densities/onebody-4He.168MeV-073deg.dens-chiralsmsN4LO+3nfN2LO-lambda500-lambdaSRG0.000setNtotmax00omegaH00.denshash=acebce17ba36d4b668d8d80cc9f68807c36ed20b34c853260972107017f0d082.v2.0.h5'
cm_symmetry_verbos_extQnumlimit=3            frame, symmetry, verbosity of STDOUT
Odelta3                                      Calctype
32                                           number of points Nx in Feynman parameter integration
COMMENTS:_v2.0
# in output & density filenames, code replaces ORDER, XXX and YYY automatically by order, energy and angle.
```

### Two-body kernels (`varsub-*.twobody/`)

| Line | Fields | Meaning |
|---|---|---|
| 7 | `NP12A NP12B` | Radial-grid points in the two bins of the (12)-subsystem momentum integral |
| 8 | `P12A P12B P12C` | Bin boundaries in `fm⁻¹`: `[0,P12A/2,P12A]` (hyperbolic map, `NP12A` pts) and `[P12B,P12C]` (linear map, `NP12B` pts) |
| 9 | `AngularType12 Nordth12 Nthbins12 Nordphi12 Nphibins12` | Angular quadrature in (12) subsystem. `AngularType12=1`: Gauss–Legendre on theta and phi separately. `AngularType12=2`: Lebedev–Laikov combined grid; `Nordth12` then sets `Nanggrid12` and the remaining fields are ignored |
| 10 | `COMMENTS:_<tag>` | Descriptor appended to output filename |

Example (from `varsub-PionPion.twobody/input.dat`):

```
131.857 131.857 15                           omegaLow, omegaHigh, omegaStep
000 000 30                                   thetaLow, thetaHigh, thetaStep
out-test.dat
'/path/to/densities/twobody-3He.132MeV-000deg.dens-chiralsmsN4LO+3nfN2LO-lambda400-lambdaSRG0.000setNtotmax00omegaH00.denshash=d306fae448811a1161f693eea1ff4b5cf7d9f9aa523e81e3fef26d4446fcc761.v2.0.h5'
nosymmetry_verbos_extQnumlimit=3             frame, symmetry, verbosity of STDOUT
Odelta4_j12max=1_numDiagrams=1               Calctype, maximal total ang mom in (12) subsystem, and number of diagrams
2 2                                          NP12A, NP12B
1.0 5.0 15                                   P12A, P12B, P12C
2 2                                          AngularType12,(Nordth12 OR Nanggrid12),Nthbins12,Nordphi12,Nphibins12
COMMENTS:_v2.0
# in density filename, code replaces XXX and YYY automatically by energy and angle.
```

### calcstring 1: frame, symmetry, verbosity, `extQnumlimit`

This line is parsed by substring match, so flags may appear in any order
and the leading label (e.g. `cm_symmetry_verbos_extQnumlimit=3`) is only
a mnemonic. Recognised tokens:

- **Frame**: only the centre-of-mass frame is currently implemented; the
  `cm` token is a label for the reader.
- **Symmetry**: `usesymmetry1` … `usesymmetry6` request a symmetry
  reduction (not yet implemented). Omit, or use `nosymmetry`, to compute
  all amplitudes independently.
- **Verbosity**: `verbose2`, `verbose3`, `verbose4` raise the STDOUT
  verbosity. Default is level 1.
- **`extQnumlimit=N`** (`N=1…16`): number of independent external
  quantum-number combinations to compute. Defaults to `4` (Compton: two
  in × two out circular polarisations). Threshold pion photoproduction
  and pion scattering use `extQnumlimit=3`.

### calcstring 2: `Calctype` and two-body flags

`Calctype` selects the chiral order of the kernel. Synonyms in
parentheses are accepted:

| Token | Meaning |
|---|---|
| `Odelta0` (`OQ2`) | `O(q²)` chiral perturbation theory |
| `Odelta2` (`OQ3`) | Full `O(q³)` chiral perturbation theory |
| `Odelta3` (`Oepsilon3`) | Adds explicit Δ(1232) diagrams |
| `Odelta4` (`OQ4`) | `O(q⁴)`; meson-exchange currents are identical for Δ-less and Δ-ful |
| `VaryAXN` (one-body, experimental) | Varies amplitude `AX` of nucleon `N` (`N=p,n`; `X=1…6`) for sensitivity studies |

The `VaryAXN` tokens are really just left over from the Compton
scattering code, so we don't use them for the processes in this
repository.

For two-body kernels the same line may carry additional underscore
separated flags parsed by substring:

- **`j12max=J`** (`J=0…5`): maximum total angular momentum in the (12)
  subsystem. Default `2` for one-body, `1` for two-body (twobody
  converges below 1% by `j12max=1`).
- **`numDiagrams=N`** (`N=1…7`): number of diagrams to include. Default `1`.

### Filename placeholders

Inside `outfile` and `densityFileName`, the strings `ORDER`, `XXX`, and
`YYY` are replaced at runtime by the calculation order, beam energy, and
scattering angle, allowing a single `input.dat` to drive an entire scan.

## Documentation

The `documentation/` folder collects the theoretical notes, derivations,
and writeups that motivate the kernels implemented in the Fortran
modules. The "Related code" column in each table below points to the
implementing files.

### `PionPhotoProduction/`

Threshold and above-threshold neutral pion photoproduction (γ N → π⁰ N).

| File | Purpose | Related code |
|---|---|---|
| `Kernel/pionphoto-kernel.tex/.pdf` | Derivation of the photoproduction kernel | `PionPhotoProd.onebody/`, `PionPhotoProdThresh.onebody/`, `varsub-PionPhotoProdThresh.twobody/` |
| `Kernel/Diag-eval.nb` | Mathematica evaluation of the Feynman diagrams | same as above |
| `Kernel/Finite-Transform.nb` | Finite-momentum transform used in the kernel | same as above |
| `Kernel/SimpleKinematics.nb` | Sanity check of the kinematics | — |
| `Kernel/figs/{1B,2B}-diag*.pdf` | Diagram figures included by `pionphoto-kernel.tex` | — |
| `Kernel/Wirzba-Hanhart.pdf` | Reference paper used in the derivation (third-party) | — |
| `pionphotoproduction-kinematics.pdf` | Kinematics summary | `common-densities/photonmomenta.f` |
| `translating-momentum-transfers-in-densities.v1.0.nb` | Momentum-transfer conventions across density formats | `common-densities/readinput-densities.f` |
| `notes-on-pionproduction.fromAndreas20230320.pdf` | External notes from Andreas Wirzba (third-party) | — |

### `PionScattering/`

Elastic π N scattering, including the pion-pion kernel and isospin algebra.

| File | Purpose | Related code |
|---|---|---|
| `writeup/pi-pi-kernel.tex/.pdf` | Full writeup: kinematics, cross section, scattering length, one-body and two-body diagram contributions | `PionPion.onebody/`, `varsub-PionPion.twobody/` |
| `writeup/[12][a-i].pdf`, `A29-rule.png` | Diagram figures included by `pi-pi-kernel.tex` | — |
| `isospin.nb` | Isospin algebra used in the πN scattering length | `varsub-PionPion.twobody/varsub-2Bkernel.PionPion.f` |

### `twobody-documentation/`

General two-body integration and units, including the variable substitution that handles the pion-pole singularity.

| File | Purpose | Related code |
|---|---|---|
| `twobody-structure-integrations-and-units.20240331.tex/.pdf` | Conventions for the two-body radial/angular integrations and the variable substitution | `varsub-*.twobody/`, `common-densities/varsub-density-modules/` |
| `figuresfeyn.mp` | METAPOST source for the two-body diagrams (gitignored; regenerate with `mpost`) | — |

### `poles-pionphotoprod/`

Notes and checks on the pole approach to photoproduction multipoles (Rijneveen et al.).

| File | Purpose | Related code |
|---|---|---|
| `InvertT.nb` | Inversion of the multipole expansion | `PionPhotoProd.onebody/poles.f`, `onebodyvia1Ndensity/calculateAs.OQ3poles.f` |
| `LegendrePCheck.nb` | Cross-check of the Legendre polynomial expansion | same as above |
| `papers.txt` | Citation list for the references in this folder | — |
| (Other PDFs and `Transfer_PiPh_Rijneveen.zip`) | Third-party reference materials, kept locally; not redistributed | — |

### Top-level documentation files

| File | Purpose | Related code |
|---|---|---|
| `check-trafo1Namps-to-helicity.nb` | Mathematica verification of the single-nucleon amplitude → helicity-basis transformation | `onebodyvia1Ndensity/trafo1Namps-to-helicity.f` |

### Two-body convolution speedup (November 2026)

The two-body mantle in `varsub-twobodyvia2Ndensity/` was reorganised in
November 2026 because production runs had become too slow to be practical.
A single energy and angle for ³He on the usual grid (16 radial points, 50
Lebedev-Laikov angles, `j12max=1`, `Odelta2`) took about 49 minutes. The
same run now takes about one minute, on one core, with the same result up
to round-off. Nothing was parallelised.

**Why it was slow.** The old loop order was channels on the outside and
integration points on the inside. `main` looped over the incoming (12)
quantum numbers and `p12`, `finalstatesums` over the outgoing ones, `p12p`,
`Mzp` and `Mz`, and `Calculate2BIntegralI2` over `msp`, `ms` and finally
the angles. The kernel sat at the very bottom of all of that, so at every
integration point it was evaluated again for every channel pair, every
`(Mzp,Mz)` and every `(msp,ms)`, about 3700 times per point for
`j12max=1`. Most of those calls had exactly the same arguments. On top of
that, each density interpolation first converted a whole slice of
`rhoDensity` (56×56 numbers for the ³He densities) to double precision,
and a NaN check ran over the full integral array in the innermost loop.

**What changed.**

- The loop order is inverted: integration points outside, channels inside.
- At each integration point the kernel is called once for each distinct
  set of its quantum-number arguments (`ml12, ml12p, t12, mt12, t12p,
  mt12p, l12, s12, l12p, s12p`). For `j12max=1` that is 201 calls. `Mz`,
  `Mzp`, `j12` and `m12` are never passed to the kernel, so sharing one
  call between them is exact and not an approximation.
- The spherical harmonic of `p12p` and the position of `|p12p|` on the
  density grid are worked out once per distinct `ppVecs` the kernel
  returns, not once per channel.
- The bilinear interpolation of the density is split in two. The `p12`
  half is done once per radial point, the `p12p` half once per integration
  point. The arithmetic is the same as before, so the interpolated values
  are bit-identical, and the array copy is gone.
- Channel pairs that have no entry in the density are dropped at setup,
  and kernel blocks that come back exactly zero at a point are skipped.
- The NaN check now runs once per `(p12,p12p)`.

**What this assumes about a kernel.** Only that `Calc2Bspinisospintrans`
gives the same output for the same arguments, i.e. it keeps no memory
between calls. Every argument is treated as something the kernel might
depend on, including `l12, ml12, l12p, ml12p`, which no kernel uses at
the moment. Kernels for other processes or higher orders need no changes
in the mantle. The interface is untouched and no kernel file was edited.
One side effect: debug `write` statements inside a kernel now fire far
less often than they used to.

**Where things live now.**

| File | Before | Now |
|---|---|---|
| `varsub-main.twobodyvia2Ndensity.f` | Looped over `mt12, j12, s12, l12, m12, ip12` and called `twobodyfinalstatesumsvia2Ndensity` inside | One call to `twobodyfinalstatesumsvia2Ndensity` per energy and angle |
| `varsub-finalstatesums.twobodyvia2Ndensity.f` | Sum over outgoing channels, `p12p`, `Mzp`, `Mz` | Builds the tables of (12) states, channel pairs and kernel argument sets once, then does both radial integrations. `Calculate2BIntegralI2` is now an internal subroutine here (after `contains`) and does the angular integrals and spin sums for all channel pairs at fixed `(p12,p12p)` |
| `varsub-calculate2BI2.f` | `Calculate2BIntegralI2`, `interpolate`, `getHarmonicCart` | Helpers only: `densitybracket` (replaces `interpolate`) and `getHarmonicCart` |

The header comments of `varsub-finalstatesums.twobodyvia2Ndensity.f`
explain the bookkeeping (states, pairs, kernel argument sets) in more
detail. The older `twobodyvia2Ndensity/` folder was not touched.

**Timings and checks.** Wall time on one core of an i9-14900KS; roughly
3.5 s of every run is spent reading the density file.

| Run | Before | Now |
|---|---|---|
| `varsub-PionPhotoProdThresh.twobody`, 16 radial × 50 angles, `Odelta2` | 2925 s | 58 s |
| same, `Odelta4` | not measured | 130 s |
| `varsub-PionPhotoProdThresh.twobody`, 8 radial × 14 angles, `Odelta2` | 61 s | 4.9 s |
| same, `Odelta4` | 84 s | 6.3 s |
| `varsub-PionPion.twobody`, 8 radial × 14 angles | 279 s | 16 s |

The new code was compared with the old one on six grids for pion
photoproduction (up to 16 radial points or up to 50 angles, both orders)
and one for pion scattering. Differences are below 2×10⁻¹⁴ of the largest
amplitude. The results are not bit-identical because the order of
summation changed. The full 16 × 50 grid was not re-run with the old
code; the 2925 s is the time recorded in an earlier production output.

### Building the writeups

`.tex` writeups are built with `latexmk`:

```bash
cd documentation/PionScattering/writeup
latexmk -pdf pi-pi-kernel.tex
```

The figure dependencies (`[12][a-i].pdf` in `PionScattering/writeup/` and
`figs/*.pdf` in `PionPhotoProduction/Kernel/`) are tracked alongside the
`.tex` so a fresh clone compiles without extra setup. The compiled
writeup PDFs are also tracked for convenience.

## Citation

If you use this code, please cite:

> A. P. Long, *Scattering Observables from Few-Body Densities and Application
> in Light Nuclei*, PhD thesis, The George Washington University (2026).

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

This project is licensed under the [2-Clause BSD License](LICENSE.md).

This code builds on earlier work; see [ACKNOWLEDGMENTS.md](ACKNOWLEDGMENTS.md).
