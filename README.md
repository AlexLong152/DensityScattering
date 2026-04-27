# DensityScattering

This project implements the Transition Density Amplitude (TDA) method for
calculating scattering observables off light nuclei (³H, ³He, ⁴He, ⁶Li).
It comprises the work of Alexander P. Long's PhD thesis at The George Washington University.

The TDA approach factorises the scattering amplitude into an irreducible
few-body kernel (encoding probe–nucleon interactions) and a transition
density amplitude (encapsulating nuclear structure). The densities are
computed once per nucleus and then combined with the appropriate kernel for
any elastic reaction, providing a modular framework within chiral effective
field theory (χEFT).

Processes implemented:
- **Compton scattering** — differential cross sections at various photon energies
- **Threshold neutral pion photoproduction** — S-wave multipole amplitudes E₀₊ and L₀₊ (one-body and two-body)
- **Elastic pion–nucleus scattering** — scattering lengths from SAID partial-wave amplitudes and leading two-body diagrams

## Prerequisites

- Fortran compiler with MPI support (`mpif90`, e.g. via OpenMPI or MPICH)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with Fortran bindings
- Python 3 (for analysis scripts and run automation)
- GNU Make

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

Elastic π N scattering, including the pion–pion kernel and isospin algebra.

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

## Usage

Each process folder contains an `input.dat` file showing the expected input
format. A typical run:

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

### Two-body kernels (`varsub-*.twobody/`)

| Line | Fields | Meaning |
|---|---|---|
| 7 | `NP12A NP12B` | Radial-grid points in the two bins of the (12)-subsystem momentum integral |
| 8 | `P12A P12B P12C` | Bin boundaries in `fm⁻¹`: `[0,P12A/2,P12A]` (hyperbolic map, `NP12A` pts) and `[P12B,P12C]` (linear map, `NP12B` pts) |
| 9 | `AngularType12 Nordth12 Nthbins12 Nordphi12 Nphibins12` | Angular quadrature in (12) subsystem. `AngularType12=1`: Gauss–Legendre on theta and phi separately. `AngularType12=2`: Lebedev–Laikov combined grid; `Nordth12` then sets `Nanggrid12` and the remaining fields are ignored |
| 10 | `COMMENTS:_<tag>` | Descriptor appended to output filename |

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

## Citation

If you use this code, please cite:

> A. P. Long, *Scattering Observables from Few-Body Densities and Application
> in Light Nuclei*, PhD thesis, The George Washington University (2026).

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

This project is licensed under the [MIT License](LICENSE.md).
