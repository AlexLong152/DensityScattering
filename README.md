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

The Makefiles assume HDF5 is installed at `/usr/local/hdf5/openmpi`.
Adjust the `HDF5_DIR` variable in each `Makefile` if your installation differs.

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
| `documentation/` | Derivations and notes on kernels, kinematics, and code structure |
| `tools/` | Python analysis and plotting scripts |
| `tests/` | Verification tests |

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

## Citation

If you use this code, please cite:

> A. P. Long, *Scattering Observables from Few-Body Densities and Application
> in Light Nuclei*, PhD thesis, The George Washington University (2026).

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

This project is licensed under the [MIT License](LICENSE.md).
