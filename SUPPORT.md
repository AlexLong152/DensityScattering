# Support

## Maintenance status

`DensityScattering` is actively maintained through completion of the
corresponding PhD thesis at The George Washington University (2026).
After that, the repository is expected to receive only occasional
fixes or to be archived in a frozen state, with the Zenodo archive
serving as the canonical citation.

Substantial new development — new reactions, new targets, new
formalism extensions — is expected to live in **forks**, not upstream.
See [`GOVERNANCE.md`](GOVERNANCE.md) and
[`CONTRIBUTING.md`](CONTRIBUTING.md) for the rationale.

## Reporting bugs and build issues

Please open a GitHub issue including:

- a minimal `input.dat` or run script that reproduces the problem,
- the density file (or a pointer to which file from
  [`nucdens`](https://datapub.fz-juelich.de/anogga/files/densities/)
  was used),
- the full compiler/runtime error message or, for numerical issues,
  the observed and expected output.

Response time is best-effort and may lengthen significantly after
thesis completion.

## Known issues

A list of catalogued imperfections shipped with the current release
is maintained in [`TODO.md`](TODO.md). These items (a Jacobian sign
convention compensated in the parsing layer, duplicated `varsub-*`
modules, unexploited matrix-element symmetries) do not affect the
correctness of the shipped results but are documented openly so a
future maintainer or forker can address them.

## In scope for support

- Build and install failures
- Suspected numerical bugs in the shipped kernels
- Documentation gaps and unclear instructions

## Out of scope

- Porting to platforms or compilers other than those documented in
  `INSTALL.md`
- New-physics or new-target extensions — please fork (see
  [`CONTRIBUTING.md`](CONTRIBUTING.md))
- Debugging user-written kernels in downstream forks

## Citation

See [`CITATION.cff`](CITATION.cff) and the Zenodo archive referenced
in the project metadata. Please cite both the software and the
underlying TDA references listed in the README.

## Scientific collaboration

For collaboration proposals or scientific questions that go beyond
issue-tracker scope, contact the corresponding author. Contact
information appears in `CITATION.cff` and the accompanying JOSS paper.
