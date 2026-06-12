# Contributing

Thank you for your interest in `DensityScattering`. Before opening a
pull request, please read this document together with
[`GOVERNANCE.md`](GOVERNANCE.md) to understand the project's scope
and contribution model.

## Two contribution paths

**1. Small upstream contributions — pull requests welcome.** Bug
fixes, documentation corrections, Makefile and build-system fixes,
small improvements to the Python helpers in `tools/`, and similar
focused changes are welcome via pull request.

**2. New-physics contributions — please fork.** Adding a new
reaction kernel, a new nuclear target, a new formalism extension, or
a substantial performance rewrite should happen in a fork rather
than upstream. The maintainer is glad to link to scientifically
relevant forks from the README. See [`GOVERNANCE.md`](GOVERNANCE.md)
for the rationale.

If you are unsure which path applies to your change, open a GitHub
issue describing it before writing code.

## Submitting a small contribution

1. Open a GitHub issue describing the proposed change, unless it is
   a clear, self-contained fix (typo, obvious build error, etc.).
2. Fork the repository and create a topic branch.
3. Make the change. Keep the diff minimal and scoped to the issue.
4. Verify that affected kernels still build (`make` in the relevant
   directory) and, where applicable, that a representative run
   reproduces the prior numerical output to within the documented
   precision.
5. Open a pull request referencing the issue. Briefly describe the
   change, the motivation, and any caveats.

## Coding conventions

- **Fortran:** match the surrounding style. The codebase uses
  fixed-form Fortran organized into modules; new code should follow
  the conventions of the file it sits in.
- **Python helpers (`tools/`):** follow PEP 8.
- **Commits:** one logical change per commit; informative commit
  messages preferred over `fix stuff`.

## Testing

There is no automated unit-test suite. Regression testing is
performed by reproducing published benchmark results — most notably
the threshold pion photoproduction amplitudes of Lenkewitz and
collaborators for $^{3}$H and $^{3}$He, and the corrected $^{6}$Li
threshold results that fix the $\sqrt{2}$ convention error in the
prior literature. If your change touches kernel arithmetic, please
confirm that a representative benchmark still reproduces.

The `tests/` directory contains a small set of focused checks
(`sigmaYcheck`, `IsospinMap.py`) that exercise specific physics
identities used in the kernels.

## Issues and questions

For build problems, suspected bugs, or scientific questions, see
[`SUPPORT.md`](SUPPORT.md).

## Conduct

Contributors are expected to follow ordinary academic norms of
courtesy and scientific honesty.
