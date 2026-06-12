# Governance

## Maintainer

`DensityScattering` is maintained by Alexander P. Long (corresponding
author), with physics oversight from Harald W. Grießhammer (thesis
supervisor). All decisions about project scope, architecture, and
merges rest with the maintainer.

## Decision model

The project follows a single-maintainer (BDFL) model. This is the
honest description of how decisions are made and reflects the
project's origin as research software supporting a single PhD thesis.

## Scope

`DensityScattering` exists to support calculations associated with
the corresponding PhD thesis and ongoing research at The George
Washington University Institute for Nuclear Studies (GWU INS). The
public release exists to enable:

- reproduction of published results,
- citation by other researchers,
- forking by groups extending the framework to new physics or new
  targets.

The maintainer will not adopt upstream maintenance burden for kernels
or extensions outside the GWU INS scientific stake.

## How changes land

- **Small bug fixes** — typos, Makefile and build-system fixes,
  documentation corrections, small fixes in the Python `tools/`
  helpers — are accepted via pull request. See
  [`CONTRIBUTING.md`](CONTRIBUTING.md) for the workflow.
- **New kernels, new targets, new formalism extensions, or
  substantial performance rewrites** belong in a fork. The maintainer
  is happy to link to scientifically relevant forks from the README.

## Long-term outlook

Maintenance is tied to the maintainer's active research. After
thesis completion (2026), the repository is expected either to
receive only occasional updates or to be archived in a frozen state,
with the Zenodo archive serving as the canonical citation.
Substantial future development is anticipated to occur in forks led
by subsequent students or collaborators in the field.

## Conduct

Contributors are expected to follow ordinary academic norms of
courtesy and scientific honesty. Discussion in GitHub issues and
pull requests should remain respectful and on-topic.
