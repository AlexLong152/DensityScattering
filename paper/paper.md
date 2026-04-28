---
title: 'DensityScattering: A modular Fortran framework for elastic scattering off light nuclei via transition density amplitudes'
tags:
  - Fortran
  - nuclear physics
  - chiral effective field theory
  - few-body systems
  - light nuclei
  - Compton scattering
  - pion photoproduction
  - pion-nucleus scattering
authors:
  - name: Alexander P. Long
    orcid: 0009-0009-5890-2713
    corresponding: true
    affiliation: 1
affiliations:
  - name: Institute for Nuclear Studies, Department of Physics, The George Washington University, Washington DC, USA
    index: 1
date: 28 April 2026
bibliography: paper.bib
---

# Summary

`DensityScattering` is a Fortran code suite that computes elastic
scattering observables off light nuclei within the Transition Density
Amplitude (TDA) formalism of [@Hammer3He]. The TDA approach
factorizes the nuclear scattering amplitude into a process-specific
irreducible few-body kernel, encoding the interaction of an external
probe with the active nucleons, and a target-specific transition
density amplitude that encapsulates the nuclear structure of the
spectator nucleons. Densities are computed once per nucleus and then
convolved with any kernel of interest, so that a researcher
implementing a new reaction supplies only the process-specific kernel.
The package currently provides kernel implementations for Compton
scattering, neutral pion photoproduction (at and above threshold),
and elastic pion-nucleus scattering, with applications to ${}^3$H,
${}^3$He, ${}^4$He, and ${}^6$Li.

# Statement of need

Few-body calculations of scattering observables in chiral effective
field theory ($\chi$EFT) [@chiEFT] are computationally expensive and
are typically tied to a single reaction and a single nuclear target.
Each new probe and target combination historically required a bespoke code that recomputes
the nuclear-structure ingredients alongside the reaction-specific
amplitude. The TDA factorization removes this redundancy: nuclear
structure enters only through precomputed one-body and two-body
densities, and the same densities serve every reaction. This division
also matches how few-body theory groups divide labor, with one team
producing densities from a chiral nucleon-nucleon interaction and
another producing reaction kernels from $\chi$EFT or partial-wave
analyses.

`DensityScattering` provides the missing reusable infrastructure:
density input/output, quantum-number summation, angular and radial
quadrature, and the convolution itself. Implementing a new reaction
reduces to writing a kernel subroutine of order 500 lines of Fortran;
all remaining machinery is inherited. The framework supports the
calculations of the corresponding PhD thesis on Compton scattering,
neutral pion photoproduction, and elastic pion scattering on
${}^3$H, ${}^3$He, ${}^4$He, and ${}^6$Li. The modularity also supports targeted
sensitivity studies, since a single kernel can be evaluated against
densities generated from different chiral orders, cutoffs, or
similarity renormalization group resolution scales [@Reinert2018]
without rebuilding either ingredient.

# State of the field

Few publicly available codes target the same problem. The original
TDA application to $\gamma\,{}^3$He [@Hammer3He] and the subsequent
extension to ${}^4$He [@hammer4He] used in-house implementations not
released as community software. Lenkewitz and collaborators
[@Lenke2011; @Lenke2013; @LenkeThesis] and Braun [@BraunThesis]
calculated threshold pion photoproduction on light nuclei using
private codes. Liebig and collaborators [@Liebig_2011] provided the
diagrammatic basis for elastic pion-nucleus scattering, again with
private code. `DensityScattering` collects, generalizes, and
publicly releases the infrastructure required to convolve any of
these kernels with TDAs computed externally, and adds new kernels for
Compton scattering on ${}^6$Li, threshold and above-threshold pion
photoproduction, and elastic pion scattering. The transition
densities themselves are obtained from the `nucdens` package
[@nucdens] of A. Nogga (Forschungszentrum Julich), with which
`DensityScattering` integrates directly.

# Software design

The architecture separates process-independent infrastructure (the
*mantle*) from process-specific physics (the *kernel*). The mantle
reads input parameters, loads the densities from a compressed HDF5
[@hdf5; @zfp; @h5z_zfp] archive, performs quantum-number summations,
carries out the convolution, and writes the resulting nuclear
amplitude. The kernel returns the single-nucleon scattering matrix
(one-body case) or the two-nucleon kernel matrix (two-body case)
evaluated at a momentum point and a set of quantum numbers.

For the one-body sector, the convolution accumulates
$$
  \mathcal{M}_{1N}(e, M', M)
  \mathrel{+}= A \,
  \rho^{(1N)}_{\alpha,\alpha'}\,
  \mathcal{M}^e_N\!\bigl(m_3^{s\prime}, m_3^{s}, m_3^{t}, m_3^{t\prime}\bigr),
$$
where $e$ labels external quantum numbers (e.g., photon polarization
or pion charge), $M, M'$ are nuclear spin projections,
$m_3^{s,t}$ and $m_3^{s\prime,t\prime}$ the single-nucleon spin and
isospin projections, $\rho^{(1N)}$ the one-body density, and
$A = \binom{A}{1}$ counts equivalent active nucleons. The two-body
sector additionally integrates over the relative momenta of the
$(12)$ subsystem and sums over its orbital, spin, and isospin
quantum numbers.

Two-body kernels with pion-pole intermediate propagators of the form
$1/\vec{q}^{\,2}$ exhibit a removable moving singularity in the
integration domain. Prior calculations
[@Lenke2013; @BraunThesis] regularize this with an extrapolation
$\Lambda \to 0$ in a Yukawa regulator. `DensityScattering`
implements an alternative variable substitution that maps the
singularity to the origin of the integration variable, where it is
canceled by the Jacobian of the spherical measure. This avoids the
extrapolation step and the associated systematic uncertainty.
Two-body kernel directories are prefixed `varsub-` to indicate this
treatment.

The repository contains three process-independent components:
`common-densities/` (shared modules), `onebodyvia1Ndensity/`
(one-body mantle), and `varsub-twobodyvia2Ndensity/` (two-body
mantle). Reaction kernels currently shipped are
`PionPhotoProd.onebody/`, `PionPhotoProdThresh.onebody/`,
`PionPion.onebody/`, `varsub-PionPhotoProdThresh.twobody/`, and
`varsub-PionPion.twobody/`. A `tools/` directory provides Python
helpers for downloading densities, running parameter scans, and
analyzing output; these helpers also handle the conversion between
photon energy plus scattering angle and the kinematic variables of
each reaction. The code is compiled with `gfortran` [@gfortran] and
links against parallel HDF5 [@hdf5] with the H5Z-ZFP filter plugin
[@h5z_zfp; @zfp] for reading the compressed density archives.

# Mathematics

The full elastic amplitude is the sum of one-body and two-body pieces,
$\mathcal{M}_X = \mathcal{M}_{1N} + \mathcal{M}_{2N}$, with
differential cross section
$\mathrm{d}\sigma/\mathrm{d}\Omega = (64\pi^2 s)^{-1}\,|\mathcal{M}_X|^2$,
where $s$ is the Mandelstam invariant. The two-body convolution
integrates the kernel
$\mathcal{K}_{2N}(d, e; s_{12}', m_{s_{12}}', s_{12}, m_{s_{12}})$
against the corresponding two-body density, summed over the diagram
index $d$, the external channel $e$, and the $(12)$-subsystem
quantum numbers, including the orbital, spin, and isospin projections
constrained by the antisymmetry condition
$(-1)^{s_{12}+l_{12}+t_{12}} = -1$.

# Research impact

`DensityScattering` underpins the first calculation of Compton
scattering on ${}^6$Li in $\chi$EFT and the associated PhD thesis
at The George Washington University. It is the
computational basis for the calculations supporting the HI$\gamma$S
proposal [@Godagama2025] for measuring Compton scattering off
${}^6$Li. The framework reproduces the threshold pion photoproduction
amplitudes of Lenkewitz and collaborators [@Lenke2011; @Lenke2013] for
$^{3}$H and $^{3}$He, and resolves a $\sqrt{2}$ normalization
discrepancy in the $^{6}$Li results of Braun [@BraunThesis] arising
from the convention for spin-1 angular-momentum matrices. Cutoff
selection follows the Bayesian convergence analysis of
[@Millican2024]. Above-threshold one-body amplitudes are constructed
from the SAID partial-wave analyses [@SAID_website; @SAID2023piphoto].
The package is provided with a citable Zenodo archive [@long_2026].

# Related publications

The accompanying PhD thesis (in preparation, The George Washington
University, 2026) contains the full set of results for Compton
scattering, threshold and above-threshold neutral pion
photoproduction, and elastic pion-nucleus scattering on ${}^3$H,
${}^3$He, ${}^4$He, and ${}^6$Li. No journal-version manuscript is
currently in review.

# AI usage disclosure

Anthropic's Claude (Opus and Sonnet model families, accessed through
the Claude Code CLI between 2024 and 2026) was used to assist with
the following: drafting and copy-editing of this manuscript; debugging
of Fortran and Python source code in the `DensityScattering`
repository; and translating Python prototypes into Fortran for the
above-threshold one-body kernels (`PionPhotoProd.onebody/` and
`PionPion.onebody/`), which evaluate single-nucleon amplitudes from
SAID partial-wave inputs. The author reviewed, edited, tested, and
validated all AI-assisted output. All scientific content, design
decisions, kernel derivations, and physics validation are the
author's own. AI tools were not used in interactions with editors or
reviewers.

# Conflict of interest

The author declares no conflict of interest.

# Acknowledgements

The author thanks Harald W. Grieβhammer for supervision of the
underlying thesis work and for guidance on the TDA formalism.
The author also thanks Andreas Nogga and Xiang-Xiang Sun
(Forschungszentrum Ju}lich) for providing the
nuclear transition densities through the `nucdens` package
[@nucdens]. Financial support for this work was provided by [FUNDING SOURCES TO BE
COMPLETED].

# References
