# Index

Routing file: link + one-line summary per page, nothing more. Read this first on
any query. Updated on every ingest.

On GitHub the double-bracket names below are plain text: each one is a file of
that name under `wiki/` (for example `wiki/notes/finite-shear-null-experiment.md`).

## Projects

- [[godt]] — Gradient Optical Diffraction Tomography: core paper *Communications Physics* 8, 510 (2025); multi-shear simulation study

## Sources

- [[s42005-025-02416-3]] — **canonical GODT reference** (Comms Phys 2025): published version, adds finite-shear theory + shear-value study
- [[2411-08423-godt]] — GODT arXiv preprint, superseded by the published version; retained to document what changed
- [[vishnyakov2016]] — DIC tomography (2016): reconstructed the 3D RI derivative before GODT, but via Radon/straight-ray rather than diffraction tomography
- [[wolf1969]] — origin of the ODT projection theorem (Born approximation, Fourier diffraction theorem)
- [[chen1998]] — quantitative Born-vs-Rytov validity criteria: Rytov's range is 3x Born's, but fails via phase unwrapping
- [[pan1983]] — direct inversion (bilinear interpolation) matches filtered backpropagation at 1/14th the cost
- [[krauze2016]] — GTVIC regularization for limited-angle ODT; quantifies the phase-only vs. complex-data reconstruction quality gap
- [[ziemczonok2019]] — fabrication/verification of the calibrated 3D-printed cell phantom used in GODT's experimental validation
- [[brenner1993]] — wave-propagation method (WPM): angle-unrestricted forward light-scattering simulation, GODT's forward-model engine
- [[lim2015]] — GP vs. EP vs. TV regularization for missing-cone mitigation: no universal winner, depends on sample structure
- [[2024-06-14-598984-polarization-hilbert-holotomography]] — same PG shearing hardware as GODT run in total-shear mode; single-shot Hilbert-transform phase demodulation; names the EWALD reconstruction codebase
- [[optica-8-1-6]] — qOBM extended to 3D RI tomography of opaque samples in epi-mode via a non-paraxial 3D OTF; shares the missing-cone limitation of limited-angle tomography
- [[nguyen2017]] — GLIM: closest predecessor to GODT; same gradient measurement, but integrates + z-scans instead of tomographic inversion
- [[kandel2019]] — epi-GLIM for opaque samples; independently confirms shear makes the transfer function non-symmetric, and sidesteps it with fiducial-PSF deconvolution
- [[kwon2019]] — 1 mm³ metasurface phase-gradient microscope; single-shot via spatial multiplexing of three DIC phase offsets
- [[sciadv-adl0501]] — single metalens captures both orthogonal shear directions and four phase shifts in one exposure
- [[nguyen2025-deep-learning-lsi]] — two-direction shear built from polarization gratings; neural demultiplexing of a composite interferogram
- [[2509-15382-osi-flex]] — GODT authors' open-source multi-shear phase-reconstruction framework; joint phase + shear estimation via autodiff
- [[ma2017]] — WPM as a nonlinear forward model inside tomographic reconstruction, for multiply-scattering objects
- [[guizar-sicairos2011]] — X-ray phase tomography: direct reconstruction from phase derivative, wrapping-insensitive derivative estimation, fiducial-free alignment
- [[mehta-sheppard2008]] — DIC partially coherent image formation: the origin of the finite-shear filter and its null comb, plus the shear-vs-coherence coupling
- [[ledwig-robles2023]] — broadband partially coherent 3D OTFs for arbitrary sources; Ewald *plane* replaces Ewald sphere, $O(n^3)$
- [[chung2021-projectiongan]] — unsupervised CycleGAN missing-cone filling in projection space; 27 s vs TV's 605 s, but fills by generation not measurement
- [[zhu2024-appl-phys-b]] — multi-direction shearing reconstruction via differential Zernike fitting at arbitrary angles; where direction gains saturate
- [[zhu2024-sensors]] — the same hardware applied to aspheric testing: three shear directions buy ~0.7 % RMS
- [[zhu2025-sci-rep]] — single-frame four-step polarization phase shifting; 4→16 shear directions changes RMS by 0.0000 λ
- [[2605-04297-double-sagnac-ls-qpm]] — cascaded Sagnacs give both shear directions in one shot; reaches 3D by z-scanning rather than diffraction tomography
- [[servin2007]] — 2D LSI frequency response, inverse-filter poles at the shear nulls, regularised LS inverse, spatial support analysis, $n$-direction generalisation; establishes null comb as prior art

## Concepts

- [[optical-diffraction-tomography]] — Wolf's projection theorem, Ewald sphere, direct inversion / filtered backpropagation
- [[rytov-approximation]] — first-order weak-scattering linearization behind ODT/GODT reconstruction
- [[quantitative-phase-imaging]] — 2D phase-delay contrast, coherence tradeoffs, predecessor techniques
- [[common-path-shearing-interferometry]] — self-referencing shearing holography, total-shear vs. gradient regimes
- [[gradient-phase-imaging]] — DIC/qDPC/GLIM/GROM/QWLI/metasurface family GODT extends into tomography
- [[missing-cone-problem]] — NA-bounded illumination angle range leaves a cone of 3D spectrum unfilled
- [[multi-shear-acquisition]] — the project's core extension: why multi-shear, acquisition strategies, hardware feasibility, what's still open
- [[optimal-experimental-design]] — convex A/D/maximin design over time shares, the equivalence-theorem certificate, and why optimal designs here are not sparse
- [[partially-coherent-transfer-function]] — TCC/Hopkins theory, coherence ratio, WOTF/PGTF, and the broadband Ewald-plane extension
- [[pseudo-thermal-illumination]] — rotating-diffuser source: kills spatial coherence only, which reopens single-shot carrier multiplexing for multi-shear

## Entities

- [[shepp-logan-phantom]] — standard synthetic tomography test object
- [[nano-printed-cell-phantom]] — calibrated 3D-printed phantom used in GODT experimental validation

## Notes

- [[forward-citation-sweep-2026-07]] — OpenAlex forward sweep over 9 anchors: 2024-26 leads, two preprint-to-published corrections, three closely related papers from the GODT authors' collaboration
- [[finite-shear-null-experiment]] — confirms GODT's finite shear nulls 3 in-band frequencies at the published 0.8 μm (simulated, monochromatic); multi-shear fills them if ratios avoid small integers
- [[multishear-fusion-design-problem]] — concept glossary + precise formulation: per-voxel Fisher information under Ewald sampling, BLUE/Wiener fusion, convex optimal shear design
- [[multishear-fusion-design-solution]] — results: J separates, 3.02 % true lost volume (0.09 % with a second direction), rule Δ ≈ 2π/K_max, Õ-space fusion beats equal-weight averaging 21×, sparsity conjecture refuted (audited 2026-09-25, 2026-10-20)
- [[multishear-fusion-derivation]] — the theory from Helmholtz up: Fourier diffraction theorem, separation, closed-form M(K), nulls, fusion = precision-weighted averaging, OED certificate, E1/E2
- [[multiple-scattering-check]] — WPM re-test of the design solution: null planes stay put, margins shrink with contrast, gridding is a large error term, shear data skip the unwrapping step
- [[colour-multiplexed-shear]] — colour as a free shear channel (Δ ∝ λ): equals a mechanical shear pair, blue + green is the right pair, unnecessary at the design-rule shear; related 2D hardware in a 2026 preprint by the GODT authors' group
