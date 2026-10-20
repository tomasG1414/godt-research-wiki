---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
source: "https://doi.org/10.1038/s41598-019-55330-4"
---

# 3D-Printed Biological Cell Phantom for Testing 3D Quantitative Phase Imaging Systems

**Citation:** M. Ziemczonok, A. Kuś, P. Wasylczyk, M. Kujawińska, "3D-printed biological cell phantom for testing 3D quantitative phase imaging systems," Scientific Reports 9, 18872 (2019). Warsaw University of Technology / University of Warsaw.

**Gist:** Fabrication and verification paper for the calibrated [[nano-printed-cell-phantom]] used experimentally by [[2411-08423-godt]]. Uses two-photon laser photolithography (direct laser writing) with refractive-index engineering (via laser power / exposure dose) in a single fabrication step to build a life-scale phantom mimicking real mammalian-cell sub-structures with precisely known, independently-verified geometry and RI — filling a gap where prior QPI test targets were either purely numerical or, experimentally, unknown/simple uniform microspheres.

## Key results

- Phantom design (Fig. 1): truncated-ellipsoid cell body, 30×25 μm lateral × 11 μm height, containing (1) lateral (XY) resolution test lines at 300/500/700 nm widths, ΔRI=0.03; (2) three spherical "nucleoli" (2.6 μm diameter) suspended in a lower-RI "nucleus" region, ΔRI=0.015; (3) vertical (Z) resolution test lines at 800/1100/1400 nm, ΔRI=0.03; (4) a cylinder with a slow radial RI gradient (ΔRI=0.015) mimicking natural within-cell RI variation, useful for testing whether a reconstruction algorithm's piecewise-constant assumptions hold; (5) a thin wall separating the low-RI nucleus from the surrounding medium, testing boundary-retrieval accuracy.
- Achievable RI contrast via exposure dose alone: high dose → n=1.5167±0.0008, low dose → n=1.5021±0.0012, unexposed → n=1.4882±0.0006 (max ΔRI ≈ 2×10⁻² via cured-polymer dose, extendable to 4×10⁻² using uncured liquid regions).
- Independently verified via SEM (geometry, <5% relative error on resolution-line widths) and white-light interferometry (height, ±200nm) and cross-checked against isotropic (rotation + illumination-scanning) holographic tomography for RI (accuracy ~10⁻⁴) — a three-modality verification loop giving a trustworthy ground truth distinct from the phantom's own design file.
- Applied as a benchmark for a limited-angle holographic tomography (LAHT) system (direct Fourier-diffraction-theorem/Rytov reconstruction + iterative Gerchberg-Papoulis + [[krauze2016|GTVIC]]-style TV-constrained 3D mask — same reconstruction lineage as [[krauze2016]]): external boundaries recovered accurately (height error 200nm); resolves 300nm features at 600nm period laterally, 800nm at 2500nm period axially, but with degraded contrast for the smallest features (P-V ratio badly reduced, e.g. lateral-300nm P-V drops to 0.17 of model in experiment) — attributed to insufficient sampling, iterative-replenishment behavior, or high-RI-gradient regions violating the weak-scattering assumption.
- Nucleoli (spherical, known ground truth) reconstructed with accurate lateral diameter (~130nm error) but ~20% vertical elongation — the paper explicitly attributes this to "elongation artifacts... typical in the limited angle version of holographic tomography," i.e. exactly the [[missing-cone-problem]] signature. Measured RI of the spheres was underestimated experimentally (8% cumulative RI loss vs. <0.5% in the numerically-simulated version of the same reconstruction pipeline) — separates hardware/experimental error from pure reconstruction-algorithm error, a useful methodological pattern (compare model → synthetic reconstruction → experimental reconstruction, three-way).

## Method

Two-photon photolithography (Nanoscribe Photonic Professional GT, IP-L 780 photoresist, 100× NA 1.4 objective); RI controlled solely by laser power at fixed scan parameters; phantoms printed either at an optical fiber tip (for rotation-based isotropic reference tomography) or on a coverslip (for illumination-scanning LAHT testing) to keep verification and target-system measurements compatible. Metrology combines SEM (lateral/geometry), white-light interferometry (height), and isotropic holographic tomography (RI ground truth), then compares model vs. numerically-simulated reconstruction vs. experimental reconstruction to separate algorithmic from hardware error sources.

## Relevance

Per [[godt]]: the fabrication and verification source for the [[nano-printed-cell-phantom]] that [[2411-08423-godt]] uses for its experimental validation (its total-shear vs. gradient-mode comparison in Fig. 3).

## Raises

- The three-way model/synthetic/experimental comparison methodology here (isolating hardware error from algorithm error) is a good pattern for the multi-shear simulation thread.
