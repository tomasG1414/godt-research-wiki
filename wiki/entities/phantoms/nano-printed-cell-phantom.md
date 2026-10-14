---
type: entity
projects: [godt]
created: 2026-07-18
updated: 2026-07-18
status: current
---

# Nano-Printed Cell Phantom

Calibrated 3D-printed phantom fabricated via two-photon laser photolithography (direct laser writing, RI controlled by exposure dose in a single fabrication step) — see [[ziemczonok2019]] for fabrication and independent verification (SEM, white-light interferometry, isotropic holographic tomography). Truncated-ellipsoid cell body (30×25×11 μm) containing lateral/vertical resolution test lines (300–1400nm, ΔRI=0.03), three spherical "nucleoli" (2.6 μm, ΔRI=0.015) in a low-RI "nucleus," a radial-RI-gradient cylinder (ΔRI=0.015, tests piecewise-constant reconstruction assumptions), and a thin nucleus-boundary wall.

Used by [[2411-08423-godt]] as the primary experimental validation target: imaged in both total-shear (ODT) and gradient (GODT) mode on the same polarization-grating common-path system, showing GODT gives higher contrast for subtle RI features at the cost of some residual diffraction-artifact fringes (phase-only/uniform-amplitude assumption). Independently, [[ziemczonok2019]] benchmarks a limited-angle holographic tomography system against this same phantom and finds the ground-truth-spherical nucleoli come out ~20% vertically elongated — a documented instance of the [[missing-cone-problem]] signature.

## Relevance

Per [[godt]]: a known-ground-truth physical phantom for validating reconstructions against real data.
