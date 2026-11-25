---
type: concept
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# Convergent Born series

An iterative solver of the full inhomogeneous Helmholtz equation $\nabla^2\psi+k(\mathbf r)^2\psi=-S$
that, unlike the ordinary Born series, converges for objects of any size and contrast.
Introduced by Osnabrugge, Leedumrongwatthanakun & Vellekoop, *J. Comput. Phys.* **322**,
113 (2016), arXiv:1601.05997 *(not ingested; the formulas below were checked against the
arXiv text)*. This project uses it as an **exact referee**: the one forward model here
that makes no weak-scattering, forward-only, paraxial or thin-slice approximation
(`lab/shearlib.py`, `CBS`).

## The iteration

Split $k(\mathbf r)^2=k_0^2+i\varepsilon+V(\mathbf r)$, choose $k_0^2$ at the midpoint of the
range of $\operatorname{Re}k^2$ and $\varepsilon\ge\max|k^2-k_0^2|$, and precondition with
$\gamma=(i/\varepsilon)V$:

$$\psi\leftarrow\psi+\gamma\big[g_0*(V\psi+S)-\psi\big],\qquad \tilde g_0(\mathbf p)=\frac{1}{|\mathbf p|^2-k_0^2-i\varepsilon}.$$

Each step costs one forward and one inverse 3D FFT, and advances the field by about
$2k_0/\varepsilon$, so the cost grows with the contrast and with the absorbers'
strength. Here the unknown is the *scattered* field of a plane wave, with source
$S=-O\psi_i$, so the incident wave is never absorbed.

## Setting it up for transmission ODT

Three things learned the hard way, each recorded in `lab/README.md`:

- **A periodic transverse box does not converge in any useful time.** Without absorption
  at the sides, near-grazing waves and guided modes of the periodic array of objects never
  reach the axial absorbing layers. The update shrank only algebraically (still $10^{-4}$
  after 3000 iterations on a small sphere). An absorbing frame acting on the scattered
  field restores geometric convergence: about 200 iterations to $10^{-7}$.
- **Add the frame and the axial layers where they meet, do not take the larger.** Taking
  the maximum lowers $\varepsilon$ but left a weakly damped mode in small boxes (1500
  iterations, unconverged, against about 200).
- **The window decides what is measured.** At 60° illumination the forward lobe of a
  6 μm cell walks about 11 μm sideways. Measured on a central 12.8 μm crop, even the exact
  solution's linear part is off first Born by tens of percent. Measured on the whole
  window at the reconstruction lattice (`shearlib.lattice_samples`) it converges on
  Born's as the window grows ([[full-wave-referee]] §3).

## Validation

Against the exact scalar partial-wave (Mie) series for a sphere of radius 1.5 μm
([[full-wave-referee]] §1): 0.4 % relative error in the far-field amplitude at normal
incidence and 0.9 % at 60°, for spheres where first Born is off by 16 % ($n=1.34$) and
80 % ($n=1.38$, phase delay 1.9 rad). The forward amplitude agrees to 0.1–0.2 % (normal
incidence) and 0.5–0.6 % (60°).

## What it costs

On a $512\times512\times132$ grid in single precision, one illumination took about seven
minutes on four cores (about 350 iterations, while sharing the machine with another job). That is fine for a referee on one illumination, and too slow for
the 180 illuminations of a design study. [[gridding-and-forward-model]] therefore uses
[[multi-slice-forward-models|multi-slice Born]], after checking it against this solver.
