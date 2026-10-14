---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-07-22
status: current
source: "https://arxiv.org/abs/2605.04297"
---

# High-Sensitivity, High-Throughput Double Sagnac Lateral Shearing Quantitative Phase Microscopy and Tomography with Pseudo-Thermal Illumination

**Citation:** Paweł Gocłowski, Maciej Trusiak, Balpreet S. Ahluwalia, Azeem Ahmad. arXiv:2605.04297v1, 2026-05. UiT The Arctic University of Norway + Warsaw University of Technology. (18 pp. incl. supplement).

**Gist:** Two Sagnac shearing interferometers cascaded in orthogonal orientations, bolted onto the output port of a standard bright-field microscope, produce four laterally sheared wavefronts whose mesh-like interferogram encodes $\partial\phi/\partial x$ and $\partial\phi/\partial y$ simultaneously in **one shot**. Pseudo-thermal (dynamic-speckle) illumination gives low spatial coherence for speckle suppression while keeping temporal coherence high enough for dense carrier fringes. Fourier-transform demodulation yields two gradient maps, integrated in 2D to a quantitative phase map. Validated on polystyrene beads, fixed and live HeLa cells, and a z-scanned bead-in-PDMS volume.

## Key results

**Sensitivity and speed.** Spatial phase sensitivity 11.7 mrad ($x$) and 13.3 mrad ($y$); temporal 8.7 mrad ($x$) and 9.7 mrad ($y$), against the ~30 mrad the authors quote as typical for non-common-path systems (§3.1, Fig. 2, from 1200 sample-free frames at 20 Hz over 1 min). Single-shot operation means frame rate is the only speed limit; they claim ~1 kHz is reachable with a faster camera.

**Resolution.** With a 60×/1.2 NA objective, the integrated phase map resolves ~690 nm along $x$ and ~776 nm along $y$, versus a theoretical limit near 500 nm (Supplement §3). The anisotropy is attributed to unequal shear in the two axes. The stated selling point against QWLSI is avoiding its grating-pitch resolution penalty, which the introduction puts at "up to fourfold" relative to the diffraction limit.

**Systematic shear-magnitude study (§3.2, Supplement §1, §3).** The most directly useful part of the paper for this project. Shear is tuned by translating the camera on a 12 mm stage; 25 positions were scanned at 0.5 mm steps with 1 µm and 300 nm beads, then refined at 0.1 mm steps over 4.5–5.5 mm with a phase USAF target. The behaviour is a clean optimum:

- **Too much shear** (positions 0–4.0 mm): beads distorted, individual beads unresolvable. Attributed to shear exceeding "the effective diffraction-limited sampling of the system," i.e. resolution loss.
- **Optimum** (~4.8–5.0 mm): beads well resolved with defined boundaries.
- **Too little shear** (5.5–6.0 mm): phase contrast collapses toward the zero-shear condition, where the sheared copies cancel and no phase is recoverable.

**Shear is neither isotropic nor spatially uniform.** The zero-shear camera position differs between axes — ~5.7 mm for $x$, ~6.4 mm for $y$ — and the integrated map's apparent zero sits near 6.0 mm, roughly their average. Correct reconstruction required integrating with a $dx\!:\!dy$ ratio of **1.5:1** rather than 1:1; equal values elongated the reconstructed beads (Supplement §2, Fig. S2). Shear also varies *across the field of view* from tube-lens distortion, degrading the right side of the image, and this — not raw resolution — is why the final operating point was chosen conservatively.

**Volumetric imaging (§3.6).** 1 µm beads in a ~100 µm PDMS layer, z-scanned at 0.5 µm steps. Axial resolution ~4.5 µm (FWHM ≈ 9 px on a bead's $z$ profile).

**Biology.** Fixed HeLa (subcellular structure in both gradient and integrated maps) and live HeLa over 2 h at 1 min intervals at room temperature, capturing membrane contraction and intracellular motion. Measured differential phase stays within ±0.3 rad, well inside ±π, so phase unwrapping is largely unnecessary.

## Method

Four-beam interference (§2.1). Two beams sheared by $\pm d_1/2$ along $x$ carrying tilt $\pm\theta_x$, two by $\pm d_2/2$ along $y$ with tilt $\pm\theta_y$. Keeping only the $(E_1,E_2)$ and $(E_3,E_4)$ cross terms and Taylor-expanding for small shear gives

$$I \approx I_0\left[1 + b_x\cos\left(2\pi f_x x + d_1\frac{\partial\phi}{\partial x}\right) + b_y\cos\left(2\pi f_y y + d_2\frac{\partial\phi}{\partial y}\right)\right]$$

— two orthogonal carriers, each phase-modulated by one gradient component, in a single frame. Demodulation is the standard Takeda Fourier method (their ref. 22): isolate the two first-order peaks with circular masks, inverse transform, then integrate the two gradient maps in 2D using D'Errico's off-the-shelf MATLAB `inverse-integrated-gradient` package (ref. 24).

**Optical layout (§2.2, Fig. 1).** After the microscope's tube lens forms an image plane, a 45° polarizer feeds Sagnac 1 (PBS1 + mirrors M3/M4), which produces two orthogonally polarized $x$-sheared copies. A second 45° polarizer (P2) equalizes their polarization so they can interfere, then the composite enters Sagnac 2 (PBS2 + M5/M6) oriented orthogonally, producing the $y$ shear. P3 before the camera brings all four beams to a common polarization. Shear magnitude and fringe spacing are set by mirror tilt (M3–M4 for $x$, M5–M6 for $y$) plus camera translation.

**Illumination.** Pseudo-thermal — laser through a rotating diffuser (RD in Fig. 1a) — see [[pseudo-thermal-illumination]]. This is the design choice that makes single-shot Fourier demodulation work under low spatial coherence, and it is a different tradeoff from [[s42005-025-02416-3|GODT]]'s LED.

**Prior version.** Their ref. 21 is the single-Sagnac predecessor: Gocłowski et al., "High-Fidelity Single-Shot Quantitative Differential Phase Microscopy Using Pseudothermal Sagnac Interferometer," arXiv:2604.20739 (2026).

## Relevance

### For [[godt]]: how its tomography differs

The paper reaches 3D by **z-scanning** rather than diffraction tomography. Section 3.6 stacks single-plane integrated phase maps acquired at 0.5 µm axial steps and renders the volume, with gradients integrated to 2D phase first. The authors name diffraction tomography as a natural extension: the approach "can be extended to optical diffraction tomography if multiple interferograms corresponding to different illumination directions are acquired and used for 3D refractive index reconstruction" (§3.6).

So this sits in the same position [[nguyen2017|GLIM]] does relative to GODT — a strong 2D gradient method reaching 3D by defocus stacking. The ~4.5 µm axial resolution is the signature of that: sectioning by depth of field, not by tomographic inversion.

Multi-shear *fusion in $\tilde O$-space*, the subject of this project's simulations, is therefore a different question from the one this paper answers.

### What it does contribute

1. **Two orthogonal shear directions, single shot, under low spatial coherence — a route [[multi-shear-acquisition]] had marked closed.** That page lists carrier-frequency multiplexing as incompatible with GODT because it needs spatial coherence. This paper shows the binding requirement is really *temporal* coherence for dense fringes, and that a pseudo-thermal source satisfies it while still suppressing speckle. The closure was over-broad; see [[pseudo-thermal-illumination]].
2. **Cascading two shearing interferometers works.** [[multi-shear-acquisition]] asked whether two transmissive PG pairs at different orientations can cascade without the first corrupting the second. This is the Sagnac analogue and the answer is yes, with a 45° polarizer between stages to reset the polarization basis.
3. **An experimental shear-magnitude optimum**, the empirical counterpart to [[finite-shear-null-experiment]]'s theoretical treatment — with the caveat below.
4. **Anisotropic and spatially varying shear is a first-order practical problem**, not a detail. The 1.5:1 $dx\!:\!dy$ correction and the FOV-dependent degradation are concrete evidence for the calibration warnings already on [[multi-shear-acquisition]].

### Two notes

**The shear study reports no transfer-null structure.** Degradation at large shear is attributed entirely to resolution loss, and at small shear to vanishing contrast. Nothing resembling the $K = 2\pi m/\Delta$ null comb of [[finite-shear-null-experiment]] is reported. This is not evidence against the nulls — the study reads bead sphericity and USAF element visibility, which are not sensitive probes of a zero at a single spatial frequency, and the scanned range is narrow — but it is the closest experimental shear-magnitude scan in the corpus, so it is worth noting.

**Operating point.** §3.3 gives a camera position of 5.0 mm as the optimal shear condition; Supplement §3 gives 4.8 mm. Worth confirming which applies before reusing a number quantitatively.

## Raises

- **Beam geometry of the cascade.** A serial cascade of two orthogonal Sagnacs could produce beams sheared in both axes, at $(\pm d_1/2, \pm d_2/2)$, which would add diagonal cross-terms to the spectrum alongside the two gradient peaks used *(inference)*. Worth understanding before building on this architecture.
- **Read the single-Sagnac predecessor** (arXiv:2604.20739) for the line-artifact claim: §2.2 asserts that integrating a single-direction differential phase "introduces line artifacts," an independent in-group argument for direction diversity arising from the *integration* step rather than the transfer function.
- **Patent application 2513325.7** covers LS-QPM (Ahmad and Ahluwalia). Relevant to any work built on this architecture.
- The collaboration now has three related shearing lines (GODT gradient tomography, LS-ODT partial shear, double-Sagnac LS-QPM); how they relate is worth following.
