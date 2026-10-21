---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-10-20
status: current
source: "https://doi.org/10.3390/s24237714"
---

# Aspherical Surface Wavefront Testing Based on Multi-Directional Orthogonal Lateral Shearing Interferometry

**Citation:** Yahui Zhu, Ailing Tian, Hongjun Wang, Bingcai Liu, *Sensors* **24**(23), 7714 (2024). DOI `10.3390/s24237714`. Huaiyin Normal University + Xi'an Technological University.

**Gist:** Optical metrology, not microscopy — testing aspheric mirrors, not cells. It earns its place here because it is the corpus's only source that actually *measures* what extra shear **directions** buy you, and the honest answer is: **almost nothing beyond noise averaging**. That is a genuinely useful negative result for [[multi-shear-acquisition]], and it points the project's argument away from direction count, towards what directions beyond one orthogonal pair cannot add. **[audit 2026-10-20]** Superseded in 3D: a second direction also covers most of the first direction's null planes (0.8 x + 0.8 y loses 0.09 % of the support vs 3.02 % for 0.8 x); see [[multi-shear-acquisition]].

## Key results

**Four sheared replicas in one exposure, from two stacked birefringent crystals.** Two calcite polarization beam displacers (PBD1 fast axis along $x$; PBD2 rotated 45°, 9.28 mm thick) form a "double birefringent crystal beam displacer" (DBCs-BD) producing $E_{x\pm}$ and $E_{y\pm}$ simultaneously. A 45°/135° analyser makes them interfere, and one CCD frame carries **both orthogonal shear directions at once**. Rotating the whole module gives further direction pairs (0°/90°, 15°/105°, 45°/135°, 60°/150°). Shear $s$ = 1 mm = **1/10 of the beam aperture** — very large in fractional terms compared with microscopy shears.

**Extra directions barely move the result.** This is the finding worth carrying away. Experimental reconstruction of a 4 μm-asphericity paraboloid (Table 2, $\lambda$ = 632.8 nm):

| Directions | PV (λ) | RMS (λ) |
|---|---|---|
| 0° and 45° | 7.9848 | 2.3744 |
| 0° and 120° | 7.7962 | 2.3748 |
| 45° and 120° | 7.9856 | 2.3746 |
| 0°, 45°, 120° | 7.9834 | 2.3721 |

Deviation from design value likewise moves from PV 0.6075λ / RMS 0.0543λ (one direction pair) to PV 0.6065λ / RMS 0.0539λ (three) — a 0.7 % change in RMS for three times the data. The 10.3185 μm hyperboloid behaves the same way. In simulation the gain is real but modest: residual ΔRMS improves from $1.68\times10^{-5}\lambda$ (0°/90°) to $9.90\times10^{-6}\lambda$ (0°, 90°, 45°, 135°), a factor of 1.7 — consistent with averaging independent random errors, not with recovering missing information.

**The authors' own explanation matches that reading**: averaging differential wavefronts and then averaging Zernike coefficients "effectively reduc[es] the impact of system random errors."

**Validation is solid.** Against a LuphoScan profilometer: RMS agreement 0.0085λ (4 μm part) and 0.0078λ (10.3 μm part), i.e. better than λ/100; better than λ/180 after adding integral-minimisation fitting. PBD orthogonality error within ±0.5° costs at most $1.4\times10^{-4}\lambda$ RMS.

## Method

Modal, not spectral. The wavefront is expanded in the first $N$ Zernike polynomials; each shear direction contributes a *differential* Zernike design matrix $\Delta Z$; coefficients come from ordinary least squares $a = (\Delta Z^T\Delta Z)^{-1}\Delta Z^T\Delta W$ (Eq. 21), stacking direction pairs as extra rows, then averaging coefficients across direction sets (Eq. 22). Differential phase per direction is extracted by Fourier-transform demodulation of a single carrier-fringe interferogram. Because the aspheric test is non-null, ray tracing in ZEMAX plus iterative reverse optimisation removes retrace error, and an "optimal spherical surface" reference is subtracted analytically.

## Relevance

Per [[godt]] and [[multi-shear-acquisition]]:

- **The negative result clarifies where the gain is.** Direction count alone gives small gains — at least in 2D wavefront metrology, once you already have one orthogonal pair. That is consistent with the theory: an orthogonal pair already spans the full transverse gradient, so further directions are largely redundant, and this hardware delivers a full pair in *every* exposure. It leaves **magnitude diversity** — which fills the transfer nulls of [[finite-shear-null-experiment]] — as the load-bearing part of the multi-shear argument, and direction diversity as the fix for the narrower $K_x=0$ blind line. The two motivations should not be presented as equally strong. **[audit 2026-10-20]** Superseded in 3D: a second direction also covers most of the first direction's null planes (0.8 x + 0.8 y loses 0.09 % of the support vs 3.02 % for 0.8 x); see [[multi-shear-acquisition]].
- **Caveat before over-reading it.** This is a smooth, low-order, modally-representable wavefront fitted in a Zernike basis. A global modal fit is largely *immune* to transfer nulls, because a null at one spatial frequency is compensated by the basis's global support — so this experiment could not have revealed null-filling gains even if they existed. It is evidence about direction redundancy, not evidence against magnitude diversity. *(inference)*
- **Modal versus spectral fusion is a real design choice for the extension.** This paper fuses shear channels by stacking rows in one least-squares Zernike problem; [[2509-15382-osi-flex]] fuses them by summing data-fidelity terms over channels in an autodiff optimisation; GODT works in $\tilde O$-space. The modal route trades high-frequency capability (capped by $N$) for null-immunity and heavy regularisation — probably wrong for tomography of textured biological samples, but it is the third option on the table and should be named as such.
- **Hardware note.** The DBCs-BD is a fourth route to multi-direction shear, alongside [[nguyen2025-deep-learning-lsi]]'s polarization gratings and [[sciadv-adl0501]]'s metalens: passive, cheap, single-exposure, no diffractive elements. It needs a spatially coherent source (He-Ne here) for Fourier-transform carrier demodulation, so as written it is **incompatible with GODT's LED** for the same reason the carrier-multiplexing route already noted in [[multi-shear-acquisition]] is. The crystals themselves would work under low coherence; the demodulation scheme would not.

## Raises

- Zhu et al., *Appl. Opt.* **62**, 654–664 (2023) — QWLSI based on double birefringent crystals; the source of Eqs. 1–4 here and the primary hardware paper.
- Zhu, Tian, Wang & Liu, *Appl. Phys. B* **130**, 218 (2024) — "Wavefront Reconstruction Based on Multi-directional Orthogonal Lateral Shearing Interferometry"; the companion reconstruction paper, likely with the derivation this one compresses.
- Shanmugam, Light, Turley & Falaggis, *Light: Advanced Manufacturing* **3**, 193–210 (2022), "Variable shearing holography with applications to phase imaging and metrology" — a **variable**-shear system, i.e. magnitude diversity in hardware.
- Does the direction-redundancy result survive when the reconstruction is spectral rather than modal, and when the object is high-frequency? Testable in `lab/` with the existing null-experiment machinery.
