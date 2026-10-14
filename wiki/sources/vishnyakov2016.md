---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-07-22
status: current
source: "https://doi.org/10.1364/OL.41.003037"
---

# Differential Interference Contrast Tomography

**Citation:** G. Vishnyakov, G. Levin, V. Minaev, M. Latushko, N. Nekrasov, V. Pickalov, "Differential interference contrast tomography," Optics Letters 41(13), 3037–3040 (2016). All-Russian Research Institute for Optical and Physical Measurements (Moscow) / Khristianovich Institute (Novosibirsk). Cited by the published GODT paper ([[s42005-025-02416-3]], ref. 58).

**Gist:** Presents DICT (differential interference contrast tomography) — tomographic reconstruction whose output is a 3D DIC image, explicitly defined as the partial derivative of the 3D refractive-index distribution in the shear direction, $\text{DIC}(x,y,z) = \partial f(x,y,z)/\partial y$ where $f = n(x,y,z)-n_0$. Built on a common-path lateral-shearing phase-shifting interference microscope with low-coherence LED illumination and angular illumination scanning. **It reconstructs the same quantity as GODT, a tomographic 3D RI derivative**, but within a non-diffractive (straight-ray) framework.

## Key results

- **Core relation** (their Eq. 1-2): defines the 3D DIC tomogram as $\partial f/\partial y$, and shows that differentiating the Radon transform with respect to $y$ gives $\partial R(p,y,\theta)/\partial y = \iint \frac{\partial f}{\partial y}\delta(p-x\cos\theta-z\sin\theta)\,dx\,dz$ — i.e., **DIC projective data is itself a Radon transform of the derivative field**, so standard tomographic inversion applies directly to gradient data without integrating it first. Structurally the same insight as GODT's $\tilde{O}_\Delta = iK_x\tilde{O}$, but expressed in Radon rather than Fourier/Ewald terms.
- Because the DIC projective data is a valid Radon transform, iterative reconstruction algorithms that require this property (ART, Gerchberg-Papoulis) remain applicable — explicitly noted by the authors as the reason gradient data can be fed straight into standard tomographic solvers.
- **Crucially, DICT is non-diffractive.** The authors justify using a Radon-transform (straight-ray) algorithm rather than a diffraction-tomography algorithm on the grounds that their sample's RI contrast is small (Δn ≈ 0.02) so "diffraction at the edges of the cell is minor," citing Kujawinska et al. that filtered backprojection and filtered backpropagation give practically identical results at such contrast. **This is the decisive difference from GODT**, which builds on the first-order [[rytov-approximation|Rytov]] approximation and [[wolf1969|Wolf's]] Ewald-sphere theory and therefore does compensate diffraction.
- Experimental: point-source LED (650nm, 25nm bandwidth — low coherence, comparable to GODT's 530nm/35nm LED), two-axis galvanometer mirror for angular scanning over a 140° range, 40–125 projections, lateral shear 440 nm (small-shear/gradient regime, comparable to GODT's 0.33–0.8 μm effective shear), **five phase-shifted lateral shearing interferograms per DIC image** (same 5-frame phase-shifting count as GODT), least-squares phase retrieval, ART reconstruction with Gerchberg-Papoulis filtering. Demonstrated on a white blood cell (6 μm), 256³ voxels at 63 nm.
- Validation approach: compares reconstructed 3D DIC tomogram slices against *conventional* DIC images of the same sample at the same depths, showing the tomogram resolves nucleus and cytoplasm granules sharply where conventional out-of-focus DIC images are blurred — i.e., demonstrating the depth-sectioning benefit of tomographic reconstruction over z-scanned DIC.
- Explicitly acknowledges the [[missing-cone-problem]] (angular range limited to ≤140° by objective NA) and uses iterative ART + Gerchberg-Papoulis to mitigate — consistent with [[lim2015]].
- **Notable methodological remark**: the authors point out that the refractive index itself can be recovered by integrating the reconstructed 3D DIC tomogram along the shear axis, and argue this is more accurate than integrating the projection data before reconstruction (Vishnyakov et al. 2016) — i.e., if you want RI rather than its derivative, integrate *after* reconstruction, not before. GODT sidesteps this entirely by presenting the derivative itself as the contrast domain of interest.

## Method

Common-path lateral shearing Michelson interferometer in convergent light (beam splitter + two plane mirrors at the Fourier lens focal plane; shear produced by a small tilt of one mirror, shear direction set by tilt axis; no pinhole required, unlike competing common-path designs). Phase shifting via piezo-driven mirror. Transmission geometry, matched oil-immersion condenser and imaging objectives (100×, NA 1.30).

## Relevance

Per [[godt]]: the closest predecessor in the corpus. DICT reconstructed the 3D index derivative in 2016 from lateral-shearing, low-coherence gradient data under angular illumination scanning, an experimental configuration very similar to GODT's. The two differ in an important way:

- DICT reconstructs via **Radon transform / straight-ray tomography** and explicitly declines a diffraction algorithm because its particular sample had negligible diffraction. GODT derives the gradient-tomography relation within the **Rytov/Ewald diffraction** framework, so it compensates diffraction and is not restricted to the low-contrast, diffraction-negligible regime. GODT's Eq. 19 ($\tilde{O}_\Delta = iK_x\tilde{O}$) is a genuinely different and more general result than DICT's Radon-derivative relation.
- GODT demonstrates on dense/confluent samples (neurons) where total-shear common-path fails, and on a calibrated phantom ([[nano-printed-cell-phantom]]) with known ground truth; DICT demonstrates on a single white blood cell against conventional DIC images.

The published GODT paper draws exactly this distinction: DICT does not incorporate diffraction in its reconstruction, GODT does.


## Raises

- For GODT's actual samples, how much does diffraction compensation matter in practice? DICT's argument (Δn≈0.02 → diffraction negligible → straight-ray adequate) may or may not hold for neurons and the phantom. No ingested source quantifies this crossover.
