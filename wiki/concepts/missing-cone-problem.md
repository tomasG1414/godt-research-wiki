---
type: concept
projects: [godt]
created: 2026-07-18
updated: 2026-09-21
status: current
---

# Missing Cone / Missing Frequency Problem

In illumination-scanning [[optical-diffraction-tomography]], the achievable illumination directions are bounded by the numerical aperture (NA) acceptance angle, so the Ewald-sphere caps swept during acquisition never cover the full 3D object spectrum $\tilde{O}$ — a cone-shaped region around the optical (z) axis stays unfilled ("missing cone"). Consequence: **anisotropic, axially-limited 3D resolution** — reconstructions are systematically blurred/elongated along z relative to the transverse plane, regardless of reconstruction algorithm.

Traceable to [[wolf1969]]'s original projection theorem: only *homogeneous* (propagating) plane-wave components of the angular spectrum carry retrievable information ($p^2+q^2<1$); evanescent components are inaccessible in principle, and the *range* of homogeneous components sampled is further restricted in practice to whatever illumination/detection angles the NA admits.

## Overlap with the shear blind spot — settled

[[finite-shear-null-experiment]] found that no set of *in-plane* shear vectors recovers the $K_x = K_y = 0$ axis — a shear along $\mathbf{d}$ is blind wherever $K\!\cdot\!\mathbf{d} = 0$ — and asked whether that axis is already inside the missing cone.

**Answered 2026-09-21 against the real 3D sampling ([[multishear-fusion-design-solution]]): yes, and the question was pointed at the wrong set.** The $K_z$ axis carries exactly zero cap samples except at $\mathbf{K}=0$, because a cap reaches $\mathbf{K}_\perp=0$ only where $\mathbf{k}_\perp = \mathbf{k}_{i\perp}$. So the blind set shared by *every* in-plane shear set costs nothing — those frequencies were unmeasured anyway, and the two limitations do not compound there.

What does cost is the blind set of a **single** shear, which is not the axis but the whole plane $\mathbf{K}_\perp\!\cdot\!\boldsymbol{\Delta}=0$ — and that plane *is* measured: 0.50 % of the sampled support exactly on it, 4.52 % within a 1 rad/μm slab, against 0.019 % near the axis. A second, non-parallel shear direction buys the plane back and leaves only the free residue. This is a structural argument for direction diversity in 3D that is independent of the error-averaging one the Zhu group measured, and its size — half a percent of the support — is the right order for the modest gains they report.

## The cone follows from the cap geometry

The same study derives the cone rather than asserting it. Writing $c(\mathbf{K}) = -(|\mathbf{K}|^2 + 2k_mK_z\cos\theta)/(2k_mK_\perp\sin\theta)$, the illumination azimuths whose cap passes through $\mathbf{K}$ are the roots of $\cos(\phi-\psi)=c$, so the sampled support is exactly $\{|c|<1\}$. Near the origin that boundary reduces to $K_z = \pm K_\perp\tan\theta$ — half-angle $90^\circ-\theta$, i.e. **30° for GODT's 60° illumination cone** — narrowing to 19° at $K_z = 5.9$ rad/μm as the caps curve. Verified against binned digital caps.

A second, separable budget sits on top of it: **180 projections reach only 76.6 % of the support the geometry allows** on a 25.6 μm field under nearest-neighbour gridding (88.4 % on 12.8 μm). That is angular undersampling, not the missing cone, and it is fixed by more projections or better interpolation ([[pan1983]]) rather than by priors.

## Mitigation approaches

- **Comparative regularizer study ([[lim2015]])**: no single algorithm wins universally. Gerchberg-Papoulis (non-negativity only) is cheapest but weakest (most elongation, most noise-sensitive). Edge-preserving (EP, half-quadratic) regularization suits samples with genuine smooth internal texture (RBC/hepatocyte cytoplasm) but is more NA-dependent. Total variation (TV, split-Bregman) gives the best elongation reduction and least NA-dependence on piecewise-constant samples (microspheres) but introduces cartoon-like/stair-step artifacts on samples with real smooth gradients. [[krauze2016]]'s GTVIC (below) is effectively a response to this exact tradeoff — using TV only to extract shape, not content.
- **GTVIC** ([[krauze2016]]): two-stage regularization — a preliminary TV-minimized (Chambolle-Pock) reconstruction gives a geometrically-clean but structurally-corrupted object silhouette; Otsu-threshold + morphological dilation turns this into a binary 3D mask; the mask is then fed as an external-geometry + nonnegativity constraint into *any* iterative solver (tested on both a diffraction-blind straight-ray solver, GDRA, and a proper Rytov/Fourier solver, FTRA). Cuts iteration count 6-10x and measurably reduces axial elongation/blur, on both synthetic cell-phantom data (Paramecium-like, 45° cone, NA=1.3) and a real fixed C2C12 muscle cell. Diffraction-aware FTRA consistently outperforms diffraction-blind GDRA regardless of regularization — reinforcing Rytov/Fourier-based reconstruction (what GODT also uses) as the better base solver.
- **Unsupervised deep learning in projection space** ([[chung2021-projectiongan]]): GP-initialise, forward-project to 2D at all angles, then learn an optimal-transport map from the blurry (missing-angle) projection distribution to the sharp (measured-angle) one with a CycleGAN, and filtered-backproject. 27 s against TV's 605 s, best RI recovery of the four methods compared (1.424 against a true 1.46 on a silica bead, versus TV 1.388 and GP 1.398), and the Ewald sampling curves visibly disappear from the reconstructed spectrum. Two caveats worth keeping attached to that result: the cone is filled by *generation*, not measurement, which sits awkwardly with a quantitative-RI claim, and the residual 2.5 % bead error is the honest measure of the difference; and the FBP back-end abandons diffraction compensation.
- **Aperture synthesis** ([[kandel2019]]): with a fully open condenser (objective used as its own condenser), frequency coverage becomes the *autocorrelation of the objective pupil function*, capturing in a single measurement the coverage other instruments synthesize through multiple angular measurements. An axis for improving k-space coverage orthogonal to both illumination-angle scanning and shear-direction diversity; whether it composes with GODT's approach is unexamined.
- Not fixable by more illumination angles alone within a single-objective epi/trans geometry — fundamentally a support-of-transfer-function limitation set by NA, wavelength, and illumination-angle range. [[krauze2016]]'s few-view study additionally shows reconstruction quality doesn't improve monotonically with projection count below ~60-120 views — sensitive to exact projection subset due to the iterative solver settling into different local optima.

## Relevance

Per [[godt]]: [[2411-08423-godt]] attributes the blur in its reconstructions (vs. the ideal finite-difference derivative) to this effect (Sec 4.3). The missing cone is a property of the illumination geometry, so shear diversity does not fill it; [[multishear-fusion-design-solution]] §2 and [[colour-multiplexed-shear]] §3 confirm this for shear magnitude, direction and colour.
