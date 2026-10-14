---
type: project
projects: [godt]
created: 2026-07-18
updated: 2026-10-14
status: current
---


# GODT — Gradient Optical Diffraction Tomography

## Goal

Understand gradient ODT ([[s42005-025-02416-3]], *Communications Physics* 8, 510 (2025);
preprint [[2411-08423-godt]]) in depth, then study one extension numerically: acquiring
and fusing **several shear channels** (different shear magnitudes and directions) in a
single tomographic reconstruction. Everything here runs on synthetic data, with no lab needed.

## Background

**Core method.** GODT feeds common-path [[common-path-shearing-interferometry|shearing-interferometry]] phase-gradient data directly into a first-order [[rytov-approximation|Rytov]] tomographic reconstruction ([[optical-diffraction-tomography]]'s Ewald-sphere / direct-inversion pipeline), instead of first numerically integrating the gradient into a phase map. The central result is $\tilde{O}_\Delta = iK_x\tilde{O}$: shear-derivative data maps onto the exact derivative of the ODT scattering-potential spectrum, so the diffraction-compensating machinery carries over unchanged. The published version generalises this to finite shear, with exact factor $\frac{e^{j\Delta K_x}-1}{\Delta}$ and $iK_x$ as its $\Delta\to0$ limit. Validated numerically ([[shepp-logan-phantom]]) and experimentally ([[nano-printed-cell-phantom]], hippocampal neurons) on a polarization-grating common-path system.

**Why it matters.** Total-shear common-path ODT needs a sample-free region for self-referencing, so it fails on dense, confluent samples. GODT's small-shear regime needs no such region and works with low-coherence LED illumination, so it images confluent samples (neurons) that total-shear mode cannot.

**Lineage.** The gradient *measurement* is shared with the whole [[gradient-phase-imaging]] family ([[nguyen2017]], [[kandel2019]], [[kwon2019]], [[optica-8-1-6]]). Reconstructing a 3D index derivative straight from gradient data also has precedent with straight-ray tomography ([[vishnyakov2016]]) and in X-ray phase tomography ([[guizar-sicairos2011]]). GODT does this inside a rigorous *diffraction* framework, which is the distinction its published version draws itself.

## Why several shear channels

See [[multi-shear-acquisition]] for the full picture. The two motivations are not equally strong:

1. **Directional contrast.** A single shear sees index variation along one direction only ([[kandel2019]]; [[mehta-sheppard2008]] shows the transfer perpendicular to the shear is brightfield). Measured gains from extra *directions* saturate quickly, at about one orthogonal pair ([[zhu2024-appl-phys-b]], [[zhu2024-sensors]], [[zhu2025-sci-rep]]).
2. **Filling transfer nulls by magnitude diversity.** The finite-shear factor vanishes at $K=2\pi m/\Delta$. The null comb itself is known in 2D ([[servin2007]], [[mehta-sheppard2008]]). This project quantifies it in the 3D Ewald setting and shows that a second, incommensurate shear magnitude fills it ([[finite-shear-null-experiment]]). That recovers information that is otherwise absent at any SNR, which is stronger than noise averaging.

Acquisition hardware for two shear directions already exists in several forms ([[nguyen2025-deep-learning-lsi]], [[sciadv-adl0501]], [[2605-04297-double-sagnac-ls-qpm]]), and colour can act as a shear magnitude ([[colour-multiplexed-shear]]). The question studied here is the reconstruction side: how to fuse shear channels in $\tilde O$-space, and how to choose them.

## Results of the simulation thread

Code lives in `lab/` (see `lab/README.md`). Every number below is reproduced by a script.

**Step 1 — transfer nulls** ([[finite-shear-null-experiment]]). At 0.8 μm shear, three nulls fall inside the measurable band and cost 3.56 % of it along the shear axis. At 0.1 μm there are none. Two design rules follow: avoid commensurate magnitude ratios, and use magnitudes and directions for different defects.

**Step 2 — fusion and optimal design in 3D** ([[multishear-fusion-design-problem]], [[multishear-fusion-design-solution]], derivation in [[multishear-fusion-derivation]]):

1. **The problem separates exactly**: $J(\mathbf K)=\rho\,S_{\mathbf t}(\mathbf K_\perp)\,M(\mathbf K)$, with a closed form for the cap density $M$. Shear design and tomographic geometry never mix, so the whole study is minutes of compute.
2. **The real blind spot is a plane, not the axis.** The $K_z$ axis is unsampled by transmission ODT anyway. A single shear's blind *plane* costs 0.50 % of the support, and a second direction buys it back.
3. **Lost volume at 0.8 μm is 3.02 %**, split 83.6 % null planes / 16.4 % blind plane. Single-shear quality is a trade between a null comb and a blind slab, not monotone in $\Delta$.
4. **Design rule: $\Delta\approx2\pi/K_{\max}$** (0.20 μm for the simulation optics), as large as possible without admitting the null comb.
5. **Per-frequency weighting is worth 21×** over equal-weight averaging of reconstructions at equal photons in the linear model. Fusion equals precision-weighted averaging, so this is a standard estimator applied well, not a new one.
6. **Calibration can limit the high frequencies**: assuming 10 nm rms random shear error, the band-edge SNR ceiling is 20 dB below the photon limit. A broadband LED already partly fills the higher-order nulls.

**Step 2d — colour as a shear channel** ([[colour-multiplexed-shear]]). Three LED colours act like a mechanical pair of shear magnitudes, blue + green is the best pair, and at the design-rule shear the trick is unnecessary. The separation theorem breaks across colours, and the missing cone is untouched.

**Step 3 — multiple scattering** ([[multiple-scattering-check]]), with the same estimator on data from a [[brenner1993]] WPM forward model:

1. **The nulls survive multiple scattering exactly.** Shearing acts on the measured field after all propagation, so the geometry of step 2 is unconditional (checked to 4e-16).
2. **The margins shrink with contrast**: 9.6× under the ideal linear model, 1.8× with realistic gridding, 1.1× with multiple scattering at [[chen1998]]'s validity bound.
3. **Gridding matters.** Moving from binned to continuum $K_z$ data with the real pupil costs a factor 6.2 under nearest-neighbour deposition, the largest correctable error in the pipeline ([[pan1983]]); a bilinear ablation will separate gridding from pupil effects.
4. **WPM vs. first Born (provisional).** Multi-slice WPM takes the obliquity factor at the incident ray; on thick volumes this implementation shows a larger residual that may partly be a simulation-window effect and needs a wider-window re-run.
5. **Shear data need no unwrapping step.** Against ODT *without* unwrapping, at high contrast, this is worth 24× in a strictly single-scattering pipeline; the comparison against unwrapped ODT is still to run, and since the shear difference has its own wrap limit (2.2 % of pixels at 0.8 μm vs 0.45 % one-pixel residues) the net benefit there may be small.

## Next technical steps

- **Bilinear (or NUFFT) gridding.** The cheapest fix with the largest measured effect (the factor 6.2 above), and an ablation that separates gridding from pupil effects.
- **Fairer baselines and checks.** An unwrapped-ODT baseline for the wrapping comparison, and a wider-window WPM run for the forward-model gap.
- **Spatially varying shear.** Real two-direction instruments show shear varying across the field ([[2605-04297-double-sagnac-ls-qpm]]). A position-dependent $\Delta$ moves the nulls, which could average them away or blur any inverse filter.
- **Born-series / Lippmann–Schwinger forward model**, whose linear limit is first Born by construction, to remove the splice that [[multiple-scattering-check]] needed.
- **Partial coherence.** Whether the shear factor stays a simple multiplier under [[ledwig-robles2023]]'s broadband Ewald-plane treatment.
- **Real samples.** Placing GODT's targets (neurons, [[nano-printed-cell-phantom]]) on the $a\,\delta n/\lambda$ axis, which decides how much of the design study applies to them.
