---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-09-22
status: current
source: "https://doi.org/10.1364/AO.32.004984"
---

# Light Propagation Through Microlenses: A New Simulation Method

**Citation:** K.-H. Brenner and W. Singer, "Light propagation through microlenses: a new simulation method," Applied Optics 32(26), 4984–4988 (1993). Angewandte Optik, Universität Erlangen-Nürnberg.

**Gist:** Introduces the **wave-propagation method (WPM)**, a forward-simulation technique for light propagating through inhomogeneous-refractive-index optical elements (microlenses) that removes the paraxial/small-angle restriction of the standard beam-propagation method (BPM), at the cost of higher computational effort. [[2411-08423-godt]] cites this as the method used to generate its synthetic scattering dataset (180 object waves through the Shepp-Logan phantom) — i.e., this is GODT's forward-model engine, not a tomography reconstruction method itself.

## Key results

- **BPM's core limitation**: BPM solves the scalar Helmholtz equation by splitting each propagation step into (1) homogeneous-medium propagation in Fourier space using an *averaged* refractive index, and (2) an index-inhomogeneity phase correction applied in the spatial domain with no angular sensitivity. Both steps introduce errors proportional to the propagation angle — BPM is numerically reliable only for propagation angles roughly ±10° (NA ≲ 0.17 in air); beyond that it manifests as artificial spherical aberration.
- **WPM's fix**: instead of averaging the index before propagating, WPM keeps the full plane-wave angular-spectrum decomposition and applies the *local, angle-dependent* phase term $\exp[-in(x,z)k_0\cos(\theta)\delta z]$ to each spatial-frequency component separately, with the propagation angle $\theta$ tied to each component via $n(x,z)\sin(\theta) = \lambda_0 v_x$ (a per-frequency-component Snell's law). This removes the paraxial restriction entirely — validated against exact Snell's-law refraction up to propagation angles beyond 70° (only Fresnel-amplitude effects, which neither method models, become significant beyond that).
- Demonstrated on Maxwell's fisheye lens (ideal point-to-point focusing gradient-index profile): BPM shows visible spherical-aberration-like defocus at high-NA parts of the field; WPM reproduces the theoretically correct focal position with near-symmetric intensity distribution.
- **Cost**: WPM is ~100× slower than BPM for the same grid (1024 samples, δz=0.1 μm test case took minutes with BPM, ~100× longer with WPM) — a genuine speed/angle-range tradeoff, not a strict improvement.
- Method limitations shared with BPM: only forward propagation (no back-reflection), no Fresnel amplitude transmission/reflection at interfaces (differential rather than integral Maxwell's-equations formulation), evanescent-field contributions dropped (plane-wave expansion only propagates homogeneous components) — the same fundamental restriction to propagating waves that underlies the [[missing-cone-problem]] in [[wolf1969]]'s original theorem.

## Method

2D scalar-field simulations (extendable to 3D), angular-spectrum decomposition via Fourier transform, per-component phase propagation with local (position- and angle-dependent) refractive index rather than a globally-averaged one. Tested on tilted-interface refraction (vs. Snell's law), Maxwell's fisheye, hyperboloidal lenses, and multi-plane-wave superposition inputs.

## Relevance

Per [[godt]]: this is the forward-simulation engine behind [[2411-08423-godt]]'s numerical Shepp-Logan validation (its ref 78) — not a reconstruction algorithm, but the tool that generates the 180 synthetic object waves later fed into [[optical-diffraction-tomography|DI]] reconstruction. For the multi-shear extension's simulation thread, WPM (or an equivalent full-angular-spectrum forward solver) is the natural choice for generating synthetic multi-shear data, since GODT's conical illumination scheme (60° inclination in the core paper) is well outside BPM's ~10° validity range — BPM would not be an acceptable substitute here.

## Later use as an inverse-problem forward model

[[ma2017]] carries WPM from pure simulation into tomographic *reconstruction*, using it as a nonlinear forward operator inside a FISTA optimization loop (with an efficient scheme for the field derivative) to handle multiply-scattering objects that break the Born/Rytov linearization. Relevant if a multi-shear GODT extension takes an optimization-based form and wants a forward model valid beyond single scattering.

## Measured here: WPM vs. first Born (provisional)

Implemented and exercised in `lab/shearlib.py` for [[multiple-scattering-check]].
Two things that matter and are not in the paper.

**The obliquity factor is taken at the incident ray.** For a single slice WPM
reduces analytically to the transmission
$t(\mathbf r)=\exp[i\,k_z(n(\mathbf r),\mathbf k_{i\perp})\,\delta z]$ — the local
axial wavenumber evaluated at the *incident* transverse frequency — so its
scattered amplitude carries $1/k_{iz}$ where the first Born approximation carries
$1/k_z(\mathbf k_s)$, the *scattered* one. The ratio is exactly $k_z/k_{iz}$,
verified numerically to 1e-3 (and the transmission itself to 1e-12) in
`shearlib.selftest`. The two agree for forward scattering and part company as the
scattering angle grows. This is a property of multi-slice propagation, shared with
BPM, and is separate from the paraxial limitation the paper *does* discuss.

**Consequence for this project.** On a thick cell-like object at 60° illumination
this implementation shows a WPM-vs-Born gap of ~0.67 relative that converges to a floor
as $\delta z\to0$. Part of it may be a window/absorber effect (see
[[multiple-scattering-check]] §3), so it is provisional. It does not affect GODT's own
Shepp–Logan check, which compares GODT with differentiated ODT on the same WPM data. So validating a [[rytov-approximation|Rytov]]/Ewald pipeline on
WPM-generated data tests the pipeline and this forward-model gap together. Anyone
validating a multi-shear design this way needs a forward model whose linear limit
*is* the inverse's model — which is the argument for [[ma2017]]'s route, or for a
Born-series solver.

## Raises

- WPM's ~100× compute cost vs. BPM is a real consideration if a multi-shear extension needs to generate large synthetic datasets (multiple shear directions × many illumination angles) — worth checking whether [[2411-08423-godt]]'s actual simulation runtime/grid size (1000×1000, dx=0.1μm) is documented anywhere, to estimate compute budget for a multi-shear synthetic study.
