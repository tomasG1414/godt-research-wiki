---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-07-22
status: current
source: "https://doi.org/10.1364/OE.19.021345"
---

# Phase Tomography from X-Ray Coherent Diffractive Imaging Projections

**Citation:** M. Guizar-Sicairos, A. Diaz, M. Holler, M. S. Lucas, A. Menzel, R. A. Wepf, O. Bunk, "Phase tomography from x-ray coherent diffractive imaging projections," Optics Express 19(22), 21345–21357 (2011). Paul Scherrer Institut / ETH Zurich.

**Gist:** Post-processing and reconstruction algorithms for X-ray ptychographic phase tomography: removing constant/linear phase offsets, fiducial-free automatic projection alignment, and — most relevant here — **tomographic reconstruction directly from the phase derivative, computed in a way that is inherently insensitive to phase wrapping**, eliminating the phase-unwrapping step entirely. Demonstrated on bacteria populations in legume root-nodule cells.

## Key results

- **Direct reconstruction from the phase derivative** (Eq. 12): a modified Fourier filter — previously introduced for beam-deflection tomography and X-ray differential-phase-contrast tomography — allows the tomogram to be computed straight from $\partial\varphi/\partial x$:
$$-\frac{2\pi\Delta}{\lambda}\delta(\mathbf{r}) = \int_0^\pi \mathcal{F}^{-1}_u\left\{\frac{1}{i2\pi\,\mathrm{sgn}(u)}\mathcal{F}_x\left\{\frac{\partial\varphi(x,y,\theta)}{\partial x}\right\}\right\}d\theta$$
"removing the need for an intermediate integration step that may introduce noise or artifacts." The derivation follows from the same FT differentiation identity GODT uses, $\mathcal{F}_x\{\partial\varphi/\partial x\} = i2\pi u\,\mathcal{F}_x\{\varphi\}$ — here inverted as a $1/(i2\pi\,\mathrm{sgn}(u))$ filter within filtered backprojection, rather than mapped onto Ewald caps.
- **Wrapping-insensitive derivative estimation** (Eqs. 13–14) — the key transferable idea. Instead of differentiating the wrapped phase directly (which produces artifacts at wrap discontinuities), define $g(x,y,\theta) = \exp[i\varphi(x,y,\theta)]$ and estimate
$$\frac{\partial\varphi}{\partial x} \approx \frac{1}{h}\arg\left[g(x+h/2,y,\theta)\,g^*(x-h/2,y,\theta)\right]$$
Because the phase is mapped to the unit circle, wrapping discontinuities simply disappear and $g$ is well-behaved. **Caveat retained by the authors**: this still requires the phase to change by less than π between adjacent pixels, which is the sampling requirement for any complex-valued function.
- Fiducial-free automatic alignment of tomographic projections, reported as robust and applicable beyond X-ray to lens-based nanoscale tomography; validated against independent interferometer measurements (RMS agreement 40–55 nm). Phase unwrapping is needed *only* for this alignment step, not for the reconstruction itself.
- Simple method for removing undesired constant and linear phase terms from ptychographic reconstructions (which recover phase only up to such an offset).

## Relevance

Per [[godt]]: **this substantially answers the open question filed against [[chen1998]]** about whether GODT's gradient formulation escapes the phase-unwrapping bottleneck that caps standard Rytov ODT's validity range.

[[chen1998]] showed first-Rytov ODT's ~3× validity advantage over first-Born collapses abruptly once its criterion is violated, and identified phase unwrapping as the mechanism. This paper demonstrates that a phase *derivative* estimated as $\frac{1}{h}\arg[g(x+h/2)g^*(x-h/2)]$ is structurally immune to wrapping. That expression is essentially what a shearing interferometer measures physically: GODT's fringe phase is $\arg$ of the product of the object field with its sheared conjugate replica, i.e. the same complex-product construction, realized in optics rather than numerics. So GODT plausibly inherits the unwrapping immunity *by construction* — the wrapping problem never arises because absolute phase is never formed.

Two caveats keep this from being a clean win. First, the π-per-pixel sampling requirement survives — GODT can only tolerate phase gradients up to that bound, so the constraint is relaxed rather than removed. Second, no ingested source states this connection explicitly; it is *(inference)* from combining this paper's Eq. 13 with GODT's measurement model, and should be verified before being relied on.

Independently useful: the direct-from-derivative filter (Eq. 12) is a *non-diffractive* (filtered-backprojection) analogue of GODT's Ewald-space relation, with X-ray precedent predating both GODT and [[vishnyakov2016]]. This reinforces that "reconstruct straight from the gradient without integrating" is a recurring idea across imaging modalities, and that GODT's specific contribution is doing it inside a rigorous *diffraction* framework — the same conclusion reached from the [[vishnyakov2016]] comparison, now from a third direction.

The fiducial-free alignment algorithm would also matter for experimental multi-shear data, where per-shear-channel registration errors would matter — note [[2509-15382-osi-flex]]'s finding that 0.1-pixel shear error already degrades reconstruction.

## Raises

- The π-per-pixel sampling bound sets a hard limit on how large a phase gradient GODT can measure regardless of unwrapping immunity. Whether GODT's samples approach it is unquantified by any ingested source.
