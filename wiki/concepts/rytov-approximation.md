---
type: concept
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Rytov Approximation (first-order)

Weak-scattering approximation used to linearize the relationship between a measured scattered field and a sample's 3D scattering potential, underpinning [[optical-diffraction-tomography]]'s Wolf projection theorem. Defines the Rytov-transformed field:

$$u_R = \ln(u/u_i) = \ln(A/A_i) + i(\varphi - \varphi_i)$$

(log of complex amplitude ratio, rather than the raw field difference used in the first Born approximation). Its Fourier transform maps directly onto the Ewald-sphere cap of the object spectrum $\tilde{O}$.

## Rytov vs. Born

The first-order Born approximation instead linearizes the raw scattered field difference $u - u_i$. Rytov is known to hold under **less strict weak-scattering conditions** than Born (larger phase excursions tolerated) — this is why ODT/GODT papers in this project default to Rytov.

**Quantitative validity ([[chen1998]], via cylinder simulations with exactly-known scattering)**: validity of each approximation scales with the product of object size $a$ and refractive-index difference $\delta n$, not $\delta n$ alone. Criteria for 25%-error reconstructions: first Born (plain filtered backpropagation) $a\,\delta n \le 0.08\lambda$; first Rytov (*hybrid* FBP — requires exact field backpropagation to the image center before Rytov-transforming) $a\,\delta n \le 0.24\lambda$, a **3× larger** validity range. This overturned an earlier claim (Slaney et al. 1984) that Rytov validity is size-independent. Crucially, Rytov's failure mode past its validity bound is abrupt collapse via **phase-unwrapping breakdown**, not graceful degradation like Born — extending Rytov's practical range further is a phase-unwrapping problem, not a weak-scattering-theory problem.

### What shear data avoid: the unwrapping step, not the linearisation

[[guizar-sicairos2011]] shows a phase derivative computed as $\frac{1}{h}\arg[g(x+h/2)g^*(x-h/2)]$, with $g=\exp(i\varphi)$, is inherently wrapping-insensitive — the phase lives on the unit circle, so wrap discontinuities never appear. A shearing interferometer measures exactly this construction optically (fringe phase = argument of the field times its sheared conjugate replica), so GODT never forms the absolute phase that would need unwrapping.

**Confirmed numerically** in [[multiple-scattering-check]], on multiply-scattered fields rather than by argument. Both routes to the same shear difference were computed and compared: $\log[u(\mathbf r+\boldsymbol\Delta)/u(\mathbf r)]$, which the interferogram delivers, against $\log u(\mathbf r+\boldsymbol\Delta)-\log u(\mathbf r)$, which is what you hold if the absolute phase is formed first. They disagree exactly where a $2\pi$ jump separates the two sampled points, and that fraction is what the shearing geometry gets for free.

**But the surviving constraint is more binding than the page previously implied.** It is $|\Delta\varphi|<\pi$ across the **shear distance**, not across a pixel, because $\arg[\cdot]$ is itself a principal value — and at GODT's published 0.8 μm shear that is the *tighter* of the two conditions. Measured: at $a\,\delta n=0.92\lambda$ the absolute phase has residues at 0.45 % of pixel pairs, while the 0.8 μm shear difference exceeds $\pi$ at 2.2 % and the 1.6 μm one at 4.2 %, roughly doubling with $\Delta$. Below $a\,\delta n\approx0.46\lambda$ neither wraps.

So GODT does not need an unwrapping step at all: violations stay on the offending pixels instead of propagating through an unwrap. It still relies on the Rytov linearisation, and the shear difference has its own wrap limit; a comparison with properly unwrapped ODT has not been run. It is also a third independent reason to keep $\Delta$ small, alongside $\Delta\approx2\pi/K_{\max}$ from [[multishear-fusion-design-solution]]. See [[multiple-scattering-check]] §6.

## Role in GODT

[[2411-08423-godt]]'s core derivation differentiates the Rytov field in the shear direction ($u_{\Delta R} = \tfrac{1}{A}\tfrac{dA}{dx} + i\tfrac{d\varphi}{dx}$) and shows, via the FT differentiation property ($\tfrac{d}{dx}f(x) \leftrightarrow iK_x\tilde f(K_x)$), that the shear-derivative Rytov field maps onto $iK_x\tilde{O}$ rather than $\tilde{O}$ — i.e., gradient interferometry data is compatible with the *same* Rytov/Ewald-sphere machinery as full-field ODT, just with a known linear filter ($iK_x$) applied. This is the mathematical foundation for GODT.

## Where it actually breaks, for this project's phantom

[[multiple-scattering-check]] measures the linearisation error directly, by running a [[brenner1993]] WPM forward model against its own linearisation — so forward-model differences cancel and only the nonlinearity is left — on a 6 μm cell-like object under GODT's 60° conical illumination:

| $a\,\delta n/\lambda$ | linearisation error | at $\lvert K\rvert>20$ rad/μm |
|---|---|---|
| 0.11 | 0.070 | 0.228 |
| 0.23 | 0.149 | 0.481 |
| 0.46 | 0.330 | 1.077 |
| 0.92 | 1.045 | 4.724 |

Two things worth carrying. The error is **about 3× worse at the band edge** than at low frequency (4.7× at the highest contrast), so multiple scattering damages precisely the high transverse frequencies a multi-shear design exists to recover. And at [[chen1998]]'s bound $a\,\delta n=0.24\lambda$ the resulting *reconstruction* error is 0.288, against chen1998's stated ~25 %. The geometry and the metric differ, so this is consistency in magnitude rather than a precise cross-check.

## Relevance

Per [[godt]]: load-bearing approximation for the entire reconstruction pipeline (classical and gradient). [[chen1998]] and [[wolf1969]] provide the theoretical grounding; [[multiple-scattering-check]] prices its failure.
