---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-07-19
status: current
source: "https://doi.org/10.1364/OE.16.019462"
---

# Partially Coherent Image Formation in the DIC Microscope

**Citation:** Shalin B. Mehta, Colin J. R. Sheppard, *Optics Express* **16**(24), 19462–19479 (2008).
Companion conference abstract by the same authors — *Transfer function analysis of partially coherent phase imaging methods*, OSA FiO 2009, paper FMK6 — applies the same machinery to phase contrast and Hoffman modulation contrast and concludes DIC/DPC have the transfer properties best suited to quantitative imaging. Nothing in the abstract is load-bearing for [[godt]] beyond what the full paper establishes.

**Gist:** The rigorous partially coherent image-formation theory for [[gradient-phase-imaging|DIC]], and the paper that pins down *why* the shear factor takes the form it does. Its central object, $\sin(2\pi m\Delta - \phi)$, is the same finite-shear filter that [[s42005-025-02416-3|GODT]] writes as $(e^{j\Delta K_x}-1)/\Delta$ — so the transfer-null structure was already explicit in 2D DIC theory in 2008 (and in Cogswell & Sheppard 1992). This project's contribution is the tomographic setting.

## Key results

**DIC images the coherent difference of two sheared replicas, but each replica is imaged partially coherently.** This is the paper's headline correction to the prior literature, which had assumed Nomarski-DIC illumination is ordinary Köhler illumination. The effective object spectrum seen by the microscope is

$$T_N(m,n) = 2i\,T(m,n)\sin(2\pi m\Delta - \phi) \qquad \text{(Eq. 8)}$$

with $2\Delta$ the shear and $2\phi$ the bias, and the effective transmission cross-coefficient

$$C_N(m,n;p,q) = 4\,C_{BF}(m,n;p,q)\sin(2\pi m\Delta - \phi)\sin(2\pi p\Delta - \phi) \qquad \text{(Eq. 11)}$$

where $C_{BF}$ is the ordinary brightfield TCC. See [[partially-coherent-transfer-function]] for the general framework.

**Three configurations, distinguished only by how coherence at the shear distance is obtained.** *Nomarski-DIC* uses a matched condenser-side Wollaston prism to make the two illumination beams mutually coherent at the shear separation, so it keeps full condenser aperture. *Köhler-DIC* omits the condenser prism, putting the shear entirely in the objective pupil: $P_K(\xi,\eta) = iP_{BF}\sin(2\pi\xi\Delta - \phi)$ (Eq. 3). *PlasDIC* (Zeiss) is Köhler-DIC plus a condenser slit narrowed in the shear direction.

**The coherence–shear coupling (§5, the paper's own "key conclusion").** With coherence ratio $S = \mathrm{NA}_\text{cond}/\mathrm{NA}_\text{obj}$ and shear $2\Delta = \lambda/4\mathrm{NA}$, all three configurations become equivalent below $S \approx 0.4$ — the condenser prism is redundant there. Above it they diverge sharply: at matched illumination ($S=1$) Köhler-DIC and PlasDIC **cannot image weak phase at all**, because the odd part of $C(m;0)$ vanishes identically. Crucially, that $S \approx 0.4$ threshold *scales with shear*: "For bigger shear, the condenser will need to be stopped down further... systems employing smaller shear can work with a more open condenser."

**Perpendicular to the shear, DIC is just a brightfield microscope.** The TCC slice at $m=p=0$ has even symmetry about $m=-p$ in every configuration (Fig. 3b), and even symmetry transfers amplitude, odd symmetry transfers differential phase. So a single shear direction carries no in-band phase information at all along the perpendicular axis.

**At $S=1$ Nomarski-DIC has the broadest support with no zero crossings** in the odd part of $C(m;p)$, for $2\Delta = 1/4m_0$, $2\phi = \pi/2$. Zero crossings are called out explicitly as harmful: they "could lead to imaging of certain frequencies with reversed contrast and hence artifacts."

## Method

Paraxial Hopkins transmission-cross-coefficient theory, computed as the overlap area of three pupils (condenser intensity, shifted objective pupil, conjugate shifted objective pupil) on a discrete 4D $(m,n;p,q)$ grid in MATLAB, without assuming a weak or slowly varying object — the novelty relative to Cogswell & Sheppard 1992. Reduced to 2D by assuming a 1D specimen. Verified experimentally on a 16 μm mouse intestine section (20×/0.75 NA objective, 0.9 NA condenser, 550 nm filtered halogen) with a motorised condenser sweeping $S = 0.3\ldots1.2$, and on a 50 μm optical fibre as a slowly varying phase object across biases $2\phi = 0, \pi/2, \pi$.

## Relevance

For [[godt]] and [[multi-shear-acquisition]], four things:

1. **The transfer nulls are not a new observation.** $\sin(2\pi m\Delta - \phi)$ has the same zero comb that [[finite-shear-null-experiment]] found in GODT's $(e^{j\Delta K_x}-1)/\Delta$, and this paper (via Cogswell & Sheppard 1992, its Eq. 13) predates it by decades, naming contrast reversal at zero crossings as an artifact source. [[finite-shear-null-experiment]] takes this known 2D DIC result into the 3D tomographic setting and quantifies it there.

2. **Bias is a second knob that moves the nulls.** DIC's zeros sit at $2\pi m\Delta - \phi = k\pi$, so the bias $\phi$ translates the whole comb. GODT has no bias term because its phase-shifting recovers the complex gradient directly, but this raises the question of whether bias diversity could substitute for, or cheapen, magnitude diversity. *(inference — no source combines the two.)*

3. **A possible constraint on multi-shear.** Larger shear demands a more closed condenser to retain coherence across the shear distance, so a scheme mixing magnitudes might need different illumination apertures per channel. This bites only for Köhler-like geometries; Nomarski-DIC escapes it via the matched condenser prism.

4. **Quantitative backing for the directional-contrast motivation.** "Single shear gives orientation-dependent contrast" is usually stated qualitatively; here the perpendicular-direction TCC is shown to be *brightfield*, i.e. zero in-band phase transfer, not merely reduced.

## Raises

- Cogswell & Sheppard, *J. Microsc.* **165**, 81–101 (1992), Eq. 13 — the origin of the DIC partially coherent transfer function, and the primary citation for the null structure.
- Arnison, Larkin, Sheppard, Smith & Cogswell, *J. Microsc.* **214**, 7–12 (2004), "Linear phase imaging using DIC" — the two-orthogonal-shear 2D reconstruction that the multi-shear extension generalises.
