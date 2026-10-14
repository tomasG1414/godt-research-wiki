---
type: concept
projects: [godt]
created: 2026-07-19
updated: 2026-07-19
status: current
---

# Partially Coherent Transfer Functions

The framework describing how a microscope transfers object spatial frequencies to image intensity when the illumination is neither fully coherent nor fully incoherent — the regime every LED- or lamp-illuminated instrument in this corpus actually operates in, including [[s42005-025-02416-3|GODT]].

## The core object

Under Hopkins' theory the image intensity is *bilinear* in the object spectrum, not linear:

$$I(x,y) = \iiiint T(m,n)\,T^*(p,q)\,C(m,n;p,q)\,e^{2\pi i[(m-p)x+(n-q)y]}\,dm\,dn\,dp\,dq$$

The **transmission cross-coefficient** (TCC) $C(m,n;p,q)$ is the overlap area of three pupils — the condenser intensity, the objective pupil shifted by $(m,n)$, and its conjugate shifted by $(p,q)$ ([[mehta-sheppard2008]], Eqs. 2–4). The consequence that matters: the strength with which a frequency appears in the image depends on *pairs* of object frequencies, so there is no single "OTF" unless the object is restricted.

Two 1D reductions recover a linear picture, both named in [[mehta-sheppard2008]]:

- **Weak object transfer function (WOTF)**, $C(m;0)$ — valid for weak objects. Its *even* part transfers amplitude, its *odd* part transfers phase gradient.
- **Phase gradient transfer function (PGTF)**, $C(m;m)$ — valid for slowly varying objects. Its shape predicts how image intensity tracks the local phase slope, and matches measured intensity across a 50 μm fibre.

## The symmetry rule

Symmetry of $C(m;p)$ about the $m=-p$ axis determines what kind of information reaches the image: **even symmetry transfers amplitude, odd symmetry transfers differential phase**. This is the cleanest available statement of why shearing produces phase contrast at all, and it immediately gives the directional-contrast result — perpendicular to the shear the TCC is even, i.e. brightfield, so a single shear direction transfers no in-band phase information along that axis at all. See [[gradient-phase-imaging]] and [[multi-shear-acquisition]].

## Coherence ratio

$S = \mathrm{NA}_\text{cond}/\mathrm{NA}_\text{obj}$ parameterises the whole family. Lowering $S$ raises the odd (phase) part of the TCC relative to the even (amplitude) part, so **smaller $S$ gives better phase-gradient contrast** — while costing lateral resolution and depth sectioning. For DIC configurations the practical consequence is a coupling between $S$ and shear magnitude, discussed in [[mehta-sheppard2008]]: bigger shear forces a more closed condenser. This is the main reason to care about the framework in this project.

## The broadband extension

[[ledwig-robles2023]] generalises the picture to arbitrary sources by folding temporal frequency into the radial direction of a 3D reciprocal space, so angle and wavelength become one distribution $S(\vec u)$. The sampling surface then becomes an **Ewald plane** rather than the Ewald sphere cap of coherent [[optical-diffraction-tomography]] — the set of source-field pairs contributing incoherently to a given image frequency, independent of wavelength. Amplitude signal is the sum of the two sampling planes, phase the difference. This removes the symmetry and discreteness assumptions that earlier partially coherent 3D treatments needed, and makes computation $O(n^3)$.

## Why this matters for GODT

GODT measures a phase gradient under broadband LED illumination and reconstructs with a coherent monochromatic [[rytov-approximation|Rytov]] model. Everything on this page describes the gap between those two statements. Three specific consequences are tracked elsewhere:

1. The finite-shear factor $\sin(2\pi m\Delta - \phi)$ appearing inside the DIC TCC is the same filter as GODT's $(e^{j\Delta K_x}-1)/\Delta$ — so the transfer nulls of [[finite-shear-null-experiment]] have 2D DIC precedent going back to Cogswell & Sheppard (1992).
2. The bias $\phi$ translates the null comb, which GODT's phase-shifting scheme discards. Whether bias diversity is a cheaper substitute for shear-magnitude diversity is open. *(inference)*
3. Whether GODT's actual broadband transfer function differs materially from its assumed coherent one is computable with [[ledwig-robles2023]]'s method and currently unknown.
