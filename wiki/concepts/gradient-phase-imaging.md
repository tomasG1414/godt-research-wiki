---
type: concept
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Gradient-Based Phase Imaging (family)

Techniques that measure a spatial *derivative* of phase rather than phase itself, tracing back to Nomarski differential interference contrast (DIC) microscopy. Gradient contrast suppresses slowly-varying background, giving high SNR/contrast for subtle features. Differentiating a phase map *numerically* amplifies high-frequency noise; an optical shear measures the phase difference directly, so it avoids that step ([[2411-08423-godt|GODT]] draws this contrast with numerical gradient operations).

## Members (as surveyed by [[2411-08423-godt]])

- Programmable 2D/3D quantitative differential phase contrast (qDPC)
- Gradient Light Interference Microscopy (GLIM) — [[nguyen2017]], epi-illumination variant [[kandel2019]]. Closest methodological predecessor to GODT: same shear-direction phase-gradient measurement (Wollaston prism + SLM phase shifting instead of polarization gratings + rotating polarizer), same low-coherence-for-SNR philosophy, same thick-sample motivation. Suppresses multiple-scattering background because both sheared copies accumulate the *same* incoherent background, which cancels in the difference. Reaches 3D via gradient integration + z-scanning/optical sectioning, **not** tomographic inversion — precisely the gap GODT closes.
- Gradient Retardance Optical Microscopy (GROM)
- Quadriwave Lateral Shearing Interferometric Microscopy (QWLI)
- Metasurface-enabled single-shot phase-gradient microscopy — [[kwon2019]]: two cascaded dielectric metasurface layers (~1 mm³ total) capture three phase-offset DIC images simultaneously via polarization + spatial multiplexing, giving single-shot gradient recovery with no moving parts. Still single-shear, still 2D, still proposes z-scanning/TIE for depth
- Oblique Plane Back-Illuminated Microscopy / quantitative oblique back-illumination microscopy (qOBM) — [[optica-8-1-6]], extended to opaque-sample 3D RI tomography via a non-paraxial 3D OTF; gradient signal comes from subtracting diametrically-opposed LED illumination pairs rather than beam shearing, and shares the missing-cone limitation of any limited-angle tomography

- Double Sagnac lateral shearing QPM — [[2605-04297-double-sagnac-ls-qpm]]: two orthogonal shear directions captured in a **single shot** as a two-carrier mesh interferogram, under [[pseudo-thermal-illumination|pseudo-thermal]] illumination. The family's most recent member and the one that most directly demonstrates two-direction acquisition; still integrates to 2D phase and still reaches 3D by z-scanning (§3.6, ~4.5 µm axial resolution), so it does not close the gap below.

## The gap GODT fills

Until [[2411-08423-godt|GODT]], techniques in this family produced **2D integrated** quasi-quantitative phase images: 3D, where attempted, came from z-scanning (refocusing at different depths) rather than tomographic reconstruction, which loses correlations in the sample's 3D structure and isn't diffraction-compensated. Getting to a 3D RI map also traditionally required numerically **integrating** the measured gradient back into a phase map first, a step prone to experimental and numerical error. A contemporaneous QWLI-based ODT concept (Yuan et al., ref 61 in the GODT paper) does attempt tomography from gradient data but integrates the gradient first, then runs standard ODT.

GODT's distinguishing move: skip integration entirely and feed the gradient (shearing-interferometry) data directly into the [[rytov-approximation|Rytov]]-based tomographic reconstruction, producing a genuinely tomographic (diffraction-compensated, full 3D) reconstruction of the RI-derivative field.

Reading [[nguyen2017]] alongside [[2411-08423-godt]] clarifies GODT's contribution: not the gradient *measurement* (GLIM and this whole family already did that, with essentially the same sheared-copy interferometry), but feeding gradient data straight into tomographic *reconstruction*.

## The underlying theory (DIC, partially coherent)

[[mehta-sheppard2008]] supplies the rigorous image-formation theory this whole family rests on, and corrects a standing misconception: Nomarski-DIC illumination is *not* ordinary Köhler illumination, so DIC images the coherent difference of two sheared replicas while each replica is imaged **partially coherently**. The effective object spectrum is $T_N(m,n) = 2iT(m,n)\sin(2\pi m\Delta - \phi)$ — the finite-shear filter, with $2\Delta$ the shear and $2\phi$ the bias — and the resulting transmission cross-coefficient is the brightfield one multiplied by that factor twice over. See [[partially-coherent-transfer-function]].

Three consequences the family shares:

- **Symmetry decides what is transferred.** Even symmetry of the TCC about $m=-p$ carries amplitude, odd symmetry carries differential phase. Perpendicular to the shear the TCC is even — so a single shear direction is, along that axis, a plain brightfield microscope with **zero** in-band phase transfer.
- **Shear magnitude is coupled to illumination coherence.** All three DIC configurations (Nomarski, Köhler, PlasDIC) coincide below $S \approx 0.4$, and that threshold scales with shear: bigger shear demands a more closed condenser. At $S=1$, Köhler-DIC and PlasDIC cannot image weak phase at all.
- **Zero crossings are a known artifact source**, causing reversed contrast at affected frequencies — the 2D ancestor of the transfer nulls quantified in [[finite-shear-null-experiment]].

## Shear-induced transfer-function asymmetry

Independently confirmed as a real obstacle by [[kandel2019]], which states outright that "the shear in GLIM renders the transfer function non-radially symmetric" and that this blocks straightforward 3D reconstruction. This is the same anisotropy [[2411-08423-godt]] reports as its single-shear contrast limitation. The two lineages respond differently:

- **GLIM/epi-GLIM** sidesteps analytic inversion entirely — measures the PSF empirically using co-localized fiducial microspheres, then Wiener-deconvolves.
- **GODT** derives the exact $\tilde{O}_\Delta = iK_x\tilde{O}$ relation and accepts the resulting directional contrast.
- **A multi-shear extension** would be a third response: combine shear directions to restore symmetry. The GODT paper names this as future work.

[[2605-04297-double-sagnac-ls-qpm]] does take it, but only in 2D and for a different stated reason — removing the line artifacts that 1D integration of a single-direction gradient produces, rather than restoring transfer-function symmetry. Since GODT never integrates, that particular motivation does not transfer; the hardware demonstration does.

## Multi-shear precedent in 2D

DIC-family work (Cogswell, Shribak, and the two-orthogonal-shear DIC lineage cited in [[nguyen2017]]) already established that two orthogonal shear directions recover the full transverse phase gradient in **2D**. A multi-shear extension of GODT (named as future work in the paper) would lift this into the 3D tomographic (Ewald-space) setting.

## Relevance

Per [[godt]]: this family is GODT's direct lineage and comparison set; [[nguyen2017]], [[kandel2019]] and [[kwon2019]] cover specific members.
