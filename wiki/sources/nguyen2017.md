---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
source: "https://doi.org/10.1038/s41467-017-00190-7"
---

# Gradient Light Interference Microscopy for 3D Imaging of Unlabeled Specimens

**Citation:** T. H. Nguyen, M. E. Kandel, M. Rubessa, M. B. Wheeler, G. Popescu, "Gradient light interference microscopy for 3D imaging of unlabeled specimens," Nature Communications 8, 210 (2017). University of Illinois at Urbana-Champaign (Popescu group). Cited as ref. 49 in [[2411-08423-godt]].

**Gist:** Introduces GLIM (Gradient Light Interference Microscopy) — a DIC-based, low-coherence, phase-gradient QPI method that suppresses the multiple-scattering incoherent background in thick specimens by interfering two transversely-sheared copies of the field that suffer *equal* scattering degradation. A foundational member of the [[gradient-phase-imaging]] family and one of GODT's most direct conceptual predecessors: it measures the same shear-direction phase gradient GODT does, but recovers 3D structure via z-scanning + optical sectioning rather than true diffraction tomography.

## Key results

- **Measurement principle nearly identical to GODT's**: a Wollaston prism (in a DIC microscope) generates two cross-polarized field replicas shifted by less than a diffraction spot ($\delta r$); an SLM applies controlled phase shifts $\phi_n = n\pi/2$ ($n=0,1,2,3$) between them; from the 4 intensity images the phase gradient along the shear direction is recovered, $\nabla_x\phi(r) \approx \Delta\phi(r)/\delta x$. This is the same finite-difference-of-sheared-copies gradient measurement as [[common-path-shearing-interferometry|GODT's shearing interferometry]], just implemented with a Wollaston prism + SLM rather than polarization gratings + rotating polarizer.
- **Key contrast mechanism**: because the two interfering fields traverse nearly the same path, they accumulate the *same* incoherent multiple-scattering background, which cancels in the gradient/difference, while the coherent (structural) contribution differs and survives — this is what makes GLIM work on thick, multiply-scattering samples. Combined with a fully-open condenser aperture (high illumination NA), this gives strong optical sectioning.
- **The critical limitation GODT overcomes**: GLIM's route to 3D is (1) integrate the measured phase gradient back into a phase map via $\phi(x,y)=\int_0^x \nabla_x\phi\,dx' + \phi(0,y)$ (Eq. 3, requiring a known boundary value) and then (2) z-scan + optical-sectioning + spatial high-pass filtering to build a depth stack. This is **not** a rigorous diffraction-tomographic reconstruction — it relies on sectioning rather than inverse-scattering, loses the diffraction-compensation of true ODT, and requires the error-prone gradient-integration step. [[2411-08423-godt]] explicitly positions GODT as removing both of these: feeding the gradient *directly* into Rytov tomographic reconstruction, no integration, no z-scanning.
- Demonstrated: quantitative phase of 4.5μm polystyrene beads (matching diffraction-aware ground truth); time-lapse dry-mass/growth measurement of HeLa cultures over mm-scale fields of view (mass doubling time DT≈36h); dispersion-relation transport analysis computed *directly from the phase gradient* (not the integrated phase); time-resolved tomography of single HeLa cells through mitosis (63×/1.4NA, Δz=0.07μm over 28μm depth).
- Notable methodological point: the intracellular-transport dispersion analysis (Eq. 5-6) was performed on the raw phase gradient $\nabla\phi$ rather than the integrated phase $\phi$ — an early hint that gradient data is directly usable for quantitative analysis without integration, which GODT generalizes to full tomographic reconstruction.

## Method

Add-on module to a commercial DIC microscope: removes the DIC analyzer, Fourier-transforms the two cross-polarized sheared fields (lens L1), applies SLM phase retardation to one polarization at the Fourier plane, re-images (lens L2), and recombines via a 45° polarizer. Low-coherence (broadband) illumination. 4-frame phase shifting per gradient measurement; z-scanning + spatial high-pass filtering for depth sectioning.

## Relevance

Per [[godt]]: the closest methodological predecessor to GODT in the corpus — same shear-direction phase-gradient measurement, same low-coherence-for-SNR philosophy, same thick-sample motivation. The difference is the reconstruction path: GLIM integrates + z-scans (quasi-3D, sectioning-based), GODT skips integration and does Rytov diffraction tomography (3D, diffraction-compensated). Also relevant to the multi-shear extension: GLIM/DIC-family work already established that two orthogonal shear directions recover the full 2D phase gradient (cited here via Cogswell, Shribak), a 2D precedent for the multi-shear idea GODT names as future work in 3D.

## Raises

- GLIM (and DIC generally) uses orthogonal-shear pairs to recover the full transverse gradient in 2D. GODT's core paper does single-shear 3D and names multi-shear as future work — worth checking whether the GLIM/DIC orthogonal-shear 2D reconstruction math offers any reusable structure for combining multiple shear directions in the 3D tomographic (Ewald-space) setting, or whether the 3D case is fundamentally different. Open question, not resolved by any ingested source.
- The epi-illumination variant of GLIM (for opaque samples) is a separate paper, [[kandel2019]].
