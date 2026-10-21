---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
source: "https://doi.org/10.1038/s41467-019-12634-3"
---

# Epi-Illumination Gradient Light Interference Microscopy for Imaging Opaque Structures

**Citation:** M. E. Kandel, C. Hu, G. Naseri Kouzehgarani, E. Min, K. M. Sullivan, H. Kong, J. M. Li, D. N. Robson, M. U. Gillette, C. Best-Popescu, G. Popescu, "Epi-illumination gradient light interference microscopy for imaging opaque structures," Nature Communications 10, 4691 (2019). University of Illinois at Urbana-Champaign (Popescu group). Cited as ref. 50 in [[2411-08423-godt]].

**Gist:** Reflection-geometry version of [[nguyen2017|GLIM]] (epi-GLIM), extending phase-gradient imaging to bulk and opaque samples that transmission GLIM cannot reach — silicon wafers, cells on opaque/plastic substrates, intact mouse tendon, live zebrafish. Reaches 3D via an empirically-measured PSF (using co-localized microspheres as fiducial markers) plus Wiener-regularized deconvolution, rather than analytic tomographic inversion.

## Key results

- Same core gradient measurement as [[nguyen2017]]: DIC/Nomarski-sheared beams, 4-step phase shifting (here via a liquid-crystal variable retarder rather than an SLM), yielding $\Delta\phi = \phi(r+\delta r)-\phi(r) \approx (\nabla_x\phi)\delta r$. Pixelwise subtraction across the 4 phase-shifted frames removes the non-modulating incoherent multiple-scattering background — the same background-rejection mechanism that makes the whole GLIM family work on thick samples.
- **Reflection doubles phase sensitivity**: on a reflective substrate the measured phase shift is twice the transmission value ($\phi = 2\cdot\frac{2\pi d}{\lambda}(n-n_0)$), validated on 1.9 μm beads (n=1.605 in n₀=1.518 oil at λ=490nm).
- **Frequency-coverage advantage relevant to [[missing-cone-problem]]**: with a fully open condenser (using the objective itself as condenser, hence higher effective NA), epi-GLIM's frequency coverage is the *autocorrelation of the objective pupil function* — meaning "an open aperture simultaneously captures the frequency coverage that is synthesized in other instruments through multiple measurements." This is a genuinely different route to filling k-space than the illumination-angle-scanning approach [[2411-08423-godt|GODT]] uses, trading angular scanning for aperture-synthesized coverage.
- **3D reconstruction is empirical, not analytic**: the authors explicitly note that the edge-response method for measuring the PSF cannot support 3D reconstruction because "the shear in GLIM renders the transfer function non-radially symmetric" — the same orientation-dependence that gives GODT its single-shear contrast anisotropy. Their workaround is a fiducial-object method: 3μm polystyrene microspheres are co-localized with the sample (mixed into the fixation solution) so they share optical path and aberrations; since the sphere's scattering potential $x_0(r)$ is known and its measured signal $y_0(r)$ is observed, $y_0 = \text{PSF} \circledast x_0$ gives the PSF by deconvolution, and an unknown object is then recovered by Wiener-regularized deconvolution (Eq. 4) that merges both steps to avoid two separate regularizers.
- Demonstrated across samples impossible in transmission: 45 cm² silicon-wafer mosaic with nanoscale topography extraction, live cells in plastic microplates, spheroid cultures with per-spheroid thickness quantification, tendon in an intact mouse leg, and time-resolved tomography of a whole live zebrafish.

## Method

Add-on module to a reflection Nomarski (DIC) microscope; liquid-crystal variable retarder for 4-step phase shifting; broadband illumination for coherence sectioning; in-house automation for x/y/z scanning and mosaicking. 3D reconstruction via fiducial-microsphere-derived PSF + Wiener deconvolution. Gradient integration (for quantitative phase) by cumulative sum along the shear direction using a separately-measured shear distance δr = 0.3 μm.

## Relevance

Per [[godt]]: the epi-mode companion to [[nguyen2017]] in the [[gradient-phase-imaging]] family, and the closest thing in the ingested corpus to an *acknowledgment of GODT's core problem from the other side*. Two points bear directly on the multi-shear extension:

1. **The shear-induced transfer-function asymmetry is explicitly named here as an obstacle** ("the shear in GLIM renders the transfer function non-radially symmetric"), which is the same orientation-dependent contrast [[2411-08423-godt]] flags as its known single-shear limitation. The GLIM lineage responds by sidestepping analytic inversion entirely (measure the PSF empirically with fiducials); GODT responds by deriving the exact $iK_x$ relation instead. A multi-shear GODT extension is a third response — combining shears to restore symmetry — and neither predecessor took that route.
2. **The aperture-synthesis idea** (open condenser giving pupil-autocorrelation frequency coverage in one shot, replacing multi-measurement synthesis) is an alternative axis for improving k-space coverage that is orthogonal to both illumination-angle scanning and shear-direction diversity. Whether it composes with GODT's approach is unexamined *(inference)*.

## Raises

- Epi-GLIM's fiducial-PSF/Wiener-deconvolution reconstruction is empirical and system-specific, unlike GODT's analytic Rytov inversion. It sidesteps rather than solves the shear-asymmetry problem — worth noting as evidence that the asymmetry is a real, independently-encountered obstacle, not an artifact of GODT's particular formulation.
