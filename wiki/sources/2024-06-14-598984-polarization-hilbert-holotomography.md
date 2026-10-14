---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-07-18
status: current
source: "https://doi.org/10.1002/lpor.202401474"
---

# Polarization Gratings Aided Common-Path Hilbert Holotomography for High-Throughput Lipid Droplets Content Assay

**Citation:** P. Zdańkowski, J. Winnik, M. Rogalski, M. J. Marzejon, E. Wdowiak, W. Dudka, M. Józwik, M. Trusiak, "Polarization gratings aided common-path Hilbert holotomography for label-free lipid droplets content assay," ***Laser & Photonics Reviews*, DOI 10.1002/lpor.202401474 (published 2025-04-13)**. Warsaw University of Technology / EMBL Heidelberg. Cited as ref. 83 in [[2411-08423-godt]]. Ingested from the bioRxiv preprint (doi:10.1101/2024.06.14.598984).

> **Cite the published version.** The ingested copy was the preprint; the peer-reviewed version appeared in *Laser & Photonics Reviews* (found via the [[forward-citation-sweep-2026-07|July 2026 citation sweep]]). Note the title changed from "high-throughput" to "label-free." Content differences between the two versions have not been checked.

**Gist:** Introduces CP-HHT (common-path Hilbert holotomography) — the **same polarization-grating (PG) common-path shearing hardware line** as [[2411-08423-godt|GODT]] (same authors, Zdańkowski/Winnik/Trusiak), but run in **total-shear mode** with a novel single-shot phase-demodulation method (Hilbert spiral transform, HST) instead of GODT's gradient/small-shear mode. Validated on a 3D-printed cuttlefish-brain phantom and applied to quantify lipid-droplet (LD) refractive-index differences (cholesterol-ester-rich vs. triacylglycerol-rich) in HeLa cells.

## Key results

- **Hardware confirmed shared with GODT**: two polarization gratings (PG1, PG2) generate two opposite-circular-polarization beam copies from linearly-polarized input; shear set by grating separation; phase-shifting achieved by rotating a linear polarizer before the camera (rotation α → 2α phase shift) — this is exactly [[common-path-shearing-interferometry|GODT's shearing module]], confirming it is a reusable hardware platform the group runs in either total-shear (this paper, CP-HHT) or gradient/small-shear (GODT) configuration by adjusting grating separation and demodulation method.
- **Space-bandwidth product (SBP) problem this paper solves**: oblique conical illumination (needed for tomography) shifts the hologram's spectral orders relative to each other at every azimuth, so a single fixed Fourier-filter size can't cleanly separate the object term from its conjugate/autocorrelation terms across all illumination angles — forcing a choice between smaller FT filters (losing resolution) or higher magnification (losing field of view). Hilbert spiral transform (HST) phase retrieval tolerates *overlapping* spectral orders (unlike the FT method), recovering comparable image quality to 4-frame temporal phase shifting (TPS) in a **single shot** — a 4× acquisition speedup with no measurable SBP penalty, validated against TPS and FT side-by-side on the brain phantom.
- Reconstruction pipeline: same [[optical-diffraction-tomography|direct inversion (DI)]] + first-order [[rytov-approximation|Rytov]] + [[wolf1969|Wolf's projection theorem]] chain as GODT, citing the identical lineage ([[pan1983]], [[chen1998]], [[wolf1969]], [[krauze2016]] for missing-cone). **Cites and uses the open-source EWALD reconstruction codebase** (github.com/biopto/EWALD), a named reference implementation of the group's DI reconstruction pipeline.
- Total-shear LD application: HeLa cells with cholesterol-ester (CE) supplementation show mean LD RI 1.479 (SD 4.91×10⁻³, wider RI range 0.036) vs. oleic-acid (OA/triacylglycerol) supplementation mean RI 1.468 (SD 2.28×10⁻³, narrower range 0.022) — a clean demonstration that RI tomography alone, without fluorescent labels, can distinguish lipid-droplet biochemical content, statistically confirmed (D'Agostino, Anscombe-Glynn, Bartlett, Wilcoxon tests, all p<0.001).
- Same conical illumination geometry as GODT (58° tilt here vs. 60° in the core paper), same 180-projection/2° azimuthal step scanning pattern, same limited-angle/[[missing-cone-problem]] caveat explicitly acknowledged (citing [[krauze2016]] directly).

## Method

Supercontinuum white-light laser source, multimode-fiber-coupled, conical oblique illumination via rotating-mirror arm; PG1/PG2 shearing module (as in GODT) followed by a linear polarizer and CMOS camera. Two demodulation modes compared: TPS (4 phase-shifted frames, quasi-inline configuration, highest quality) vs. HST (single frame, slightly off-axis low-carrier configuration, comparable quality, 4× faster acquisition) vs. FT (single frame, standard off-axis, but SBP-limited under high-NA conical illumination). Brain phantom: cuttlefish-brain histology stack (44 slices) converted to laser-power-density values and 3D-printed via two-photon polymerization (same phantom-fabrication technique family as [[ziemczonok2019]]).

## Relevance

Per [[godt]]: this shows GODT's shearing hardware is one configuration of a shared, actively developed platform — the same PG shearing module, the same DI/Rytov/Wolf reconstruction chain (with a named open-source implementation, EWALD), and the same conical-illumination/missing-cone limitations, run at large shear (total-shear CP-HHT here) vs. small shear (GODT). For the multi-shear extension, two takeaways: (1) HST is a single-shot alternative to 5-frame TPS demodulation, which would keep the acquisition time of several shear channels per illumination angle tractable; (2) the EWALD codebase is a possible base implementation for reconstruction experiments on real data.

## Raises

- EWALD (github.com/biopto/EWALD) is the cited reconstruction codebase, a natural base for any experiment on real data.
- HST's single-shot capability is demonstrated here for total-shear phase retrieval; whether it carries over to GODT's small-shear gradient regime has not been reported.
