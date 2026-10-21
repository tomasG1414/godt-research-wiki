---
type: note
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Forward citation sweep, July 2026

**Caveat up front: this page is built from OpenAlex metadata and abstracts, not from reading the papers.** Nothing here has been ingested. Titles, venues, DOIs and dates are from the citation graph; where an abstract is quoted or paraphrased it is marked. Treat every entry as a lead to verify, not an established claim.

## Method (reproducible)

Nine anchors resolved to OpenAlex work IDs, then everything citing them published 2024-01-01 onward was pulled and aggregated, so that works citing *multiple* anchors surface first. A second, one-off OpenAlex query (not in the script) listed the GODT authors' 2025+ papers, to catch related work that does not cite the anchors yet.

Re-run with `python lab/citation_sweep.py` (set `OPENALEX_MAILTO`), then `python lab/report.py multi|recent|all|<keyword>`. Results cache to `lab/out/citations.json`. Worth repeating periodically — the 2026 entries below were days old when found.

Anchors and their forward-citation counts since 2024:

| Anchor | OpenAlex | Total cites | Since 2024 |
|---|---|---|---|
| [[wolf1969]] | W2002832520 | 1334 | 137 |
| [[nguyen2017]] | W4240691888 | 269 | 73 |
| [[ziemczonok2019]] | W2995240891 | 77 | 25 |
| [[krauze2016]] | W2277581074 | 65 | 16 |
| [[brenner1993]] | W2158726345 | 75 | 11 |
| [[vishnyakov2016]] | W2474858503 | 19 | 7 |
| [[s42005-025-02416-3]] | W4416280401 | 4 | 4 |
| [[2509-15382-osi-flex]] | W4414694906 | 0 | 0 |
| [[chen1998]] | — | — | **unresolved** |

225 unique citing works; 34 cite two or more anchors. **[[chen1998]] failed to resolve** — its DOI lookup returned nothing and the title fallback matched an unrelated French air-pollution paper, which was discarded. It contributed no works, so the 225 are uncontaminated, but its forward citations remain unswept.

## Immediate corrections to the wiki

Two pages described as preprints have since been published. Both are now fixed:

- **[[2509-15382-osi-flex]]** → *Journal of Physics: Photonics*, DOI `10.1088/2515-7647/ae3504`, published 2026-01-23 (accepted 2026-01-07).
- **[[2024-06-14-598984-polarization-hilbert-holotomography]]** → *Laser & Photonics Reviews*, DOI `10.1002/lpor.202401474`, 2025-04-13.

This is the third instance in this project of a preprint having a published version with a different citation — the same pattern as the GODT preprint. Worth a habit: check for a published version before citing anything from `raw/`.

## Three closely related works

All three come from the GODT authors' wider collaboration and were invisible to a backward-citation-only view.

### 1. R2D-QPI — another route to imaging confluent samples

*"Doubling the field of view and eliminating the replica overlap problem in common-path shearing quantitative phase imaging"* — Rogalski, Zdańkowski, Heto, Mierzejewska, Lenarcik, Li, Sun, Zuo, Trusiak. *J. Phys. Photonics*, `10.1088/2515-7647/ae236c`, 2025-11-24 (arXiv 2510.01762).

Per the abstract: common-path shearing DHM "suffer[s] from sheared object beam (replica) overlap... This limits their use to sparse objects only." R2D-QPI "numerically decouples object and replica fields of view through **controlled shear scanning**", "analytically separates overlapped phase images and effectively doubles the imaged area, **requiring only two measurements**", validated on "highly confluent samples with strong object–replica overlap."

**Why this matters.** GODT's small shear removes the need for an object-free region, which makes confluent samples imageable. R2D-QPI addresses the same limitation from the other direction: keep the large shear and separate the replicas numerically. It is a 2D method, so it complements GODT rather than replacing it. Of interest here: **it works by scanning shear**, acquiring at several shear settings and combining them, which makes it a multi-shear precedent in 2D. Whether its two-measurement separation math has any bearing on multi-shear tomographic fusion is unknown and worth checking.

### 2. LS-ODT — the group's other thick-sample shearing tomography

*"Lateral shearing optical diffraction tomography of brain organoid with reduced spatial coherence"* — Gocłowski, Winnik, Dubey, Zdańkowski, Trusiak, Neogi, Varshney, Ahluwalia. arXiv 2601.20082, 2026-01-27.

Per the abstract: a common-path ODT approach for "heterogeneous and highly scattering thick samples such as tissues and organoids", using "**partial** lateral shearing off-axis interferometry to suppress the effects of multiple scattering, similar to the mechanism in DIC microscopy", plus "dynamic speckle illumination" for sensitivity. Validated on a cell phantom, mouse kidney sections, and iPSC-derived brain organoids.

**Why this matters.** Shares authors with GODT (Winnik, Zdańkowski, Trusiak) and the same multiple-scattering-suppression rationale, but appears to reconstruct RI rather than its gradient, and uses *partial* shear with speckle illumination rather than GODT's small-shear gradient regime. So the group runs at least two shearing-tomography lines, and it is useful to know how they relate.

### 3. Double Sagnac LS-QPM **and tomography**

*"High-Sensitivity, High-Throughput Double Sagnac Lateral Shearing Quantitative Phase Microscopy and Tomography with Pseudo-Thermal Illumination"* — Gocłowski, Trusiak, Ahluwalia, Ahmad. arXiv 2605.04297, 2026-05-05.

Per the abstract: "double Sagnac common-path interferometric configuration combined with pseudo-thermal illumination", enabling "robust single-shot phase retrieval", and explicitly positioned against QWLSI-based QPM which "typically suffer[s] from a trade-off between spatial resolution and field of view."

**Why this matters.** A *double* Sagnac shearing architecture is multi-shear hardware, and the title mentions tomography, so it was the most relevant read of the sweep. It has since been read in full: [[2605-04297-double-sagnac-ls-qpm]].

## Secondary leads, by theme

**Gradient/derivative tomography lineage** (the [[vishnyakov2016]] neighbourhood):
- *Phase contrast tomography (PCT)* — Optics Express, `10.1364/oe.588947`, 2026-02-26 (bioRxiv version 2024-03, `10.1101/2024.03.17.585445`). Cites both [[nguyen2017]] and [[vishnyakov2016]] — same lineage as GODT.
- *Differential phase contrast intensity diffraction tomography* — Optics Letters, `10.1364/ol.531869`, 2024-11, 8 cites. Likely the DPC+IDT work the published GODT cites as ref. 59.
- *Optical diffraction tomography based on quadriwave lateral shearing interferometry* — Yuan et al., Opt. Laser Technol., `10.1016/j.optlastec.2024.111124`, 2024-05, 10 cites. **This is ref. 61 of the GODT paper** — the QWLSI-based ODT that GODT compares itself with.

**Multi-shear / single-shot acquisition** (bears on [[multi-shear-acquisition]]):
- *Single-shot ODT based on polarization multiplexing via Kramers-Kronig relations* — Opt. Laser Technol., `10.1016/j.optlastec.2026.115467`, 2026-05.
- *Single-shot intensity diffraction tomography via polarization-multiplexed LED illumination* — Optics Letters, `10.1364/ol.533430`, 2024-08.
- *Snapshot parallel refractive index tomography* — Optica, `10.1364/optica.577676`, 2026-02.
- *k-space holographic multiplexing for synthetic aperture diffraction tomography* — APL Photonics, `10.1063/5.0203117`, 2024-05, 17 cites.
- *Single-shot quantitative DPC microscope using a single calcite beam displacer* — Applied Optics, `10.1364/ao.536593`, 2024-10.

**Reviews worth skimming for orientation:**
- *Holotomography* — Nature Reviews Methods Primers, `10.1038/s43586-024-00327-1`, 2024-07, **76 cites** — the field's current primer.
- *Optical diffraction tomography for 3D refractive index imaging* — Light: Advanced Manufacturing, `10.37188/lam.2026.077`, 2026 — cites GODT.
- *Quantitative phase imaging based on holography: trends and new perspectives* — Light Sci. Appl., `10.1038/s41377-024-01453-x`, 2024-06, **246 cites**.

**Reconstruction methods:**
- *Inverse-scattering in biological samples via beam-propagation* — bioRxiv `10.1101/2025.08.17.670744` (4 anchors), plus a 2026 absorptive-sample sequel. Relevant to the [[ma2017]] nonlinear-forward-model thread.
- *Segment and support: dual-purpose deep learning for limited angle holographic tomography* — Optics Express, `10.1364/oe.541288`, 2024-11. Warsaw lineage, [[missing-cone-problem]] mitigation.
- *Analytic Fourier ptychotomography* — Nature Communications, `10.1038/s41467-025-67460-7`, 2025-12.

## What this changes

The backward-looking corpus established the lineage GODT came from. The forward picture adds two things: **(a)** the confluent-sample problem has at least one other published solution, and **(b)** shearing tomography with multiple or partial shear and single-shot acquisition is an active area, including among GODT's own co-authors.
