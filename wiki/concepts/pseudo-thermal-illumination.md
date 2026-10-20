---
type: concept
projects: [godt]
created: 2026-07-19
updated: 2026-10-20
status: current
---

# Pseudo-Thermal (Dynamic-Speckle) Illumination

A laser passed through a **rotating diffuser**, producing a source that is spatially incoherent (a moving speckle field, averaged over the exposure) while remaining **temporally coherent** — narrowband, because the diffuser scatters without broadening the spectrum. Also called dynamic speckle illumination.

The point is that spatial and temporal coherence are independent knobs, and interferometric microscopy wants them set differently:

| | spatial coherence | temporal coherence |
|---|---|---|
| **wanted low** — kills speckle, coherent noise, spurious interference from stray reflections | ✓ | |
| **wanted high** — long coherence length, so high-density carrier fringes survive at nonzero OPD | | ✓ |

A broadband LED (what [[s42005-025-02416-3|GODT]] uses) lowers *both*. A pseudo-thermal source lowers only the spatial one, which is the one responsible for speckle. This is why [[2605-04297-double-sagnac-ls-qpm]] can run **single-shot Fourier-transform demodulation** on a dense two-carrier mesh interferogram while still getting 11.7/13.3 mrad spatial phase sensitivity — the fringes need coherence length, the sensitivity needs the speckle gone, and the two requirements do not actually conflict.

Already noted in passing on [[quantitative-phase-imaging]] as one of the two standard low-coherence strategies; this page exists because the distinction turned out to be load-bearing for [[multi-shear-acquisition]].

## Why this matters for multi-shear

[[multi-shear-acquisition]] listed carrier-frequency multiplexing (the QWLSI-style route to encoding several shear channels in one frame) as **closed for GODT**, on the grounds that it requires spatial coherence which GODT's LED does not provide. [[2605-04297-double-sagnac-ls-qpm]] shows that framing was too coarse. What dense carrier fringes actually require is sufficient *temporal* coherence; low spatial coherence is compatible with — and in a common-path shearing geometry positively good for — carrier-based demodulation, because the two interfering beams are sheared copies of the *same* field: their fringe visibility is the degree of coherence at the shear distance, $|\mu(\Delta)|$, which stays high while the shear is within the illumination's transverse coherence length.

So the route is not closed. It is closed *for an LED* and open for a pseudo-thermal source, at the cost of swapping GODT's illumination. Whether that swap is free is the open part — see below.

## What swapping the source would cost GODT

Not established by any ingested source; the following are the questions, not answers.

- **Temporal coherence is not free in a tomographic geometry.** GODT's broadband LED suppresses coherent artifacts that a narrowband source reintroduces, and [[ledwig-robles2023]] reports that broadband illumination *flattens* the 3D transfer function relative to the monochromatic case. A pseudo-thermal source is quasi-monochromatic, so it forfeits any such flattening.
- **One thing it would fix.** GODT reconstructs with a coherent *monochromatic* Rytov model. A quasi-monochromatic pseudo-thermal source matches that model better than a broadband LED does, so the illumination/model mismatch flagged earlier in this wiki partly dissolves — the remaining gap is spatial partial coherence only. *(inference)*
- **The FT-demodulation obstacle is unrelated and survives.** [[common-path-shearing-interferometry]] notes that Fourier-transform demodulation breaks down under high-NA *conical tomographic* illumination, because the filter width needed to capture the full NA conflicts with the oblique-illumination carrier shift. [[2605-04297-double-sagnac-ls-qpm]] never faces this, since it is non-tomographic and illuminates on-axis. Changing the source does nothing about it, so single-shot multi-shear *tomography* is not simply unlocked by adopting pseudo-thermal illumination.

## Elsewhere in the corpus

The Warsaw/UiT shearing lineage uses this source type in its thick-sample work: the LS-ODT brain-organoid preprint (arXiv 2601.20082, per [[forward-citation-sweep-2026-07]]) pairs "dynamic speckle illumination" with partial shear for the same multiple-scattering-suppression reasons. [[2605-04297-double-sagnac-ls-qpm]]'s refs. 18–20 (Ahmad 2016, Ahmad 2022 "Demystifying speckle field interference microscopy", Choi 2014) are the methodological basis.

## Relevance

Per [[godt]]: this is the illumination choice that reopens single-shot carrier multiplexing as a multi-shear acquisition route, at the price of giving up broadband illumination. It changes one row of the strategy table on [[multi-shear-acquisition]] from "incompatible" to "incompatible with an LED, available with a source swap."
