---
type: concept
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Quantitative Phase Imaging (QPI)

2D imaging contrast mode that maps the optical phase delay a transparent sample introduces (locally varying optical thickness = RI × physical thickness, convolved together and not separable in 2D). Predecessors were qualitative: Zernike phase contrast (1942) and Nomarski differential interference contrast (DIC) microscopy convert phase to visible intensity but don't yield interpretable values. QPI proper recovers actual phase-delay values, typically via digital holographic microscopy (DHM) and Fourier-transform fringe analysis.

Full 3D structure (decoupling RI from thickness) requires tomographic reconstruction — see [[optical-diffraction-tomography]].

## Coherence tradeoffs

- High-coherence (laser) sources give strong holographic contrast but suffer speckle noise and spurious interference.
- Low-coherence strategies: pseudothermal sources (rotating scattering plate, spatial coherence reduced, temporal stays high), wide-spectral-band sources (superluminescent diodes, LEDs, supercontinuum + AOTF) reducing temporal coherence while keeping spatial coherence high via fiber/pinhole coupling. Reduced coherence length breaks straightforward off-axis hologram recording, pushing toward [[common-path-shearing-interferometry|common-path configurations]].
- Some techniques (FPM, TIE, qDPC) skip hologram recording entirely and numerically solve for phase from fringeless intensity measurements, trading experimental simplicity for numerically ill-posed low-spatial-frequency recovery.

## Relevance

Per [[godt]]: GODT is positioned as extending gradient-family QPI ([[gradient-phase-imaging]]) into a rigorous tomographic (3D, diffraction-compensated) method, while retaining low-coherence-illumination SNR benefits.
