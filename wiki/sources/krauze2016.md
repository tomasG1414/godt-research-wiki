---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
source: "https://doi.org/10.1364/OE.24.004924"
---

# Generalized Total Variation Iterative Constraint Strategy in Limited Angle Optical Diffraction Tomography

**Citation:** W. Krauze, P. Makowski, M. Kujawińska, A. Kuś, "Generalized total variation iterative constraint strategy in limited angle optical diffraction tomography," Optics Express 24(5), 4924–4936 (2016). Warsaw University of Technology, Institute of Micromechanics and Photonics.

**Gist:** Introduces GTVIC (Generalized Total Variation Iterative Constraint), a two-stage regularization add-on for limited-angle ODT (LAT-ODT) reconstruction: a first stage generates a 3D object-shape mask via TV-minimization, which is then fed as an a-priori geometric constraint into *any* iterative tomographic solver (tested on both a diffraction-blind, straight-ray solver and a proper Rytov/diffraction-aware Fourier solver). Demonstrates substantial reduction of the characteristic axial elongation/blur of limited-angle ODT reconstructions on both synthetic cell-phantom data and a real fixed C2C12 muscle cell.

## Key results

- Confirms the mechanism behind [[missing-cone-problem]] directly in the LAT-ODT context: because illumination is confined to a cone (here 45°/48° inclination from a stationary sample, not full 360° rotation), the Fourier diffraction theorem leaves regions of the 3D object spectrum unfilled, producing systematically anisotropic resolution (good in-plane, poor axial) and z-elongated reconstructions — exactly what [[2411-08423-godt]] separately reports as blur relative to the ideal finite-difference derivative.
- GTVIC mechanism: (1) run a preliminary TV-regularized reconstruction (Chambolle-Pock primal-dual algorithm minimizing $\|Ax-B\|^2 + \lambda\|x\|_{TV}$) to get an approximate but geometrically-clean object silhouette even though inner structure is corrupted by TV's piecewise-constant bias; (2) Otsu-threshold and morphologically dilate this into a binary 3D mask; (3) pass that mask as a nonnegativity + external-geometry constraint into the real iterative solver (either GDRA — a straight-ray, diffraction-blind algebraic solver, or FTRA — a proper Rytov-approximation Fourier-inversion solver, i.e. the diffractive analog of [[optical-diffraction-tomography|DI]] run iteratively).
- Quantitative gains on synthetic Paramecium-like cell phantom (NA=1.3, λ=0.633μm, 45° conical illumination, 180 projections): GDRA+GTVIC converges in 9 iterations (13 min) vs. 96 iterations (77 min) for GDRA alone — ~6× speedup — plus +2 percentage points in image quality index (QI). FTRA+GTVIC needed only 5 iterations (7 min) vs. 27 for FTRA alone, with markedly reduced axial distortion artifacts, though at ~4.5× more total compute due to the mask-generation overhead (mask stage is ~88% of FTRA+GTVIC's total time vs. only 45% for the already-slower GDRA+GTVIC).
- **FTRA (diffraction-aware, Rytov-based) consistently outperforms GDRA (diffraction-blind, straight-ray)** in QI even with GTVIC applied to both: the paper attributes 0.29 % of the gain to the change of solver alone and a further 0.44 % to using amplitude as well as phase (≈0.73 % combined). FTRA on a CPU also runs ~25× faster than GDRA on a GPU. This directly supports treating the diffraction-aware Rytov/Fourier reconstruction (the [[optical-diffraction-tomography|DI]] family GODT also uses) as the better base solver, independent of any regularization added on top.
- **Amplitude information matters quantitatively, confirming a point [[2411-08423-godt]] only demonstrates qualitatively**: reconstructing from Rytov data with amplitude patterns *equalized to unity* (phase-only, exactly GODT's practical default) loses 0.44 percentage points of QI relative to using true complex (amplitude+phase) projections — "the amplitude patterns encode most of the object information carried by the scattered field" beyond what phase alone gives. This is a directly relevant, quantified cost of the phase-only assumption GODT makes when amplitude-derivative data isn't separately captured.
- Few-view study: best FTRA reconstruction quality plateaus around 60–120 projections; below that, quality is highly sensitive to exact projection-subset choice (erratic, non-monotonic quality vs. count) because the iterative solver settles into different local optima for slightly different input — a caution for any future work assuming "more projections monotonically improves quality."
- Experimental validation on a real fixed C2C12 muscle cell (Mach-Zehnder tomographic phase microscope, opto-electronic angular scanning via SLM, 360 projections, 48° cone) reproduces the synthetic-data findings: GTVIC recovers sharp structural detail (e.g. a basal cavity) without the high-frequency noise that many iterations of plain FTRA introduce.

## Method

Two solvers compared throughout: GDRA (Generalized Data Replenishment Algorithm — SIRT3D-based, straight-ray/non-diffractive, real-valued unwrapped-phase input) and FTRA (Fourier Transform Reconstruction Algorithm — nearest-neighbor-interpolation Rytov/Fourier-diffraction-theorem solver, complex amplitude+phase input, run iteratively rather than as one-shot DI). GTVIC mask built via Chambolle-Pock TV minimization + Otsu thresholding + morphological dilation. Synthetic ground truth: layered Paramecium-like cell phantom (cytoplasm n=1.359–1.371, nucleus n=1.367–1.376, vacuole n=1.332, in water) under conical illumination.

## Relevance

Per [[godt]]: two takeaways. (1) Diffraction-aware (Rytov/Fourier) solvers beat straight-ray solvers regardless of regularization, consistent with GODT's Rytov-based approach. (2) The amplitude-vs-phase-only gap (0.44 QI points) quantifies the cost of a phase-only assumption, which GODT also makes in practice.

## Raises

