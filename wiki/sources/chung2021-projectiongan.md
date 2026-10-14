---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-07-22
status: current
source: "https://arxiv.org/abs/2103.09022"
---

# Unsupervised Missing Cone Deep Learning in Optical Diffraction Tomography (ProjectionGAN)

**Citation:** Hyungjin Chung*, Jaeyoung Huh*, Geon Kim, Yong Keun Park, Jong Chul Ye, arXiv:2103.09022v1 [eess.IV], 16 March 2021 (*equal contribution). KAIST.

**Gist:** Attacks the [[missing-cone-problem]] with an unsupervised CycleGAN that operates on **2D projections** rather than the 3D volume. The trick that makes it cheap: after a Gerchberg–Papoulis reconstruction, project the volume to 2D at all angles; projections at measured angles are sharp, projections at missing angles are blurry, and — absent a missing cone — both would be drawn from the same distribution. So learn the optimal-transport map between those two distributions. No matched ground truth needed, and 27 s per reconstruction against 605 s for TV. From Ye's group, the same lab as [[lim2015]].

## Key results

**Quantitative RI recovery on a silica microbead of true RI 1.46** (Fig. 7c line-profile means): Rytov 1.371, TV 1.388, GP 1.398, ProjectionGAN 1.424. Best of the four but still 2.5 % low — the missing cone is mitigated, not eliminated, whatever the k-space figures suggest.

**The k-space claim, stated plainly.** In Rytov/GP/TV reconstructions the 49 individual Ewald half-sphere sampling curves remain visible in the spectrum, and the missing cone is still evident in the $y$–$z$ and $x$–$z$ planes even after GP. After ProjectionGAN those curves are "no longer visible" and the cone is filled in every viewing plane. Read carefully, this is the network *synthesising* plausible spectral content, not measuring it — see Relevance.

**Reconstruction times** (GTX 2080-Ti, Table I): Rytov 1.37 s, GP 14.08 s, ProjectionGAN 27.33 s, TV 604.86 s, CryoGAN order of hours.

**Confirms the standard regularizer critiques on real data.** TV produces "cartoon-like appearances, which do not capture the true texture of the cell"; GP resolves elongation partly but leaves the support vague and RI underestimated. Both match what [[lim2015]] and [[krauze2016]] report — useful independent corroboration on NIH3T3 cells rather than phantoms.

## Method

GP reconstruction to initialise the volume → parallel-beam forward projection (a deliberate straight-ray approximation, justified via the projection-slice/projection-diffraction relation, chosen for speed) → CycleGAN enhancement of the projections → filtered backprojection at 360 equiangular views around each of the three coordinate axes, with the three resulting volumes averaged. CycleGAN is formulated as Kantorovich optimal transport with a Wasserstein-1 statistical distance, LS-GAN discriminator loss for the Lipschitz condition; U-Net generators, patchGAN discriminators, optional 4-level multi-scale discriminator when the sample has multiple particles. Trained on 67 NIH3T3 cells + 10 microbeads (of 87 + 18 acquired), and separately on 20 of 27 numerical sphere phantoms; ~1 day of training each, PyTorch, Adam.

## Relevance

Per [[godt]]: two points of context.

- **Its FBP back-end gives up diffraction compensation.** Parallel-beam projection plus filtered backprojection is straight-ray tomography — precisely the diffraction-uncompensated approach that [[optical-diffraction-tomography]] exists to replace, and the same weakness this wiki notes in [[vishnyakov2016]]. The paper is candid that this is a fast approximation, but it means the final volume is not diffraction-compensated even though the initialisation was.
- **Trust boundary.** A generative model filling the cone produces spectral content that was never measured. For quantitative RI that is a real tension, and the residual 2.5 % bead error is a fair measure of it.

## Raises

- Is there a published version? This is arXiv v1 (2021); given [[forward-citation-sweep-2026-07]] already caught two preprint-to-published corrections in this corpus, the citation should be checked before use.
- Sim, Oh, Kim, Jung & Ye, *SIAM J. Imaging Sci.* **13**(4), 2281 (2020) — the optimal-transport CycleGAN theory this rests on.
