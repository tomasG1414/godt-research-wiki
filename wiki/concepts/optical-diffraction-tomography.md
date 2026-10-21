---
type: concept
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Optical Diffraction Tomography (ODT)

Reconstructs a sample's full 3D refractive-index (RI) distribution $n(x,y,z)$ from a set of complex object waves $u=Ae^{i\varphi}$ captured under multiple illumination directions $k_i=[k_{ix},k_{iy},k_{iz}]$ (or, in some geometries, multiple viewing angles). Decouples RI from physical thickness, unlike 2D transmission [[quantitative-phase-imaging]].

## Theory (Wolf's projection theorem)

Originally derived by [[wolf1969]] under the first Born approximation for a scattering potential $F(R)=-k_0^2[n^2(R)-1]$: the 2D Fourier transform of a scattered field measured in a plane equals a scaled slice of the 3D Fourier transform of $F$, restricted to a spherical surface set by momentum transfer $\mathbf{K}=\mathbf{k}-\mathbf{k}_0$ (formally identical to the von Laue/Bragg equations of X-ray diffraction). Only *homogeneous* (propagating) wave components are accessible — the root cause of the [[missing-cone-problem]]. The modern ODT literature (including [[2411-08423-godt]]) restates this using the [[rytov-approximation]] instead of Born, since Rytov tolerates a wider weak-scattering regime; see [[chen1998]] for the Born-vs-Rytov validity comparison.

Under the first-order Born or [[rytov-approximation|Rytov approximation]] (weakly scattering samples), each captured object wave maps onto a spherical cap — the **Ewald sphere** — in the 3D Fourier spectrum $\tilde{O}$ of the object's scattering potential $O$:

$$\tilde{u}_R(K_x,K_y) = \frac{k_0^2}{2ik_z}\tilde{O}(K_x,K_y,K_z)$$

where $u_R = \ln(u/u_i)$ is the Rytov-transformed field, $(K_x,K_y)=(k_x-k_{ix},\,k_y-k_{iy})$, $K_z=k_z-k_{iz}$, and $k_z=\sqrt{k_0^2-k_x^2-k_y^2}$, with $k_0=2\pi n_0/\lambda$ the wavenumber in the medium (the GODT paper's convention). Varying illumination direction sweeps the cap across different regions of $\tilde{O}$; enough diversified views fill the 3D spectrum. Once filled, inverse FFT gives $O(x,y,z)$, and RI follows from $n = n_0\sqrt{1-O}$ (or the linearized form $n \approx n_0 - \tfrac{n_0}{2}O$).

- **Direct inversion (DI) / direct interpolation**: reposit each measured spectral cap into $\tilde{O}$ via interpolation, then inverse-FFT. [[pan1983]] compares interpolation choices: nearest-neighbor is cheapest but lowest quality (ring artifacts); **bilinear interpolation, after 8× pre-densifying the (K,φ)-plane sampling via zero-padded inverse-FFT, matches filtered-backpropagation's image quality at ~1/14th the compute cost** — the more expensive circular-sampling-theorem interpolation gives no further improvement. [[2411-08423-godt]] describes DI as most commonly realised with nearest-neighbour interpolation, the usual practical choice.
- **Filtered backpropagation (FBP)**: space-domain analog of DI; diffractive generalization of nondiffractive filtered backprojection (CT), due to Devaney (1982). Conceptually elegant but needs ~$N^2$ FFTs (vs. ~4N for bilinear DI) for an N×N reconstruction because filtering is depth-dependent. [[pan1983]]'s overall recommendation: prefer DI (bilinear) over FBP whenever the weak-scattering approximation itself holds — algorithm choice here is a cost decision, not an accuracy one.
- Because illumination/collection angle is bounded by the numerical aperture ($k_x^2+k_y^2 < k_{max}^2$, $k_{max}=2\pi NA/\lambda$), the 3D spectrum is never fully covered — see [[missing-cone-problem]].

## How densely a conical scan samples each voxel

Derived in [[multishear-fusion-design-solution]] for the illumination-scanning geometry, and useful well beyond the multi-shear question. For a conical scan at fixed inclination $\theta$ over $N_\text{fr}$ azimuths, the caps through a voxel $\mathbf{K}$ are the roots of $\cos(\phi-\psi)=c(\mathbf{K})$, where $\psi$ is the azimuth of $\mathbf{K}_\perp$ and

$$c(\mathbf{K})=-\frac{|\mathbf{K}|^2+2k_mK_z\cos\theta}{2k_mK_\perp\sin\theta}$$

Two caps pass through the voxel where $|c|<1$ and none where $|c|>1$, so **the sampled support is exactly $\{|c|<1\}$** — the [[missing-cone-problem]] is a consequence, not a separate fact. Weighting each sample by the squared Rytov factor $|k_0^2/2ik_z|^2$ gives the information density

$$M(\mathbf{K})=\frac{N_\text{fr}\,\delta k}{\pi}\cdot\frac{k_m}{k_z\,\sin\theta\;K_\perp\sqrt{1-c^2}}$$

per voxel of edge $\delta k$. Verified against direct binning of digital caps to sub-percent on totals. Two structural points fall out. The **caustic** $|c|=1$ — the envelope of the cap family — carries an integrable $1/\sqrt{\cdot}$ pile-up of samples, so voxels near the support boundary are oversampled rather than starved. And because the scattered $k_z=K_z+k_m\cos\theta$ depends on $K_z$ alone under a conical scan, the Rytov weight is independent of *which* azimuth delivered the sample, which is what makes the density separable from anything applied in the transverse plane.

Practical consequence for a GODT-like acquisition: **180 projections visit 76.6 % of the support the geometry allows** on a 25.6 μm field under nearest-neighbour gridding (88.4 % on a 12.8 μm field) in this project's simulations. Angular sampling is a budget separate from the missing cone, and it is one argument for [[pan1983]]'s bilinear interpolation.

## Configurations

- **Illumination scanning** (most common in biomedical ODT): stationary sample, scanned illumination angle. Used by [[2411-08423-godt]] and most sources in this project.
- Object/detector scanning variants exist but are less relevant here.
- **Common-path** ODT: object beam interferes with a shifted copy of itself rather than a separate reference arm — see [[common-path-shearing-interferometry]].

## Relevance

Per [[godt]]: this is the theoretical backbone GODT builds on. [[2411-08423-godt]]'s central result (Eq. 19) is that feeding shearing-interferometry phase-gradient data through the *exact same* DI/Ewald-sphere pipeline yields the derivative of $\tilde{O}$ in the shear direction, rather than $\tilde{O}$ itself — GODT changes the input, not the reconstruction algorithm.
