---
type: source
projects: [godt]
created: 2026-09-11
updated: 2026-09-11
status: current
source: "https://doi.org/10.1364/OE.15.017805"
---

# Lateral shearing interferometry: theoretical limits with practical consequences

**Citation:** M. Servín, M. Cywiak, A. Dávila, *Optics Express* **15**(26), 17805–17818 (2007). DOI: [10.1364/OE.15.017805](https://doi.org/10.1364/OE.15.017805). Centro de Investigaciones en Óptica, León, Mexico.
**Gist:** The explicit Fourier-domain analysis of the lateral shearing operator's frequency response, its inverse, and — the paper's self-stated main contribution — the spatial support (topology) of the recovered phase from bounded pupils using one, two, or $n$ sheared interferograms. This is the paper [[zhu2024-appl-phys-b]] cites as the standard reference for the frequency response of the shearing operator.

## Key results

### 1. Forward shearing filter and its nulls

The LSI forward operator in the Fourier domain (Eq. 6):

$$F[\Delta^x(x,y)] = 2i\sin(dx\,\omega_x)\,\Phi(\omega_x,\omega_y)$$

where $dx$ is the half-shear (total displacement $2dx$). Nulls at $\omega_x = m\pi/dx$ for integer $m$. The "direct" inverse filter (Eq. 7):

$$H^{-1}(\omega_x) = \frac{-i}{2\sin(dx\,\omega_x)}$$

has poles at those same frequencies and "cannot be implemented in a direct way" (their words). Fig. 2 plots the comb of poles.

This is exactly the same structure as [[s42005-025-02416-3]]'s $H_\Delta(K) = (e^{j\Delta K_x}-1)/\Delta$, rewritten: $|H| = 2|\sin(K\Delta/2)|/\Delta$, nulls at $K = 2\pi m/\Delta$. The mapping is $\omega_x = K$, $dx = \Delta/2$.

### 2. Regularized least-squares inverse

To handle the poles, they minimise a LS cost functional with a first-order difference regularisation (Eq. 10), yielding a stable filter (Eq. 11):

$$H^{-1}_{reg}(\omega_x) = \frac{i\sin(dx\,\omega_x)}{1+\eta - \cos(2dx\,\omega_x) - \eta\cos(\omega_x)}$$

The regularisation parameter $\eta < 1$ trades off stability vs. frequency content. The filter has no poles but is strongly low-pass (Fig. 3) and **cannot recover information at the forward operator's zeros** — stated explicitly.

### 3. Spatial support of the recovered phase (main contribution)

For a bounded pupil $P(x,y)$:

- **Single shear** (§4, Eq. 13): $\mathrm{supp}\,\hat\phi = P(x,y) \cap [P(x+2dx,y) \cup P(x-2dx,y)]$ — strictly smaller than $P$, losing a strip at the boundaries. Fig. 4 illustrates.
- **Two orthogonal shears** (§5, Eq. 15): $\mathrm{supp}\,\hat\phi = P \cap [P(x\pm2dx,y) \cup P(x,y\pm2dy)] \subseteq P$ — can cover the **entire pupil** even for large shears and irregular pupils.
- **Continuity** (§6): requires $\phi \in S$, a function space where $L_x^{-1}[L_x(\phi)] = \phi$ and $L_y^{-1}[L_y(\phi)] = \phi$. Notably, defocusing ($\phi = D(x^2+y^2)$) is **not** in $S$: each direction recovers only its own quadratic, and the correct phase exists only in the overlap region.

### 4. $n$-direction generalisation (§7)

Generalises to $n$ shearing directions $\mathbf{d}_i = |\mathbf{d}_i|\angle\alpha_i$ with arbitrary magnitudes and angles. The regularised cost functional (unnumbered, p. 17817) and spatial support are given. The paper notes that for $n > 2$ the support "will probably cover the whole two-dimensional wavefront pupil $P(\mathbf{r})$ even if large amounts of shearing $\mathbf{d}_i$ are used."

Cites Villa *et al.* [8] and Paez *et al.* [1] for experimentally efficient one-sided shearing interferometers with variable magnitude and direction.

## Method

Analytical Fourier analysis + regularised least-squares inverse operator. No experimental data in this paper; two synthetic wavefront reconstructions on irregular pupils (Figs. 7–10) demonstrate the spatial support theory. Phase estimation error 0.14–0.2 rad on 24-rad-amplitude test functions.

## Relevance

### For [[godt]]

This paper is the reference for prior art on the transfer nulls. The **null comb itself** — $\sin(\omega\,dx) = 0$ poles in the inverse shearing filter — is known, stated explicitly, and dated to 2007 (and implicitly earlier, since the authors call the frequency response "well known" in the abstract). The regularised inverse's inability to recover the lost information is stated. Its scope is 2D wavefront reconstruction with a focus on spatial support; [[finite-shear-null-experiment]] carries the null analysis into the 3D tomographic setting.

One further point: the $n$-direction regularised LS framework is a plausible template for multi-shear tomographic fusion if the forward model is swapped from the 2D shearing operator to the 3D $H_\Delta(K)$ Ewald relation. Ref. [20] is Servín, Malacara & Marroquín (1996), the two-orthogonal-shear recovery paper this one builds on.

## Raises

- The companion paper Servín, Malacara & Marroquín, *Appl. Opt.* **35**(22), 4343 (1996) [ref. 20] is the two-orthogonal-shear recovery this builds on.
- The $n$-shear LS framework might template the 3D fusion cost functional. Worth revisiting when designing the reconstruction half of the extension.
