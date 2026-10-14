---
type: concept
projects: [godt]
created: 2026-09-21
updated: 2026-09-21
status: current
---

# Optimal Experimental Design

Choosing *which measurements to take*, and how to split a fixed budget between them,
by optimising a scalar function of the Fisher information they deliver. Standard
statistical machinery (general knowledge, not sourced in this wiki); it enters this
project as the frame for "which shears, for how long" in
[[multishear-fusion-design-problem]] and is applied in
[[multishear-fusion-design-solution]].

## The structure that makes it work

A design is a vector of time shares $\mathbf t$ on the simplex
$\{t_j\ge0,\ \sum_jt_j=1\}$ over a candidate set of channels. When each channel
contributes additively to the information,

$$J(\mathbf K)=\sum_j t_j\,a_j(\mathbf K)$$

$J$ is **linear in $\mathbf t$**, and the three standard objectives are convex:

| criterion | objective | what it buys |
|---|---|---|
| **A-optimal** | minimise $\int w\,(J+S^{-1})^{-1}$ | mean error |
| **D-optimal** | maximise $\int w\log(1+SJ)$ | log-information, i.e. channel capacity |
| **E- / maximin** | maximise $\min_{\mathbf K} J$ | worst case; an LP after discretisation |

The Bayesian versions add the prior term $S^{-1}$, which keeps the objective finite
where $J=0$ — necessary here, since shear nulls make $J$ vanish on whole planes.

## Certifying an optimum

The reason to prefer this frame over a parameter sweep is that the answer comes with
a proof. By the general equivalence theorem, a design $\mathbf t$ minimising a convex
$f$ over the simplex is optimal **iff no vertex is a descent direction**:

$$\min_j g_j(\mathbf t)\ \ge\ \sum_j t_j\,g_j(\mathbf t),\qquad g_j=\partial f/\partial t_j$$

The normalised gap $(\sum_jt_jg_j-\min_jg_j)/|f|$ is zero at the optimum and
computable at any candidate design, so a reported optimum is verifiable rather than
asserted. For the A-criterion this is the Kiefer–Wolfowitz-type statement on the
per-channel sensitivity $\phi_j=\sum_{\mathbf K}w\,a_j/(J+S^{-1})^2$. `lab/shearlib.py`
solves by SLSQP on the full candidate set followed by Fedorov–Wynn-style
prune-and-add rounds, and reports the gap alongside every design.

## Sparsity is not automatic

Carathéodory-type results are often quoted as implying optimal designs sit on few
support points. That bound applies to the *number of distinct information matrices
needed to span a convex hull*, and it does not force sparsity when many candidates
deliver nearly interchangeable information. In the multi-shear case the optimum uses
24 of 64 candidates and stays spread across three decades of SNR — the conjecture
that it would be sparse was **refuted** ([[multishear-fusion-design-solution]]).

What did sparsify it was a *physical* penalty on the candidates: falling fringe
visibility with shear magnitude collapsed the optimum from 24 channels to 6. The
general lesson worth carrying: sparsity in an optimal design usually reflects a cost
that distinguishes the candidates, not the geometry of the information itself.

## Robustness by concave substitution

Model error that scales with the measurement — here shear miscalibration — enters as
a signal-dependent noise floor,

$$J_\text{rob}=\sum_j\frac{t_j\,a_j}{1+t_j\,c_j},$$

each term of the form $t/(1+ct)$, which is concave and increasing. So $J_\text{rob}$
stays concave in $\mathbf t$ and the A- and D-problems stay convex: robustness is
free, structurally. The qualitative effect is a ceiling on what any single channel
can deliver, which pushes the optimum to spread across channels.

## Relevance

Per [[godt]]: this is the frame that turns "how many shears, and which" from a
parameter sweep into a solved problem with a certificate. Its main payoff in this
project was negative and therefore useful — it established that a simple rule
($\Delta\approx2\pi/K_{\max}$, one orthogonal pair, two incommensurate magnitudes)
is within about 10 % of the certified optimum, so the optimisation is best presented
as the thing that *justifies* the simple rule rather than as the deliverable.
See also [[partially-coherent-transfer-function]] for the visibility cost that drives
the sparsity result.
