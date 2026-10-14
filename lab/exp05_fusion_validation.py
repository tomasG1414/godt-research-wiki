"""exp05 -- Phantom validation of BLUE/Wiener fusion and of the MSE prediction.

Establishes (design note T4):
  1. Under the linear Rytov model with digital Ewald caps (A1-A7), the
     empirical error spectrum of the fused reconstruction matches the
     predicted (J + 1/S_O)^-1 -- exactly (to Monte Carlo error) for a
     Gaussian random field drawn from the prior, and closely for a textured
     cell phantom with a shell-averaged empirical prior.
  2. At equal photon budget: single 0.8 um shear vs orthogonal pair vs
     {0.8, 0.55} x + y vs the (P-A) optimum vs the ODT benchmark, fused by
     Wiener, by BLUE without prior, and by naive averaging of per-channel
     reconstructions. Real-space RMSE against the band-limited truth and the
     full truth (the missing cone dominates the latter).
  3. The null planes of the single shear are visible as error spikes in the
     K_x profile and are gone after fusion.

Seeded (SEED). Writes: out/validation_spectra.png, out/validation_slices.png
Wiki:   wiki/notes/multishear-fusion-design-solution.md

Run: python exp05_fusion_validation.py     (~3-5 min)
"""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, KGrid, G2_vec, hermitian_sym, cap_density, cap_samples,
                      azimuths, DesignProblem, solve_design, nulls_in_band, out_path,
                      cell_phantom, shell_spectrum, hsym_c, fuse)

SEED = 3
RHO = 100.0            # per-frame-sample SNR at unit |h|, relative to S_ref (see below)
N_ILL = 720            # azimuths actually simulated (dense enough to fill the support)
N_FRAMES = 180         # exposure normalisation: photons of a 180-frame acquisition
K_LO = 2.0             # bandpass for every aggregate metric -- see below
RULE = "=" * 78

# WHY K_LO EXISTS. |G_j|^2 = 4 sin^2(K.d/2) -> 0 as K -> 0 for EVERY shear, so
# the |K| < K_LO shell is lost by any gradient method and its huge prior power
# swamps every aggregate. exp02 hit exactly this and concluded the opposite of
# the truth (lab/README.md, "Gotchas"). All aggregates below are taken on
# supp & (|K| >= K_LO); the low-K loss is reported separately, once.



# ------------------------------------------------------------ forward model
def acquire(cfg, grid, Ot, dvec, phis, sigma2, rng, chunk=32):
    """One channel (shear vector dvec, or None for the ODT benchmark).
    Returns num = sum_i h* y / sigma2 and J = sum_i |h|^2 / sigma2 (3D, not
    symmetrised). h = G(K_perp) k_m/k_z, y = h O~(K) + CN(0, sigma2)."""
    n, nz = grid.n, grid.nz
    G = (np.exp(1j * (grid.KX * dvec[0] + grid.KY * dvec[1])) - 1.0) if dvec is not None \
        else np.ones(grid.KX.shape, dtype=complex)
    Of = Ot.ravel()
    num = np.zeros(grid.size, dtype=complex)
    J = np.zeros(grid.size)
    sig = np.sqrt(sigma2 / 2)
    for start in range(0, len(phis), chunk):
        L, HY, H2 = [], [], []
        for phi in phis[start:start + chunk]:
            ix, iy, iz, kz = cap_samples(cfg, grid, phi)
            lin = (ix * n + iy) * nz + iz
            h = G[ix, iy] * (cfg.k_m / kz)
            y = h * Of[lin] + sig * (rng.standard_normal(lin.size) + 1j * rng.standard_normal(lin.size))
            L.append(lin)
            HY.append(np.conj(h) * y)
            H2.append(np.abs(h) ** 2)
        L, HY, H2 = np.concatenate(L), np.concatenate(HY), np.concatenate(H2)
        num += np.bincount(L, weights=HY.real, minlength=grid.size) \
            + 1j * np.bincount(L, weights=HY.imag, minlength=grid.size)
        J += np.bincount(L, weights=H2, minlength=grid.size)
    return num.reshape(grid.shape) / sigma2, J.reshape(grid.shape) / sigma2




def run_design(cfg, grid, Ot, channels, sigma0, phis, rng):
    """channels: list of (dvec or None, time share). Returns per-channel
    (numS, JS) symmetrised, so every fusion rule can be formed."""
    out = []
    for dvec, t in channels:
        num, J = acquire(cfg, grid, Ot, dvec, phis, sigma0 ** 2 / t, rng)
        out.append((hsym_c(num), hermitian_sym(J)))
    return out




def main():
    rng = np.random.default_rng(SEED)
    cfg = GODT_SIM
    grid = KGrid(128, 0.1, 128)          # isotropic 0.1 um voxels, 12.8 um field
    grid.check(cfg)
    Kabs = grid.Kabs()
    phis = azimuths(N_ILL)
    t0 = time.time()
    M = hermitian_sym(cap_density(cfg, grid, N_ILL, N_FRAMES))
    supp = M > 0
    print(RULE)
    print("SETUP")
    print(RULE)
    print(f"  {cfg!r}\n  {grid!r}; {N_ILL} azimuths simulated, photons of {N_FRAMES} frames;"
          f" support {supp.sum()} voxels  [{time.time()-t0:.0f} s]")

    # ------------------------------------------------ phantom and prior
    O, dn = cell_phantom(grid, rng)
    Ot = np.fft.fftn(O)
    S, edges, prof = shell_spectrum(Ot, Kabs, cfg.kmax, grid.dk)
    ref = supp & (Kabs >= 8) & (Kabs <= 24)
    S_ref = np.mean(np.abs(Ot[ref]) ** 2)
    sigma0 = np.sqrt(S_ref / RHO)          # per-frame-sample noise at full budget
    meas = supp & (Kabs >= K_LO)           # where every aggregate below is taken
    O_band = np.real(np.fft.ifftn(Ot * meas))
    print(f"  textured cell phantom: dn in [{dn.min():+.4f}, {dn.max():+.4f}]; shell spectrum "
          f"S(|K|)/S(0): |K|=5: {prof[10]/prof[0]:.1e}, 10: {prof[20]/prof[0]:.1e}, "
          f"20: {prof[41]/prof[0]:.1e}, 30: {prof[61]/prof[0]:.1e}")
    print(f"  noise: rho = {RHO:.0f} relative to the mean power in 8 <= |K| <= 24; the per-voxel SNR "
          f"then varies with S(K) M(K) S_t(K_perp)")
    print(f"  median M on support {np.median(M[supp]):.1f}; per-voxel SNR at |K|=15, |G|^2 = 2: "
          f"{RHO*np.median(M[supp])*2*prof[31]/S_ref:.0f}")
    print()
    print("  TWO LOSSES THAT NO SHEAR DESIGN TOUCHES, quoted once (relative L2 of the truth):")
    print(f"    missing cone (|K| outside the sampled support): "
          f"{np.linalg.norm(np.real(np.fft.ifftn(Ot*~supp)))/np.linalg.norm(O):.3f}")
    print(f"    low-K shell |K| < {K_LO} rad/um, where |G|^2 -> 0 for every shear: "
          f"{np.linalg.norm(np.real(np.fft.ifftn(Ot*supp*(Kabs<K_LO))))/np.linalg.norm(O):.3f}")
    print(f"  Every aggregate below is therefore taken on supp & |K| >= {K_LO} "
          f"({meas.sum()} voxels, {100*meas.sum()/supp.sum():.1f} % of the support),")
    print("  and real-space RMSE is against O_band = the truth band-limited to that set.")

    # ------------------------------------------------ designs
    X, Y = (1.0, 0.0), (0.0, 1.0)

    def ch(mag, u, t):
        return ((mag * u[0], mag * u[1]), t)
    # (P-A) optimum for THIS grid / prior / SNR
    mags = np.round(np.arange(0.1, 1.6001, 0.1), 3)
    labels, G2 = [], []
    for dn, u in (("x", X), ("y", Y)):
        for m in mags:
            labels.append((float(m), dn))
            G2.append(G2_vec(grid.KX, grid.KY, (m * u[0], m * u[1])))
    G2 = np.array(G2)
    prob = DesignProblem(np.where(meas, S * M / sigma0 ** 2, 0.0), G2, np.ones(len(labels)))
    t0 = time.time()
    tA, fA, gapA = solve_design(prob, prob.f_A)
    opt = [ch(labels[i][0], X if labels[i][1] == "x" else Y, tA[i]) for i in np.nonzero(tA > 1e-3)[0]]
    print(f"  (P-A) optimum on this grid (empirical shell prior): f_A = {fA:.5f}, gap {gapA:.1e}, "
          f"{len(opt)} channels [{time.time()-t0:.0f} s]:")
    print("    " + ", ".join(f"{np.hypot(*d):.1f}{'x' if d[1]==0 else 'y'}:{t:.3f}" for d, t in opt))

    designs = [
        ("ODT benchmark", [(None, 1.0)]),
        ("0.1 x", [ch(0.1, X, 1.0)]),
        ("0.8 x", [ch(0.8, X, 1.0)]),
        ("0.8 x + 0.8 y", [ch(0.8, X, .5), ch(0.8, Y, .5)]),
        ("{0.8,0.55} x + y", [ch(0.8, X, .25), ch(0.55, X, .25), ch(0.8, Y, .25), ch(0.55, Y, .25)]),
        ("(P-A) optimum", opt),
    ]
    fusions = {"ODT benchmark": ["wiener"], "0.1 x": ["wiener"], "0.8 x": ["wiener", "blue"],
               "0.8 x + 0.8 y": ["wiener", "naive"],
               "{0.8,0.55} x + y": ["wiener", "blue", "naive"], "(P-A) optimum": ["wiener"]}

    nulls = nulls_in_band(0.8, cfg.kmax)
    KX3 = grid.KX[:, :, None]
    nullmask = np.zeros(supp.shape, bool)
    for nu in nulls:
        nullmask |= np.abs(np.abs(KX3) - nu) < 0.5
    nullmask &= meas
    idx_sh = np.clip(np.digitize(Kabs.ravel(), edges) - 1, 0, len(edges) - 2)
    mw = meas.ravel()
    pw = np.abs(Ot) ** 2

    def spectral(est, J):
        err2 = np.abs(est - Ot) ** 2
        pred = S / (S * J + 1.0)
        return (err2[meas].sum() / pw[meas].sum(),          # power-weighted, = RMSE^2
                pred[meas].sum() / pw[meas].sum(),
                np.mean((err2 / pw)[meas]),                 # voxel-mean relative, = f_A analogue
                np.mean((pred / S)[meas]),
                err2[nullmask].sum() / pw[nullmask].sum(),
                pred[nullmask].sum() / pw[nullmask].sum(), err2, pred)

    print()
    print(RULE)
    print("RESULTS at equal photon budget. RMSE is real-space, against O_band. specMSE is")
    print(f"power-weighted over supp & |K| >= {K_LO} (= RMSE^2); voxMSE is the unweighted voxel mean")
    print("of the relative error, which is the (P-A) objective. 'emp' = measured, 'pred' = theory.")
    print(RULE)
    print(f"  {'design':20s} {'fusion':7s} {'RMSE':>8} {'specMSE':>9} {'pred':>8} {'voxMSE':>8} {'pred':>8} "
          f"{'nullplane':>10} {'pred':>8}")
    results, profiles, recs = {}, {}, {}
    for name, chans in designs:
        t0 = time.time()
        acq = run_design(cfg, grid, Ot, chans, sigma0, phis, rng)
        for how in fusions[name]:
            est, J = fuse(acq, S, how)
            est = est * meas
            Oe = np.real(np.fft.ifftn(est))
            r_band = np.linalg.norm(Oe - O_band) / np.linalg.norm(O_band)
            me, mp, ve, vp, ne, npd, err2, pred = spectral(est, J)
            results[(name, how)] = (r_band, me, mp, ve, vp, ne, npd)
            recs[(name, how)] = Oe
            profiles[(name, how)] = (
                np.fft.fftshift((err2 * meas).sum(axis=(1, 2)) / np.maximum((pw * meas).sum(axis=(1, 2)), 1e-300)),
                np.fft.fftshift((pred * meas).sum(axis=(1, 2)) / np.maximum((pw * meas).sum(axis=(1, 2)), 1e-300)))
            if how == "wiener":
                e_sh = np.bincount(idx_sh[mw], weights=err2.ravel()[mw], minlength=len(edges) - 1)
                p_sh = np.bincount(idx_sh[mw], weights=pred.ravel()[mw], minlength=len(edges) - 1)
                s_sh = np.bincount(idx_sh[mw], weights=pw.ravel()[mw], minlength=len(edges) - 1)
                profiles[(name, "shell")] = (e_sh / np.maximum(s_sh, 1e-300), p_sh / np.maximum(s_sh, 1e-300))
            flag = "  <- diverges: J = 0 on the null planes" if r_band > 10 else ""
            print(f"  {name:20s} {how:7s} {r_band:8.4f} {me:9.4f} {mp:8.4f} {ve:8.4f} {vp:8.4f} "
                  f"{ne:10.3f} {npd:8.3f}   [{time.time()-t0:.0f} s]{flag}")
    print()
    print("  BLUE is the prior-free unbiased fusion. It is defined only where J > 0, so for a single")
    print("  shear it is undefined on the null planes and blows up beside them; with an incommensurate")
    print("  pair J > 0 everywhere on the band and BLUE and Wiener agree. That is the estimator-side")
    print("  statement of why multi-shear is not merely a noise-averaging improvement.")

    # ------------------------------------------------ Gaussian random field check
    print()
    print(RULE)
    print("GAUSSIAN RANDOM FIELD DRAWN FROM THE PRIOR: here the theory is exact, not approximate")
    print(RULE)
    white = np.fft.fftn(rng.standard_normal(O.shape))
    Og = np.real(np.fft.ifftn(np.sqrt(S) * white / np.sqrt(np.mean(np.abs(white) ** 2))))
    Ogt = np.fft.fftn(Og)
    print(f"  {'design':20s} {'total emp/pred':>14} {'shell median':>13} {'within 5%':>10} {'worst shell':>24}")
    for name, chans in designs[2:6]:
        acq = run_design(cfg, grid, Ogt, chans, sigma0, phis, rng)
        est, J = fuse(acq, S, "wiener")
        err2 = np.abs(est * meas - Ogt * meas) ** 2
        pred = S / (S * J + 1.0)
        e_sh = np.bincount(idx_sh[mw], weights=err2.ravel()[mw], minlength=len(edges) - 1)
        p_sh = np.bincount(idx_sh[mw], weights=pred.ravel()[mw], minlength=len(edges) - 1)
        c_sh = np.bincount(idx_sh[mw], minlength=len(edges) - 1)
        ok = c_sh > 2000
        ratio = e_sh[ok] / p_sh[ok]
        kk = (0.5 * (edges[:-1] + edges[1:]))[ok]
        iw = int(np.argmax(np.abs(ratio - 1)))
        # how concentrated is that shell's PREDICTED error? (a few near-null voxels
        # carrying most of the predicted error make the shell a small-sample estimate)
        shell = meas & (idx_sh.reshape(meas.shape) == np.nonzero(ok)[0][iw])
        pv = np.sort(pred[shell])[::-1]
        top = pv[:max(int(0.01 * pv.size), 1)].sum() / pv.sum()
        print(f"  {name:20s} {err2[meas].sum()/pred[meas].sum():14.4f} {np.median(ratio):13.4f} "
              f"{100*np.mean(np.abs(ratio-1) < 0.05):9.0f}% "
              f"{f'{ratio[iw]:.2f} at |K|={kk[iw]:.1f}':>24}  (top 1 % of its voxels hold "
              f"{100*top:.0f} % of the predicted error)")
    print("  => shell ratios are the meaningful test: each averages thousands of voxels. Where a")
    print("     shell's predicted error concentrates on a handful of near-null voxels, the empirical")
    print("     value is a small-sample estimate of an exponential variable and scatters accordingly;")
    print("     designs that leave no near-null voxels (incommensurate pair, P-A optimum) do not.")

    # ------------------------------------------------ figures
    kc = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
    for name, c in (("0.8 x", "C3"), ("0.8 x + 0.8 y", "C1"), ("{0.8,0.55} x + y", "C0"),
                    ("(P-A) optimum", "C2"), ("ODT benchmark", "k")):
        e, p = profiles[(name, "shell")]
        ax[0].semilogy(kc, p, c, lw=1.2)
        ax[0].semilogy(kc, e, c, marker="o", ms=3, ls="none", label=name)
    ax[0].set(xlabel="|K| (rad/um)", ylabel="relative MSE per shell", xlim=(K_LO, cfg.kmax),
              title="Wiener fusion: empirical (markers) vs predicted (lines)")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3)
    kxs = np.fft.fftshift(grid.kx)
    for name, c in (("0.8 x", "C3"), ("{0.8,0.55} x + y", "C0")):
        e, p = profiles[(name, "wiener")]
        ax[1].plot(kxs, p, c, lw=1.2)
        ax[1].plot(kxs, e, c, marker=".", ms=3, ls="none", label=name)
    for nu in nulls:
        ax[1].axvline(nu, color="k", ls=":", alpha=.4)
        ax[1].axvline(-nu, color="k", ls=":", alpha=.4)
    ax[1].set(xlabel="K_x (rad/um)", ylabel="relative MSE (summed over K_y, K_z)", ylim=(0, 1.05),
              title="null planes of 0.8 um (dotted) and their removal by fusion")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out_path("validation_spectra.png"), dpi=110)

    keys = [("truth", None), ("0.8 x", "wiener"), ("0.8 x + 0.8 y", "naive"),
            ("{0.8,0.55} x + y", "naive"), ("{0.8,0.55} x + y", "wiener"), ("(P-A) optimum", "wiener")]
    fig2, ax2 = plt.subplots(2, len(keys), figsize=(3.0 * len(keys), 6.2))
    cz, cy = grid.nz // 2, grid.n // 2
    vmin, vmax = O_band.min(), O_band.max()
    for j, (name, how) in enumerate(keys):
        vol = O_band if name == "truth" else recs[(name, how)]
        ax2[0, j].imshow(vol[:, :, cz].T, cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
        ax2[1, j].imshow(vol[:, cy, :].T, cmap="gray", vmin=vmin, vmax=vmax, origin="lower", aspect=grid.dz / grid.dx)
        ttl = "band-limited truth" if name == "truth" else f"{name}\n{how}"
        ax2[0, j].set_title(ttl, fontsize=8)
        for a in ax2[:, j]:
            a.axis("off")
    ax2[0, 0].text(-8, grid.n // 2, "x-y", rotation=90, va="center", fontsize=8)
    ax2[1, 0].text(-8, grid.nz // 2, "x-z", rotation=90, va="center", fontsize=8)
    fig2.suptitle(f"Reconstructions at equal photon budget (rho = {RHO:.0f}); band-passed to the sampled support",
                  fontsize=10)
    fig2.tight_layout()
    fig2.savefig(out_path("validation_slices.png"), dpi=110)
    print()
    print(f"  wrote {out_path('validation_spectra.png')}")
    print(f"  wrote {out_path('validation_slices.png')}")


if __name__ == "__main__":
    main()
