"""exp04 -- Optimal shear design: (P-A), (P-D), (P-W) solutions, C1-C3, E1-E2.

Establishes (design note T3, T5):
  1. Solutions of the three convex design problems over a candidate grid of
     shear magnitudes x directions at equal photon budget, with the
     equivalence-theorem gap as an optimality certificate.
  2. C1 (sparsity): how many channels the optimum uses, and whether extra
     directions (45/135 deg) get any weight.
  3. C2 (coherence-limited optimum): how the optimal magnitude moves when
     fringe visibility falls with shear (jinc model, coherence length lc).
  4. C3 (incommensurability emerges): the objective vs a second magnitude
     paired with 0.8 um has its worst values at small-integer ratios.
  5. E1 (calibration): an rms shear error s caps the per-voxel SNR at
     |G|^2/(s K_perp)^2 and makes the optimal design spread over more channels.
  6. E2 (chromatic shear): under GODT's 35 nm LED the nulls survive to first
     order (spectral averaging fills them only at O(bandwidth^2)).

Writes: out/design_shares.png, out/design_scan.png, out/design_maps.png
Wiki:   wiki/notes/multishear-fusion-design-solution.md

Run: python exp04_shear_design.py     (~5-8 min; reuses out/M_*.npy from exp03)
"""
import time
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, KGrid, G2_vec, visibility, cached_density, hermitian_sym,
                      cap_density, DesignProblem, solve_design, solve_maximin, chromatic_G,
                      nulls_in_band, out_path)

RHO = 100.0
N_AZ, N_FRAMES = 1440, 180
MAGS = np.round(np.arange(0.05, 1.6001, 0.05), 3)          # 32 magnitudes, um
DIRS2 = {"x": (1.0, 0.0), "y": (0.0, 1.0)}
DIRS4 = {"x": (1.0, 0.0), "d": (np.sqrt(.5), np.sqrt(.5)), "y": (0.0, 1.0),
         "e": (-np.sqrt(.5), np.sqrt(.5))}
RULE = "=" * 78


def candidates(grid, mags, dirs):
    labels, G2, P = [], [], []
    for dn, (ux, uy) in dirs.items():
        for m in mags:
            labels.append((float(m), dn))
            G2.append(G2_vec(grid.KX, grid.KY, (m * ux, m * uy)))
            P.append((grid.KX * ux + grid.KY * uy) ** 2)
    return labels, np.array(G2), np.array(P)


def show(t, labels, thr=1e-4):
    idx = np.nonzero(t > thr)[0]
    return f"{len(idx):2d} ch: " + ", ".join(f"{labels[i][0]:.2f}{labels[i][1]}:{t[i]:.3f}" for i in idx)


def fixed(labels, spec):
    """Time-share vector from {(magnitude, direction): weight}, matched on
    magnitude to 1e-6 so float round-off cannot silently miss a candidate."""
    t = np.zeros(len(labels))
    for (m, dn), v in spec.items():
        hits = [i for i, l in enumerate(labels) if l[1] == dn and abs(l[0] - m) < 1e-6]
        assert len(hits) == 1, f"candidate {(m, dn)} not in the grid"
        t[hits[0]] = v
    return t / t.sum()


def near_fraction(r, maxden=6):
    f = Fraction(r).limit_denominator(maxden)
    return f, abs(float(f) - r) / r


def main():
    cfg = GODT_SIM
    grid = KGrid(256, 0.1, 96)
    grid.check(cfg)
    M = cached_density(cfg, grid, N_AZ, N_FRAMES)
    supp = M > 0
    A = RHO * M
    labels, G2, P = candidates(grid, MAGS, DIRS2)
    v1 = np.ones(len(labels))
    prob = DesignProblem(A, G2, v1)
    band2d = (grid.KP <= cfg.kmax) & (grid.KP > 0) & supp.any(axis=2)

    def eff_worst(t):
        E = G2
        Emax = E.max(axis=0)
        St = np.tensordot(t, E, axes=(0, 0))
        return (St / np.where(Emax > 0, Emax, 1))[band2d & (Emax > 0)].min()

    # ------------------------------------------------------- 1. baselines
    print(RULE)
    print(f"1. (P-A) OBJECTIVE = mean relative Wiener MSE over the support, rho = {RHO:.0f}, flat prior")
    print(f"   {len(labels)} candidates: {len(MAGS)} magnitudes {MAGS[0]}..{MAGS[-1]} um x directions x, y")
    print(RULE)
    base = {
        "0.1 x": {(0.1, "x"): 1},
        "0.8 x": {(0.8, "x"): 1},
        "0.8 x + 0.8 y": {(0.8, "x"): 1, (0.8, "y"): 1},
        "{0.8,0.55} x + y": {(0.8, "x"): 1, (0.55, "x"): 1, (0.8, "y"): 1, (0.55, "y"): 1},
        "{0.8,0.55,0.35} x + y": {(0.8, "x"): 1, (0.55, "x"): 1, (0.35, "x"): 1,
                                  (0.8, "y"): 1, (0.55, "y"): 1, (0.35, "y"): 1},
        "uniform over all 64": {l: 1 for l in labels},
    }
    print(f"  {'design':28s} {'f_A':>9} {'worst eff':>10}")
    f_odt = (1.0 / (A + 1.0))[supp].mean()
    print(f"  {'ODT benchmark |G|^2=1':28s} {f_odt:9.5f} {'-':>10}")
    print("  READ THE BENCHMARK WITH CARE: |G|^2 = 4 sin^2(K.d/2) peaks at 4, so on this model a")
    print("  shear design can SCORE BETTER than direct ODT. That is assumption A4 talking -- the")
    print("  phase difference is one measurement whose noise does not grow with |d| -- together")
    print("  with unit fringe visibility. Section 4 prices the visibility; the benchmark is here")
    print("  as a fixed reference level, not as a claim that gradient beats direct ODT.")
    for lbl, spec in base.items():
        t = fixed(labels, spec)
        print(f"  {lbl:28s} {prob.f_A(t)[0]:9.5f} {eff_worst(t):10.3f}")

    # the natural shear scale: the largest shear with no in-band null is
    # d = 2 pi / kmax, where the first null sits exactly at the band edge
    d_nat = 2 * np.pi / cfg.kmax
    single = np.array([prob.f_A(fixed(labels, {(float(m), "x"): 1}))[0] for m in MAGS])
    pairs = np.array([prob.f_A(fixed(labels, {(float(m), "x"): 1, (float(m), "y"): 1}))[0] for m in MAGS])
    print()
    print(f"  scan over one magnitude (kmax = {cfg.kmax:.2f} rad/um, so the first null leaves the band")
    print(f"  at d = 2 pi / kmax = {d_nat:.3f} um):")
    print(f"    best single channel: d = {MAGS[single.argmin()]:.2f} um, f_A = {single.min():.5f}"
          f"   (GODT sim 0.1: {single[np.argmin(np.abs(MAGS-0.1))]:.5f}; "
          f"published exp 0.8: {single[np.argmin(np.abs(MAGS-0.8))]:.5f})")
    print(f"    best x+y pair:       d = {MAGS[pairs.argmin()]:.2f} um, f_A = {pairs.min():.5f}"
          f"   (0.1: {pairs[np.argmin(np.abs(MAGS-0.1))]:.5f}; 0.8: {pairs[np.argmin(np.abs(MAGS-0.8))]:.5f})")

    # ------------------------------------------------------- 2. optima
    print()
    print(RULE)
    print("2. OPTIMAL DESIGNS")
    print(RULE)
    t0 = time.time()
    tA, fA, gapA = solve_design(prob, prob.f_A)
    print(f"  (P-A)  f_A = {fA:.5f}  gap = {gapA:.1e}  [{time.time()-t0:.0f} s]")
    print(f"         {show(tA, labels)}")
    print(f"         worst-case efficiency {eff_worst(tA):.3f}; vs ODT benchmark {fA/f_odt:.2f}x")
    t0 = time.time()
    tD, fD, gapD = solve_design(prob, prob.f_D)
    print(f"  (P-D)  mean log(1+SJ) = {-fD:.4f}  gap = {gapD:.1e}  [{time.time()-t0:.0f} s]")
    print(f"         {show(tD, labels)};  f_A = {prob.f_A(tD)[0]:.5f}")
    t0 = time.time()
    tW, zlp, zfull = solve_maximin(G2, v1, band2d, subsample=2)
    print(f"  (P-W)  maximin efficiency: LP {zlp:.4f}, on full grid {zfull:.4f}  [{time.time()-t0:.0f} s]")
    print(f"         {show(tW, labels, 1e-3)};  f_A = {prob.f_A(tW)[0]:.5f}")
    print("  (efficiency = S_t / max_j |G_j|^2 at that K_perp; 1 = as good as the best")
    print("   single candidate everywhere, which no fixed design reaches)")
    # SNR sensitivity
    for rho in (10.0, 1000.0):
        p2 = DesignProblem(rho * M, G2, v1)
        t2, f2, g2_ = solve_design(p2, p2.f_A)
        print(f"  (P-A) at rho = {rho:5.0f}: f_A = {f2:.5f}  {show(t2, labels)}")

    # ------------------------------------------------------- 3. C1 directions
    print()
    print(RULE)
    print("3. C1 -- SPARSITY, AND DO EXTRA DIRECTIONS HELP?")
    print(RULE)
    print(f"  support sizes: P-A {int((tA>1e-4).sum())}, P-D {int((tD>1e-4).sum())}, "
          f"P-W {int((tW>1e-3).sum())} of {len(labels)} candidates")
    labels4, G24, _ = candidates(grid, MAGS, DIRS4)
    prob4 = DesignProblem(A, G24, np.ones(len(labels4)))
    t0 = time.time()
    tA4, fA4, gapA4 = solve_design(prob4, prob4.f_A)
    wdiag = sum(tA4[i] for i, l in enumerate(labels4) if l[1] in ("d", "e"))
    print(f"  with 4 directions (x, 45deg, y, 135deg; {len(labels4)} candidates): f_A = {fA4:.5f} "
          f"(2-dir optimum {fA:.5f}), gap {gapA4:.1e}, weight on diagonals {wdiag:.3f}  [{time.time()-t0:.0f} s]")
    print(f"         {show(tA4, labels4)}")

    # ------------------------------------------------------- 4. C2 coherence
    print()
    print(RULE)
    print("4. C2 -- COHERENCE-LIMITED OPTIMUM  (visibility^2 = jinc^2(pi d / lc), first zero 1.22 lc)")
    print(RULE)
    print(f"  {'lc (um)':>8} {'best pair d*':>13} {'f_A(pair)':>10}   P-A optimum")
    pair_curves = {}
    for lc in (None, 4.0, 2.0, 1.0):
        v = visibility(np.array([l[0] for l in labels]), lc)
        pl = DesignProblem(A, G2, v)
        vals = []
        for m in MAGS:
            t = fixed(labels, {(float(m), "x"): 1, (float(m), "y"): 1})
            vals.append(pl.f_A(t)[0])
        vals = np.array(vals)
        pair_curves[lc] = vals
        i = int(np.argmin(vals))
        tl, fl, gl = solve_design(pl, pl.f_A)
        print(f"  {('inf' if lc is None else lc):>8} {MAGS[i]:13.2f} {vals[i]:10.5f}   "
              f"f_A = {fl:.5f}, gap {gl:.1e}: {show(tl, labels)}")

    # ------------------------------------------------------- 5. C3 ratios
    print()
    print(RULE)
    print("5. C3 -- INCOMMENSURABILITY: 0.8 um pair (x+y) plus a second pair at d2, equal shares")
    print(RULE)
    d2s = np.round(np.arange(0.05, 1.6001, 0.01), 3)
    fA_scan, worst_scan = [], []
    G08x, G08y = G2_vec(grid.KX, grid.KY, (0.8, 0)), G2_vec(grid.KX, grid.KY, (0, 0.8))
    band_k2 = band2d & (grid.KP > 2.0)
    for d2 in d2s:
        St = 0.25 * (G08x + G08y + G2_vec(grid.KX, grid.KY, (d2, 0)) + G2_vec(grid.KX, grid.KY, (0, d2)))
        rel = 1.0 / (A * St[:, :, None] + 1.0)
        fA_scan.append(rel[supp].mean())
        worst_scan.append(St[band_k2].min())
    fA_scan, worst_scan = np.array(fA_scan), np.array(worst_scan)
    # local maxima of f_A (worst second shears)
    peaks = [i for i in range(1, len(d2s) - 1) if fA_scan[i] >= fA_scan[i-1] and fA_scan[i] >= fA_scan[i+1]]
    peaks = sorted(peaks, key=lambda i: -fA_scan[i])[:8]
    print(f"  {'d2 (um)':>8} {'ratio 0.8/d2':>13} {'nearest p/q':>12} {'f_A':>9} {'min S_t (K>2)':>14}")
    for i in sorted(peaks):
        r = 0.8 / d2s[i]
        f, err = near_fraction(r)
        print(f"  {d2s[i]:8.2f} {r:13.3f} {str(f):>12} {fA_scan[i]:9.5f} {worst_scan[i]:14.4f}")
    ib = int(np.argmin(fA_scan))
    print(f"  best second pair: d2 = {d2s[ib]:.2f} um (ratio {0.8/d2s[ib]:.3f}), f_A = {fA_scan[ib]:.5f}; "
          f"single 0.8 pair: {fA_scan[np.argmin(np.abs(d2s-0.8))]:.5f}")
    print()
    print("  Does the rule still bind once the design has many channels? min S_t over the band")
    print("  (K_perp > 2), which is what a coincident null would drive to zero:")
    for lbl, t in (("0.8 pair alone", fixed(labels, {(0.8, "x"): 1, (0.8, "y"): 1})),
                   ("0.8 + 0.4 pairs (ratio 2)", fixed(labels, {(0.8, "x"): 1, (0.4, "x"): 1,
                                                                (0.8, "y"): 1, (0.4, "y"): 1})),
                   ("0.8 + 0.55 pairs", fixed(labels, {(0.8, "x"): 1, (0.55, "x"): 1,
                                                       (0.8, "y"): 1, (0.55, "y"): 1})),
                   ("(P-A) optimum", tA), ("(P-W) optimum", tW)):
        St = np.tensordot(t, G2, axes=(0, 0))
        xs = sorted({l[0] for i, l in enumerate(labels) if t[i] > 1e-4})
        # do any two active magnitudes share an in-band null?
        shared = sorted({round(a / b, 3) for a in xs for b in xs if b < a
                         and nulls_in_band(b, cfg.kmax) and abs(a / b - round(a / b)) < 0.02})
        print(f"    {lbl:26s} min S_t = {St[band_k2].min():.4f}   {len(xs)} magnitudes"
              f"{'; integer ratios present: ' + str(shared) if shared else '; no integer ratios'}")
    print("  => the incommensurate-ratio rule is a rule for SMALL channel sets. The optimum")
    print("     contains integer ratios and is still fine, because a null shared by two of its")
    print("     twelve magnitudes is covered by the other ten.")

    # ------------------------------------------------------- 6. E1 calibration
    print()
    print(RULE)
    print("6. E1 -- SHEAR CALIBRATION ERROR s (rms, um; 0.01 um = 0.1 px at dx = 0.1)")
    print(RULE)
    g1 = KGrid(128, 0.1, 48)
    M1 = hermitian_sym(cap_density(cfg, g1, N_AZ, N_FRAMES))
    mags1 = np.round(np.arange(0.1, 1.6001, 0.1), 3)
    labels1, G21, P1 = candidates(g1, mags1, DIRS2)
    p1 = DesignProblem(RHO * M1, G21, np.ones(len(labels1)), P=P1)
    print("  calibration-limited SNR ceiling |G|^2/(s K_perp)^2 at |G|^2 = 4:")
    for s in (0.005, 0.01, 0.02):
        print(f"    s = {s:.3f} um: K = 10 -> {4/(s*10)**2:8.0f} ({10*np.log10(4/(s*10)**2):4.1f} dB); "
              f"K = 20 -> {4/(s*20)**2:8.0f} ({10*np.log10(4/(s*20)**2):4.1f} dB); "
              f"K = 30 -> {4/(s*30)**2:8.0f} ({10*np.log10(4/(s*30)**2):4.1f} dB)")
    print(f"  photon-limited per-voxel SNR at rho = {RHO:.0f}, median M = {np.median(M1[M1>0]):.1f}, |G|^2 = 4: "
          f"{RHO*np.median(M1[M1>0])*4:.0f} ({10*np.log10(RHO*np.median(M1[M1>0])*4):.1f} dB)")
    print()
    print(f"  {'s (um)':>7} {'f_A':>9} {'#ch':>4}   design")
    tprev = None
    for s in (0.0, 0.005, 0.01, 0.02, 0.05):
        f = p1.f_A if s == 0 else (lambda t, s=s: p1.f_A_robust(t, s))
        t0 = time.time()
        ts, fs, gs = solve_design(p1, f, t0=tprev, gap_tol=1e-3, rounds=6)
        tprev = np.maximum(ts, 1e-3)
        print(f"  {s:7.3f} {fs:9.5f} {int((ts>1e-3).sum()):4d}   {show(ts, labels1, 1e-3)}   [{time.time()-t0:.0f} s]")

    # ------------------------------------------------------- 7. E2 chromatic
    print()
    print(RULE)
    print("7. E2 -- CHROMATIC SHEAR: d(lam) = d0 lam/lam0, GODT LED 530 nm, FWHM 35 nm")
    print(RULE)
    lam0, d0 = 0.53, 0.8
    nulls = nulls_in_band(d0, cfg.kmax)
    Kq = np.array(nulls + [np.pi / d0])
    print(f"  {'FWHM':>6} {'K (rad/um)':>11} {'|G_eff|^2/4':>12} {'2nd-order pred':>15} {'incoh. bound/4':>15} {'relMSE @ rho M=650':>19}")
    for fwhm in (0.0, 0.035, 0.070):
        Ge, G2i = chromatic_G(Kq, d0, lam0, max(fwhm, 1e-9))
        sig_e = fwhm / (2 * np.sqrt(2 * np.log(2))) / lam0
        for m, K in enumerate(Kq, start=1):
            pred = (2 * np.pi ** 2 * m ** 2 * sig_e ** 2) ** 2 / 4 if m <= len(nulls) else np.nan
            rel = 1.0 / (650 * np.abs(Ge[m-1]) ** 2 + 1.0)
            tag = f"null m={m}" if m <= len(nulls) else "peak"
            print(f"  {1e3*fwhm:4.0f}nm {K:8.2f} {tag:>9} {np.abs(Ge[m-1])**2/4:12.2e} "
                  f"{pred:15.2e} {G2i[m-1]/4:15.3f} {rel:19.3f}")
    print("  => the coherent-average filter keeps the nulls to first order in the bandwidth;")
    print("     the O(bandwidth^2) residue leaves them lost at any practical SNR. The incoherent")
    print("     bound (each wavelength detected separately) is not what a broadband detector sees.")
    # design under broadband G
    G2c = np.empty_like(G2)
    for j, (m, dn) in enumerate(labels):
        ux, uy = DIRS2[dn]
        Ge, _ = chromatic_G((grid.KX * ux + grid.KY * uy).ravel(), m, lam0, 0.035, n_lam=31)
        G2c[j] = (np.abs(Ge) ** 2).reshape(grid.KX.shape)
    pc = DesignProblem(A, G2c, v1)
    tc, fc, gc = solve_design(pc, pc.f_A)
    print(f"  (P-A) with the 35 nm LED filters: f_A = {fc:.5f} (mono {fA:.5f}), gap {gc:.1e}")
    print(f"         {show(tc, labels)}")
    print(f"         single 0.8 x under the LED: f_A = {pc.f_A(fixed(labels, {(0.8,'x'):1}))[0]:.5f} "
          f"(mono {prob.f_A(fixed(labels, {(0.8,'x'):1}))[0]:.5f})")

    # ------------------------------------------------------- figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), sharey=True)
    for a, (t, ttl) in zip(ax, [(tA, "(P-A) mean MSE"), (tD, "(P-D) capacity"), (tW, "(P-W) maximin efficiency")]):
        for dn, c, off in (("x", "C0", -0.01), ("y", "C3", 0.01)):
            ms = np.array([l[0] for l in labels if l[1] == dn])
            ts = np.array([t[i] for i, l in enumerate(labels) if l[1] == dn])
            a.bar(ms + off, ts, width=0.02, color=c, label=dn)
        a.set(title=ttl, xlabel="shear magnitude (um)")
        a.grid(alpha=.3)
    ax[0].set_ylabel("time share")
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(out_path("design_shares.png"), dpi=110)

    fig2, ax2 = plt.subplots(1, 2, figsize=(13, 4.2))
    ax2[0].plot(d2s, fA_scan, "C0")
    for q in (Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(1, 3), Fraction(2, 3),
              Fraction(3, 2), Fraction(3, 4), Fraction(4, 3)):
        d2 = 0.8 * float(q)
        if d2s[0] <= d2 <= d2s[-1]:
            ax2[0].axvline(d2, color="k", ls=":", alpha=.4)
            ax2[0].text(d2, ax2[0].get_ylim()[1] * 0.98, f"{q.denominator}:{q.numerator}",
                        rotation=90, va="top", ha="right", fontsize=7)
    ax2[0].set(xlabel="second magnitude d2 (um), paired with 0.8 um (x and y, equal shares)",
               ylabel="f_A", title="C3: worst second shears sit at small-integer ratios d2/0.8")
    ax2[0].grid(alpha=.3)
    for lc, vals in pair_curves.items():
        ax2[1].plot(MAGS, vals, label=f"lc = {'inf' if lc is None else lc} um")
    ax2[1].set(xlabel="pair magnitude d (um), x + y equal shares", ylabel="f_A", yscale="log",
               title="C2: orthogonal pair objective vs magnitude and coherence")
    ax2[1].legend(fontsize=8)
    ax2[1].grid(alpha=.3)
    fig2.tight_layout()
    fig2.savefig(out_path("design_scan.png"), dpi=110)

    fig3, ax3 = plt.subplots(1, 2, figsize=(11, 4.6))
    ext = [grid.kx.min(), -grid.kx.min(), grid.ky.min(), -grid.ky.min()]
    for a, (t, ttl) in zip(ax3, [(fixed(labels, base["0.8 x + 0.8 y"]), "S_t: 0.8 x + 0.8 y"),
                                 (tA, "S_t: (P-A) optimum")]):
        St = np.fft.fftshift(prob.S_t(t)).T
        im = a.imshow(St, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=4)
        a.set(title=ttl, xlabel="K_x (rad/um)", ylabel="K_y (rad/um)")
        plt.colorbar(im, ax=a, fraction=0.046)
    fig3.tight_layout()
    fig3.savefig(out_path("design_maps.png"), dpi=110)
    print()
    for f in ("design_shares.png", "design_scan.png", "design_maps.png"):
        print(f"  wrote {out_path(f)}")


if __name__ == "__main__":
    main()
