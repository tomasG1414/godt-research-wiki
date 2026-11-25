"""exp11 -- exp06's design comparison, with both implementation errors removed.

exp06 measured how the multi-shear design margins shrink from the ideal linear
model to realistic data, but two of its steps were implementation choices,
not physics:
  - nearest-voxel K_z deposition (assumption A3), which exp10 shows is ~73 %
    (in log terms) of its A -> B factor of 6.2; and
  - model D = B + (WPM excess, measured on a cropped 19.2 um window), whose
    excess exp09 shows is dominated by that window and by WPM's own errors.
This script re-runs exp06's stage 4 with
  C  MSB: multi-slice first Born (shearlib.MSB) on a 38.4 um window, measured
     on the whole window at the reconstruction lattice (shearlib.lattice_samples).
     Its linear limit is first Born by construction and, on this window, its
     measurement converges on periodic Born's to ~2 % (exp09). No splice.
  D  B + (C - (c/eps) C(eps)): the same excess on exactly Born's linear part,
     the analogue of exp06's D, for comparison.
and reconstructs every model with both estimators: exp06's nearest-voxel rule
and exp10's gridding-free column solver (shearlib.ColumnWiener).

What C is and is not: on the near-axis rays, where the exact excess is
window-converged, MSB reproduces CBS's excess to ~10 % (exp09 stage 3). On
steep rays the exact referee is not window-converged either, so C there is
MSB's answer, not a verified one. Everything else is exp06's stage 4: phantom
(SEED), 180 illuminations snapped to the 19.2 um lattice (a sublattice of
38.4 um), pupil, prior, rho, per-ray noise and equal-photon shares; the A and B
columns reproduce exp06's, which the script asserts.

Deterministic. ~45 min on 4 cores, ~6 GB.
Wiki: wiki/notes/gridding-and-forward-model.md
Run:  python exp11_exact_design.py          (contrasts 0.5, 1.0, 2.0)
      python exp11_exact_design.py --quick  (contrast 0.5 only, ~15 min)
"""
import sys
import time
import numpy as np

from shearlib import (GODT_SIM, KGrid, hermitian_sym, snap_illumination, cap_geometry,
                      born_rytov_ratio, meas_scale, cell_phantom, shell_spectrum, hsym_c,
                      fuse, MSB, ColumnWiener, lattice_samples, shear_rytov, out_path)
from exp06_multiple_scattering import (SEED, RHO, N_ILL, K_LO, PHANTOM_SCALE, PAD, EPS,
                                       shear_factor, channel_meas, bin_add)

QUICK = "--quick" in sys.argv
CONTRASTS = (0.5,) if QUICK else (0.5, 1.0, 2.0)
BIG = 384               # 38.4 um: holds both the 19.2 um and the 12.8 um lattices
FRAME = 0.1             # absorbing edge outside |x| = 15.4 um
Z0, NZC = 28, 72        # slices 28..99 of exp06's 128 hold the whole cell
RULE = "=" * 100
MODELS = ("A", "B", "C", "D")
# exp06 stage 4 (wiki/notes/multiple-scattering-check.md section 5)
EXP06 = {0.5: {("ODT benchmark", "wiener"): (0.0221, 0.1067, 0.2883),
               ("0.1 x", "wiener"): (0.2172, 0.2395, 0.3541),
               ("0.8 x", "wiener"): (0.2453, 0.2651, 0.3704),
               ("0.8 x + 0.8 y", "wiener"): (0.0408, 0.1120, 0.2898),
               ("0.8 x + 0.8 y", "naive"): (0.1804, 0.2090, 0.3287),
               ("{0.8,0.55} x + y", "wiener"): (0.0171, 0.1058, 0.2882),
               ("{0.8,0.55} x + y", "naive"): (0.1642, 0.1952, 0.3197)},
         1.0: {("{0.8,0.55} x + y", "wiener"): (0.0172, 0.1057, 0.6109),
               ("{0.8,0.55} x + y", "naive"): (None, None, 0.6025)},
         2.0: {("{0.8,0.55} x + y", "wiener"): (0.0172, 0.1057, 1.6349)}}


def main():
    t_all = time.time()
    cfg = GODT_SIM
    grid = KGrid(128, 0.1, 128)
    Kabs = grid.Kabs()
    k0, nmed = 2 * np.pi / cfg.lam, cfg.n0
    rng = np.random.default_rng(SEED)
    O1, dn1 = cell_phantom(grid, rng, scale=PHANTOM_SCALE)
    assert not np.any(dn1[:, :, :Z0]) and not np.any(dn1[:, :, Z0 + NZC:]), "z crop"
    KI, th = snap_illumination(cfg, grid, N_ILL, dk=2 * np.pi / (PAD * grid.dx))
    caps = [cap_geometry(cfg, grid, ki) for ki in KI]
    sc = meas_scale(cfg, grid)
    cover = np.zeros(grid.size, np.int32)
    for cap in caps:
        ixc, iyc = np.nonzero(cap["ok"])
        cover += np.bincount((ixc * grid.n + iyc) * grid.nz + cap["iz"][ixc, iyc],
                             minlength=grid.size).astype(np.int32)
    supp = hermitian_sym(cover.reshape(grid.shape)) > 0
    ref = supp & (Kabs >= 8) & (Kabs <= 24)
    meas = supp & (Kabs >= K_LO)
    CH = {"odt": None, "0.1x": (0.1, 0.0), "0.8x": (0.8, 0.0),
          "0.8y": (0.0, 0.8), "0.55x": (0.55, 0.0), "0.55y": (0.0, 0.55)}
    designs = [("ODT benchmark", ["odt"], ["wiener"]),
               ("0.1 x", ["0.1x"], ["wiener"]),
               ("0.8 x", ["0.8x"], ["wiener"]),
               ("0.8 x + 0.8 y", ["0.8x", "0.8y"], ["wiener", "naive"]),
               ("{0.8,0.55} x + y", ["0.8x", "0.55x", "0.8y", "0.55y"], ["wiener", "naive"])]

    def big_O(c):
        O = np.zeros((BIG, BIG, NZC))
        o = (BIG - grid.n) // 2
        O[o:o + grid.n, o:o + grid.n] = k0 ** 2 * (nmed ** 2 - (nmed + c * dn1[:, :, Z0:Z0 + NZC]) ** 2)
        return O

    def msb_meas(model, ki, ixj, iyj, phase):
        """per-channel pupil data of one MSB run, measured on the whole window"""
        r = model.ratio(ki)
        out = {}
        for k, dvec in CH.items():
            f = np.log(r) if dvec is None else shear_rytov(r, grid.dx, dvec, wrap_safe=True)
            out[k] = sc * lattice_samples(f, grid.n, ixj, iyj) * phase
        return out

    print(RULE)
    print("exp06's DESIGN COMPARISON WITHOUT ITS TWO IMPLEMENTATION ERRORS")
    print(RULE)
    print(f"  {grid!r}; {N_ILL} illuminations at {th.mean():.2f} +- {th.std():.2f} deg; rho {RHO:.0f}")
    print(f"  C: MSB on a {BIG*grid.dx:.1f} um window (absorbing edge outside |x| = "
          f"{0.5*BIG*grid.dx*(1-2*FRAME):.1f} um), measured on the whole window")
    print("  D: B + (C - (c/eps) C(eps)).  Columns marked 'x' use the gridding-free estimator.")

    t0 = time.time()
    msb_e = MSB(big_O(EPS), grid.dx, grid.dz, cfg.k_m, absorb_frac=FRAME)
    yE = []                                   # the eps run, per illumination and channel
    for j, (ki, cap) in enumerate(zip(KI, caps)):
        ixj, iyj = np.nonzero(cap["ok"])
        phase = np.exp(-1j * cap["Kz_vox"][ixj, iyj] * cap["z_c"])
        yE.append(msb_meas(msb_e, ki, ixj, iyj, phase))
    del msb_e
    print(f"  MSB at eps = {EPS}: {N_ILL} illuminations [{time.time()-t0:.0f} s]", flush=True)

    out = {}
    for c in CONTRASTS:
        Oc = k0 ** 2 * (nmed ** 2 - (nmed + c * dn1) ** 2)
        Otc = np.fft.fftn(Oc)
        Sc, _, _ = shell_spectrum(Otc, Kabs, cfg.kmax, grid.dk)
        sig0 = np.sqrt(np.mean(np.abs(Otc[ref]) ** 2) / RHO)
        Ob = np.real(np.fft.ifftn(Otc * meas))
        msb_c = MSB(big_O(c), grid.dx, grid.dz, cfg.k_m, absorb_frac=FRAME)
        numS = {(m, k): np.zeros(grid.size, np.complex128) for m in MODELS for k in CH}
        numN = {k: np.zeros(grid.size, np.complex128) for k in CH}
        J1 = {k: np.zeros(grid.size) for k in CH}
        rays = {"cols": [], "Kz": [], "c": []}
        pr = {(m, k): [] for m in ("B", "C", "D", "N", "w") for k in CH}
        t0 = time.time()
        for j, (ki, cap) in enumerate(zip(KI, caps)):
            ixj, iyj = np.nonzero(cap["ok"])
            lin = (ixj * grid.n + iyj) * grid.nz + cap["iz"][ixj, iyj]
            phase = np.exp(-1j * cap["Kz_vox"][ixj, iyj] * cap["z_c"])
            hbase = cfg.k_m / cap["kz"][ixj, iyj]
            rays["cols"].append(ixj * grid.n + iyj)
            rays["Kz"].append(cap["Kz"][ixj, iyj])
            rays["c"].append(phase)
            rb = born_rytov_ratio(Oc, cfg, grid, ki, cap)
            lgB = np.log(rb)
            yC = msb_meas(msb_c, ki, ixj, iyj, phase)
            gen = np.random.default_rng((SEED, j))           # exp06's noise draws
            for k, dvec in CH.items():
                h = shear_factor(grid, dvec)[ixj, iyj] * hbase
                w = np.abs(h) ** 2 / sig0 ** 2
                J1[k] += np.bincount(lin, weights=w, minlength=grid.size)
                nz_ = gen.standard_normal((2, lin.size))
                noise = np.sqrt(sig0 ** 2 / 2) * (nz_[0] + 1j * nz_[1])
                gN = np.conj(h) * noise / sig0 ** 2
                bin_add(numN[k], lin, gN)
                y = {"A": h * Otc.reshape(-1)[lin],
                     "B": (sc * channel_meas(lgB, rb, grid, dvec))[ixj, iyj] * phase,
                     "C": yC[k]}
                y["D"] = y["B"] + y["C"] - (c / EPS) * yE[j][k]
                for m in MODELS:
                    g = np.conj(h) * y[m] / sig0 ** 2
                    bin_add(numS[(m, k)], lin, g)
                    if m != "A":
                        pr[(m, k)].append(g)
                pr[("N", k)].append(gN)
                pr[("w", k)].append(w)
            if (j + 1) % 60 == 0:
                print(f"    contrast {c}: illumination {j+1}/{N_ILL} [{time.time()-t0:.0f} s]",
                      flush=True)
        del msb_c
        rays = {key: np.concatenate(v) for key, v in rays.items()}
        cw = ColumnWiener(grid, rays["cols"], rays["Kz"], rays["c"])
        T = {k: cw.data_term(np.concatenate(pr[("w", k)])) for k in CH}
        bX = {(m, k): cw.rhs(np.concatenate(pr[(m, k)])) for m in ("B", "C", "D", "N") for k in CH}
        del pr

        def chan(m, k, t):
            num = (t * numS[(m, k)] + np.sqrt(t) * numN[k]).reshape(grid.shape)
            return hsym_c(num), hermitian_sym((t * J1[k]).reshape(grid.shape))

        def exact(ks, t, how):
            """gridding-free estimates for B, C, D at once (one factorisation)"""
            rhs = lambda m, kk: [t * bX[(m, k)] + np.sqrt(t) * bX[("N", k)] for k in kk]
            if how == "wiener":
                return cw.solve([t * T[k] for k in ks], [rhs(m, ks) for m in "BCD"], Sc)
            per = [cw.solve([t * T[k]], [rhs(m, [k]) for m in "BCD"], Sc) for k in ks]
            return [sum(p[i] for p in per) / len(ks) for i in range(3)]

        def rmse(est):
            Oe = np.real(np.fft.ifftn(est * meas))
            return np.linalg.norm(Oe - Ob) / np.linalg.norm(Ob)

        print()
        print(f"  --- contrast {c:.2f}  (a dn/lambda = "
              f"{4.8*PHANTOM_SCALE*c*np.abs(dn1).max()/cfg.lam:.3f})   RMSE vs band-limited truth")
        print(f"  {'design':18s} {'fusion':6s} {'A':>7s} {'B':>7s} {'C':>7s} {'D':>7s} {'exp06 D':>8s}"
              f" | {'Bx':>7s} {'Cx':>7s} {'Dx':>7s}")
        for name, ks, hows in designs:
            t = 1.0 / len(ks)
            for how in hows:
                r = {m: rmse(fuse([chan(m, k, t) for k in ks], Sc, how)[0]) for m in MODELS}
                for m, e in zip(("Bx", "Cx", "Dx"), exact(ks, t, how)):
                    r[m] = rmse(e)
                out[(c, name, how)] = r
                e6 = EXP06.get(c, {}).get((name, how))
                if e6 is not None:
                    for m, v in zip("AB", e6[:2]):
                        if v is not None:
                            assert abs(r[m] - v) < 6e-5, f"exp06 {m} not reproduced ({r[m]} vs {v})"
                d6 = f"{e6[2]:8.4f}" if e6 is not None else f"{'':8s}"
                print(f"  {name:18s} {how:6s} {r['A']:7.4f} {r['B']:7.4f} {r['C']:7.4f} {r['D']:7.4f}"
                      f" {d6} | {r['Bx']:7.4f} {r['Cx']:7.4f} {r['Dx']:7.4f}", flush=True)
        print(f"  (A and B reproduce exp06 where it printed them)  [{time.time()-t0:.0f} s]")
        del numS, numN, J1, cw, T, bX

    print()
    print(RULE)
    print("THE MARGINS: fusion vs equal-weight averaging, best design vs the published single 0.8 um")
    print(RULE)
    key = "{0.8,0.55} x + y"
    for c in CONTRASTS:
        print(f"  contrast {c:.2f}:")
        for m, lab in (("A", "A ideal linear"), ("B", "B Born, nearest voxel"),
                       ("Bx", "B Born, exact grid"), ("C", "C MSB, nearest voxel"),
                       ("Cx", "C MSB, exact grid"), ("Dx", "D splice, exact grid")):
            f = out[(c, key, "wiener")][m]
            print(f"    {lab:22s} best design {f:7.4f} | fusion beats averaging "
                  f"{out[(c, key, 'naive')][m]/f:5.2f}x | beats single 0.8 "
                  f"{out[(c, '0.8 x', 'wiener')][m]/f:5.2f}x")
    import json
    with open(out_path("exp11_results.json"), "w") as fh:
        json.dump({f"{c}|{n}|{h}": r for (c, n, h), r in out.items()}, fh, indent=1)
    plot_margins(out, 0.5, 4.8 * PHANTOM_SCALE * 0.5 * np.abs(dn1).max() / cfg.lam)
    print(f"  wrote {out_path('exp11_results.json')} and {out_path('realistic_margins.png')}")
    print(f"  total {time.time()-t_all:.0f} s")


def plot_margins(out, c, adn):
    """The README figure: both margins at one contrast, per model and estimator."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    key = "{0.8,0.55} x + y"
    NEAR, EXACT, IDEAL = "#f47230", "#257de0", "#8a8a8a"   # validated pair + neutral
    TXT, TXT2 = "#191919", "#666666"
    rows = [("Ideal linear model\n(data on the voxels)", [("A", IDEAL, None)]),
            ("First Born, real\ncontinuum data", [("B", NEAR, "nearest voxel"),
                                                  ("Bx", EXACT, "exact column solve")]),
            ("+ multiple scattering\n(multi-slice Born)", [("C", NEAR, "nearest voxel"),
                                                          ("Cx", EXACT, "exact column solve")])]
    panels = (("Best design vs the published single 0.8 µm shear",
               lambda m: out[(c, "0.8 x", "wiener")][m] / out[(c, key, "wiener")][m]),
              ("Per-frequency fusion vs equal-weight averaging",
               lambda m: out[(c, key, "naive")][m] / out[(c, key, "wiener")][m]))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9), sharey=True)
    h, gap = 0.34, 0.03
    for ax, (title, ratio) in zip(axes, panels):
        yt = []
        for i, (lab, bars) in enumerate(rows):
            y0 = -i
            yt.append(y0)
            offs = [0.0] if len(bars) == 1 else [(h + gap) / 2, -(h + gap) / 2]
            for (m, col, _), dy in zip(bars, offs):
                v = ratio(m)
                ax.barh(y0 + dy, v, height=h, color=col, edgecolor="white", linewidth=0)
                ax.text(v + 0.15, y0 + dy, f"{v:.1f}×", va="center", ha="left",
                        fontsize=9, color=TXT)
        ax.axvline(1.0, color=TXT2, lw=0.8, ls=(0, (3, 3)))
        ax.text(1.05, -2.62, "no advantage", fontsize=7.5, color=TXT2, va="bottom")
        ax.set_yticks(yt)
        ax.set_yticklabels([r[0] for r in rows], fontsize=9, color=TXT)
        ax.set_title(title, fontsize=10, color=TXT, loc="left")
        ax.set_xlim(0, 16.5)
        ax.set_ylim(-2.7, 0.55)
        ax.set_xlabel("ratio of reconstruction errors (higher = larger advantage)",
                      fontsize=8.5, color=TXT2)
        ax.tick_params(colors=TXT2, labelsize=8.5)
        ax.grid(axis="x", color="#e6e6e6", lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color("#bdbdbd")
    handles = [plt.Rectangle((0, 0), 1, 1, color=col) for col in (IDEAL, NEAR, EXACT)]
    fig.legend(handles, ["on-voxel data (no gridding needed)",
                         "nearest-voxel gridding (as in exp06)",
                         "exact column solve (exp10)"],
               loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle(f"Multi-shear GODT design margins at a·δn = {adn:.2f} λ, equal photon "
                 f"budget (cell phantom, 180 illuminations)", fontsize=10.5, color=TXT, x=0.01,
                 ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(out_path("realistic_margins.png"), dpi=150, facecolor="white")


if __name__ == "__main__":
    main()
