"""exp03 -- 3D Ewald sampling: Fisher-information maps, lost volume, cone overlap.

Establishes (design note T1, T2):
  1. Separability  J(K) = rho * S_t(K_perp) * M(K): the shear design enters
     only through the 2D shear power S_t, the tomography only through the
     cap density M(K), which has a closed form (shearlib.cap_density_closed)
     verified here against direct cap binning.
  2. P1: the direction-independent blind set {K_perp = 0} is never sampled by
     transmission illumination scanning except at K = 0 -- it sits inside the
     missing cone, so with two non-parallel shears the extra loss over plain
     ODT is confined to the m != 0 null planes. The direction-DEPENDENT blind
     plane K.d = 0 of a single shear is NOT inside the cone.
  3. True 3D fractions of the sampled support lost at the GODT shear values,
     decomposed into null comb / blind plane / low SNR, and what each fusion
     design at equal photon budget recovers.
  4. Angular undersampling: how much of the support 180 projections reach on
     a 25.6 um and a 12.8 um field (nearest-neighbour direct inversion).

Writes: out/ewald_maps.png, out/lost_volume.png
Wiki:   wiki/notes/multishear-fusion-design-solution.md

Run: python exp03_ewald_fisher.py       (~1-2 min; caches M(K) in out/)
"""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, KGrid, G2_vec, hermitian_sym, cap_density,
                      cap_density_closed, cached_density, nulls_in_band, out_path)

RHO = 100.0          # per-frame, per-sample SNR at unit shear power (20 dB), flat prior
N_AZ = 1440          # dense azimuths for the information DENSITY (exposure = 180 frames)
N_FRAMES = 180       # GODT's acquisition
RULE = "=" * 78

# designs at EQUAL photon budget: list of ((dx, dy), time share)
X, Y = (1.0, 0.0), (0.0, 1.0)


def design(label, items):
    if items is None:
        return label, None
    tot = sum(t for _, t in items)
    return label, [((d * v[0], d * v[1]), t / tot) for (d, v), t in items]


DESIGNS = [
    design("ODT (|G|^2 = 1, benchmark)", None),
    design("0.1 x  (GODT sim shear)", [((0.1, X), 1)]),
    design("0.33 x (preprint exp)", [((0.33, X), 1)]),
    design("0.8 x  (published exp)", [((0.8, X), 1)]),
    design("0.8 x + 0.4 x  (ratio 2)", [((0.8, X), 1), ((0.4, X), 1)]),
    design("0.8 x + 0.55 x", [((0.8, X), 1), ((0.55, X), 1)]),
    design("0.8 x + 0.08 x", [((0.8, X), 1), ((0.08, X), 1)]),
    design("0.8 x + 0.8 y", [((0.8, X), 1), ((0.8, Y), 1)]),
    design("{0.8, 0.55} x + {0.8, 0.55} y", [((0.8, X), 1), ((0.55, X), 1),
                                            ((0.8, Y), 1), ((0.55, Y), 1)]),
]

# Plain-language labels for the README figure (the printed table keeps the short ones).
FIG_LABELS = {
    "ODT (|G|^2 = 1, benchmark)": "ODT, no shear (reference)",
    "0.1 x  (GODT sim shear)": "0.1 µm along x (paper's simulations)",
    "0.33 x (preprint exp)": "0.33 µm along x (preprint's shear)",
    "0.8 x  (published exp)": "0.8 µm along x (published shear)",
    "0.8 x + 0.4 x  (ratio 2)": "0.8 + 0.4 µm along x (ratio 2)",
    "0.8 x + 0.55 x": "0.8 + 0.55 µm along x",
    "0.8 x + 0.08 x": "0.8 + 0.08 µm along x",
    "0.8 x + 0.8 y": "0.8 µm along x and along y",
    "{0.8, 0.55} x + {0.8, 0.55} y": "0.8 + 0.55 µm along x and along y",
}


def shear_power(grid, items):
    if items is None:
        return np.ones_like(grid.KX)
    return sum(t * G2_vec(grid.KX, grid.KY, d) for d, t in items)


def main():
    cfg = GODT_SIM
    grid = KGrid(256, 0.1, 96)
    grid.check(cfg)
    KX, KY, KZ = grid.K3()
    KP = np.hypot(KX, KY)

    # ------------------------------------------------ 1. support and density
    print(RULE)
    print("1. SAMPLED SUPPORT AND CAP DENSITY  (GODT simulation geometry)")
    print(RULE)
    print(f"  {cfg!r}")
    print(f"  {grid!r}   field {grid.n*grid.dx:.1f} um")
    lo, hi = cfg.Kz_range
    print(f"  K_z range of the caps: [{lo:.2f}, {hi:.2f}] rad/um  "
          f"(axial Nyquist pitch {np.pi/max(abs(lo),hi):.3f} um)")
    t0 = time.time()
    M = cached_density(cfg, grid, N_AZ, N_FRAMES)         # Hermitian-symmetrised
    print(f"  M(K) with {N_AZ} azimuths, {N_FRAMES}-frame exposure: {time.time()-t0:.1f} s")
    supp = M > 0
    nsupp = int(supp.sum())
    band_cyl = (KP <= cfg.kmax) & (KZ >= lo - grid.dk) & (KZ <= hi + grid.dk)
    print(f"  support: {nsupp} voxels = {100*nsupp/band_cyl.sum():.1f} % of the cylinder "
          f"|K_perp| <= kmax, K_z in range")
    print(f"  M on support: median {np.median(M[supp]):.2f}, mean {M[supp].mean():.2f}, "
          f"max {M[supp].max():.1f}  (weighted caps per voxel)")

    # closed form check
    Mc = hermitian_sym(cap_density_closed(cfg, grid, N_FRAMES))
    th = np.deg2rad(cfg.theta_ill)
    with np.errstate(divide="ignore", invalid="ignore"):
        cc = -(KP ** 2 + KZ ** 2 + 2 * cfg.k_m * KZ * np.cos(th)) / (2 * cfg.k_m * KP * np.sin(th))
    safe = (np.abs(cc) < 0.8) & (Mc > 0) & (M > 0) & (KP > 3)
    dev = M[safe] / Mc[safe] - 1
    print(f"  closed form vs binning ({N_AZ} az.): total {M.sum()/Mc.sum():.4f}; voxel-wise away "
          f"from the caustic (|c|<0.8, {safe.sum()} voxels): median |dev| "
          f"{np.median(np.abs(dev)):.3f}, 90th pct {np.percentile(np.abs(dev), 90):.3f}")
    print(f"  support agreement: {100*np.mean((M>0)==(Mc>0)):.2f} % of voxels")
    # the voxel-wise scatter is the K_z rounding of the digital caps; it shrinks
    # with denser azimuths (checked on the 128 grid, where 11520 azimuths are cheap)
    g2 = KGrid(128, 0.1, 48)
    Mb = hermitian_sym(cap_density(cfg, g2, 11520, N_FRAMES))
    Mc2 = hermitian_sym(cap_density_closed(cfg, g2, N_FRAMES))
    KX2, KY2, KZ2 = g2.K3()
    KP2 = np.hypot(KX2, KY2)
    with np.errstate(divide="ignore", invalid="ignore"):
        cc2 = -(KP2 ** 2 + KZ2 ** 2 + 2 * cfg.k_m * KZ2 * np.cos(th)) / (2 * cfg.k_m * KP2 * np.sin(th))
    safe2 = (np.abs(cc2) < 0.8) & (Mc2 > 0) & (Mb > 0) & (KP2 > 3)
    dev2 = Mb[safe2] / Mc2[safe2] - 1
    print(f"  same on a 128^2x48 grid with 11520 azimuths: total {Mb.sum()/Mc2.sum():.4f}; "
          f"voxel-wise median |dev| {np.median(np.abs(dev2)):.3f}, 90th pct "
          f"{np.percentile(np.abs(dev2), 90):.3f}, mean signed dev {dev2.mean():+.3f}")

    # missing-cone boundary from the support vs the closed form |c| = 1
    km, s_, c_ = cfg.k_m, np.sin(th), np.cos(th)

    def kp_min_analytic(kz):
        """min |K_perp| on the symmetrised support at K_z: the cap at scattered
        k_z = K_z + k_m cos th has |k_perp| = r, so |K_perp| in [|r - k_m sin th|, r + k_m sin th]."""
        best = np.inf
        for z in (kz, -kz):                      # Hermitian partner
            if z + km * c_ < cfg.kz_floor:
                continue
            disc = (km * s_) ** 2 - z ** 2 - 2 * km * z * c_
            if disc > 0:
                best = min(best, abs(km * s_ - np.sqrt(disc)))
        return best
    print()
    print("  inner boundary of the symmetrised support (the missing cone): smallest sampled")
    print("  |K_perp| per K_z, binned vs closed form |c(K)| = 1. Near the origin the boundary")
    print("  is K_z = +-K_perp tan(60 deg): half-angle 30 deg; it narrows at larger |K_z|.")
    for kzv in (1.0, 2.0, 4.0, 6.0, 8.0):
        iz = int(np.rint(kzv / grid.dk))
        sl = supp[:, :, iz]
        kp_min = KP[:, :, iz][sl].min() if sl.any() else np.nan
        print(f"    K_z = {grid.kz[iz]:5.2f}: min |K_perp| binned {kp_min:5.2f}, closed form "
              f"{kp_min_analytic(grid.kz[iz]):5.2f}  -> half-angle "
              f"{np.degrees(np.arctan2(kp_min_analytic(grid.kz[iz]), grid.kz[iz])):5.1f} deg")

    # -------------------------------------------------------- 2. P1: the axis
    print()
    print(RULE)
    print("2. P1 -- IS THE UNIVERSAL SHEAR BLIND SET INSIDE THE MISSING CONE?")
    print(RULE)
    axis = M[0, 0, :]
    print(f"  samples on the K_z axis (K_perp = 0): at K_z = 0: {axis[0]:.1f};  "
          f"max over K_z != 0: {axis[1:].max():.1f}")
    print("  => exactly zero: a cap reaches K_perp = 0 only where k_perp = k_i,perp, i.e. K = 0.")
    print("     The direction-independent blind set of ANY in-plane shear set is unmeasured")
    print("     by transmission ODT anyway.  (Exact; the digital caps confirm it.)")
    r_ax = 1.0
    near_axis = supp & (KP < r_ax) & (KP > 0)
    print(f"  support within |K_perp| < {r_ax} rad/um of the axis: "
          f"{100*near_axis.sum()/nsupp:.3f} % of the support")

    # ---------------------------------------------- 3. angular undersampling
    print()
    print(RULE)
    print("3. ANGULAR COVERAGE AT FINITE PROJECTION COUNT (nearest-neighbour DI)")
    print(RULE)
    for n, dxg, nz in ((256, 0.1, 96), (128, 0.1, 48)):
        g = KGrid(n, dxg, nz)
        Md = hermitian_sym(cap_density(cfg, g, 2880, N_FRAMES))
        ref = Md > 0
        row = []
        for naz in (180, 360, 720, 1440):
            Mn = hermitian_sym(cap_density(cfg, g, naz, N_FRAMES))
            row.append(100 * np.mean((Mn > 0)[ref]))
        print(f"  field {n*dxg:5.1f} um (dk={g.dk:.3f}): support hit at 180/360/720/1440 proj: "
              + " / ".join(f"{r:5.1f} %" for r in row) + "   (ref: 2880)")
    print("  => 180 projections leave part of the reachable support unvisited on a 25.6 um")
    print("     field; the shear design is orthogonal to this (all channels share the caps),")
    print("     so the rest of this study uses the dense-azimuth information density.")

    # ----------------------------------------------------- 4. lost volume
    print()
    print(RULE)
    print(f"4. 3D VOLUME LOST TO THE SHEAR NULLS  (rho = {RHO:.0f} per frame-sample, flat prior,")
    print("   equal photon budget: time shares sum to 1)")
    print(RULE)
    print(f"  {'design':34s} {'%S_t<0.04':>10} {'mean relMSE':>12} {'%relMSE>0.25':>13} {'%relMSE>0.5':>12}")
    results = {}
    for label, items in DESIGNS:
        St = shear_power(grid, items)
        rel = 1.0 / (RHO * M * St[:, :, None] + 1.0)
        r = rel[supp]
        null3 = np.broadcast_to(St[:, :, None] < 0.04, supp.shape)
        results[label] = (100 * np.mean(null3[supp]),
                          r.mean(), 100 * np.mean(r > 0.25), 100 * np.mean(r > 0.5))
        print(f"  {label:34s} {results[label][0]:10.2f} {results[label][1]:12.4f} "
              f"{results[label][2]:13.2f} {results[label][3]:12.2f}")

    # decomposition for the published shear
    d = 0.8
    St = shear_power(grid, dict(DESIGNS)["0.8 x  (published exp)"])
    rel = 1.0 / (RHO * M * St[:, :, None] + 1.0)
    lost = supp & (rel > 0.25)
    half = 0.5
    comb = np.zeros_like(KX, dtype=bool)
    for nu in nulls_in_band(d, cfg.kmax):
        comb |= np.abs(np.abs(KX) - nu) < half
    blind = np.abs(KX) < half
    n_lost = lost.sum()
    print()
    print(f"  decomposition of the lost voxels for 0.8 x ({n_lost} voxels, "
          f"{100*n_lost/nsupp:.2f} % of support):")
    print(f"    within {half} rad/um of an m != 0 null plane : {100*(lost & comb).sum()/n_lost:5.1f} %")
    print(f"    within {half} rad/um of the blind plane K_x=0: {100*(lost & blind & ~comb).sum()/n_lost:5.1f} %")
    print(f"    elsewhere (low SNR, support edge)            : {100*(lost & ~comb & ~blind).sum()/n_lost:5.1f} %")
    print(f"  ODT loses (same criterion, no shear factor)    : "
          f"{100*np.mean((1/(RHO*M+1))[supp] > 0.25):.2f} % of support")

    # --------------------------------------- 5. blind plane vs missing cone
    print()
    print(RULE)
    print("5. THE DIRECTION-DEPENDENT BLIND PLANE IS NOT IN THE MISSING CONE")
    print(RULE)
    plane = supp & (KX == 0)
    near = supp & (np.abs(KX) < 1.0)
    print(f"  support voxels ON the plane K_x = 0 (blind for every x-shear):  "
          f"{plane.sum()} = {100*plane.sum()/nsupp:.2f} % of support")
    print(f"  support within |K_x| < 1 rad/um:                                 "
          f"{100*near.sum()/nsupp:.2f} % of support")
    print(f"  of the K_x = 0 plane's support, on the K_z axis (in the cone):    "
          f"{(plane & (KY == 0)).sum()} voxels (K = 0 only)")
    print("  => a single shear direction is blind on a whole measured plane; the second")
    print("     direction is what pays for it, and its residual blind set (the axis) is free.")

    # ----------------------------------------------------------- figures
    def sl_ky0(A):   # K_x (horizontal) vs K_z (vertical) at K_y = 0
        return np.fft.fftshift(A[:, 0, :]).T
    def sl_kz0(A):
        return np.fft.fftshift(A[:, :, 0]).T
    def sl_kx0(A):   # K_y vs K_z at K_x = 0
        return np.fft.fftshift(A[0, :, :]).T
    ext_xz = [grid.kx.min(), -grid.kx.min(), grid.kz.min(), -grid.kz.min()]
    ext_xy = [grid.kx.min(), -grid.kx.min(), grid.ky.min(), -grid.ky.min()]

    fig, ax = plt.subplots(2, 3, figsize=(15, 7.2))
    with np.errstate(divide="ignore"):
        L = np.log10(M)
    for a, (img, ext, ttl, xl, yl) in zip(ax[0], [
            (sl_ky0(L), ext_xz, "log10 M(K): K_y = 0 plane", "K_x", "K_z"),
            (sl_kz0(L), ext_xy, "log10 M(K): K_z = 0 plane", "K_x", "K_y"),
            (sl_kx0(L), ext_xz, "log10 M(K): K_x = 0 plane", "K_y", "K_z")]):
        im = a.imshow(img, origin="lower", extent=ext, aspect="auto", cmap="viridis",
                      vmin=-1, vmax=2)
        a.set(title=ttl, xlabel=f"{xl} (rad/um)", ylabel=f"{yl} (rad/um)")
        plt.colorbar(im, ax=a, fraction=0.046)
    St8 = shear_power(grid, dict(DESIGNS)["0.8 x  (published exp)"])
    Stxy = shear_power(grid, dict(DESIGNS)["{0.8, 0.55} x + {0.8, 0.55} y"])
    rel8 = 1.0 / (RHO * M * St8[:, :, None] + 1.0)
    relxy = 1.0 / (RHO * M * Stxy[:, :, None] + 1.0)
    rel8[~supp] = np.nan
    relxy[~supp] = np.nan
    for a, (img, ext, ttl, xl, yl) in zip(ax[1], [
            (sl_ky0(rel8), ext_xz, "rel. MSE, 0.8 x: K_y = 0 (null planes)", "K_x", "K_z"),
            (sl_kx0(rel8), ext_xz, "rel. MSE, 0.8 x: K_x = 0 (blind plane)", "K_y", "K_z"),
            (sl_ky0(relxy), ext_xz, "rel. MSE, {0.8,0.55} x+y: K_y = 0", "K_x", "K_z")]):
        im = a.imshow(img, origin="lower", extent=ext, aspect="auto", cmap="magma", vmin=0, vmax=1)
        a.set(title=ttl, xlabel=f"{xl} (rad/um)", ylabel=f"{yl} (rad/um)")
        plt.colorbar(im, ax=a, fraction=0.046)
    for a in ax.ravel():
        a.set_ylim(-9.5, 9.5)
    ax[0, 1].set_ylim(-32, 32)
    fig.suptitle(f"GODT simulation geometry, 180 projections (density from {N_AZ} azimuths); "
                 f"rho = {RHO:.0f}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path("ewald_maps.png"), dpi=110)

    # README figure: one panel, large type, so it stays legible at phone width
    labels = [l for l, _ in DESIGNS]
    lost = [results[l][2] for l in labels]
    with plt.rc_context({"font.size": 13}):
        fig2, ax2 = plt.subplots(figsize=(8.6, 5.4))
        ax2.barh([FIG_LABELS[l] for l in labels], lost, color="C3")
        for i, v in enumerate(lost):
            ax2.text(v + 0.04, i, f"{v:.2f} %", va="center", fontsize=12)
        ax2.set_xlim(0, 1.18 * max(lost))
        ax2.set_xlabel("% of the measurable 3D frequencies lost (relative MSE > 0.25)")
        ax2.set_title("Simulation: same total light, 20 dB, GODT simulation optics",
                      fontsize=13)
        ax2.invert_yaxis()
        ax2.grid(alpha=.3, axis="x")
        fig2.tight_layout()
        fig2.savefig(out_path("lost_volume.png"), dpi=150, bbox_inches="tight")
    print()
    print(f"  wrote {out_path('ewald_maps.png')}")
    print(f"  wrote {out_path('lost_volume.png')}")


if __name__ == "__main__":
    main()
