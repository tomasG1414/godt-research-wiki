"""exp10 -- How much of exp06's gridding loss is real?

exp06 found that moving from voxel-centred data (model A) to continuum-K_z
data with the real pupil (model B, exact first Born) multiplies the best
design's RMSE by 6.2, before any multiple scattering. Its estimator deposits
every ray at the NEAREST K_z voxel (assumption A3 of the design note). This
script replaces that one step and changes nothing else.

Along a K_perp column the continuum datum is an exact DTFT of the nz slices of
the 2D-FFT'd object, so the MAP estimate under exp06's own prior is a small
Toeplitz-plus-circulant solve per column (`shearlib.ColumnWiener`) -- no
interpolation kernel, no new approximation. With every ray placed at its
voxel's K_z it IS exp06's per-voxel rule (shearlib selftest, to 1e-12), so
the comparison below isolates the K_z rounding.

Everything else is exp06's stage 4 verbatim: same phantom (SEED), same 180
illuminations snapped to the propagation lattice, same pupil, prior, rho,
per-ray noise draws and equal-photon shares. The nearest-voxel columns
therefore reproduce exp06's A and B numbers, which the script checks.

Establishes (wiki/notes/gridding-and-forward-model.md):
  1. how much of the A -> B factor is gridding, i.e. removable;
  2. the design margins (fusion vs averaging, best pair vs single 0.8 um)
     on realistic continuum data once it is removed.

Deterministic (SEED, SEED-derived noise). ~10 min, ~6 GB.
Run: python exp10_exact_gridding.py
"""
import time
import numpy as np

from shearlib import (GODT_SIM, KGrid, hermitian_sym, snap_illumination,
                      cap_geometry, born_rytov_ratio, meas_scale, cell_phantom,
                      shell_spectrum, hsym_c, fuse, ColumnWiener)
from exp06_multiple_scattering import (SEED, RHO, N_ILL, K_LO, PHANTOM_SCALE, PAD,
                                       shear_factor, channel_meas, bin_add)

CONTRAST = 0.5          # exp06's first row, a dn / lambda = 0.229
RULE = "=" * 86
# exp06 stage 4, contrast 0.5 (wiki/notes/multiple-scattering-check.md section 5)
EXP06 = {("ODT benchmark", "wiener"): (0.0221, 0.1067), ("0.1 x", "wiener"): (0.2172, 0.2395),
         ("0.8 x", "wiener"): (0.2453, 0.2651), ("0.8 x + 0.8 y", "wiener"): (0.0408, 0.1120),
         ("0.8 x + 0.8 y", "naive"): (0.1804, 0.2090),
         ("{0.8,0.55} x + y", "wiener"): (0.0171, 0.1058),
         ("{0.8,0.55} x + y", "naive"): (0.1642, 0.1952)}


def main():
    t_all = time.time()
    rng = np.random.default_rng(SEED)
    cfg = GODT_SIM
    grid = KGrid(128, 0.1, 128)
    Kabs = grid.Kabs()
    k0, nmed = 2 * np.pi / cfg.lam, cfg.n0
    KI, th_eff = snap_illumination(cfg, grid, N_ILL, dk=2 * np.pi / (PAD * grid.dx))
    caps = [cap_geometry(cfg, grid, ki) for ki in KI]
    sc = meas_scale(cfg, grid)

    O1, dn1 = cell_phantom(grid, rng, scale=PHANTOM_SCALE)
    cover = np.zeros(grid.size, np.int32)
    for cap in caps:
        ixc, iyc = np.nonzero(cap["ok"])
        cover += np.bincount((ixc * grid.n + iyc) * grid.nz + cap["iz"][ixc, iyc],
                             minlength=grid.size).astype(np.int32)
    supp = hermitian_sym(cover.reshape(grid.shape)) > 0
    ref = supp & (Kabs >= 8) & (Kabs <= 24)
    meas = supp & (Kabs >= K_LO)

    dnc = CONTRAST * dn1
    Oc = k0 ** 2 * (nmed ** 2 - (nmed + dnc) ** 2)
    Otc = np.fft.fftn(Oc)
    Sc, _, _ = shell_spectrum(Otc, Kabs, cfg.kmax, grid.dk)
    sig0 = np.sqrt(np.mean(np.abs(Otc[ref]) ** 2) / RHO)
    Ob = np.real(np.fft.ifftn(Otc * meas))

    CH = {"odt": None, "0.1x": (0.1, 0.0), "0.8x": (0.8, 0.0),
          "0.8y": (0.0, 0.8), "0.55x": (0.55, 0.0), "0.55y": (0.0, 0.55)}
    designs = [("ODT benchmark", ["odt"], ["wiener"]),
               ("0.1 x", ["0.1x"], ["wiener"]),
               ("0.8 x", ["0.8x"], ["wiener"]),
               ("0.8 x + 0.8 y", ["0.8x", "0.8y"], ["wiener", "naive"]),
               ("{0.8,0.55} x + y", ["0.8x", "0.55x", "0.8y", "0.55y"], ["wiener", "naive"])]

    print(RULE)
    print("SETUP (exp06 stage 4, contrast %.2f, a dn/lambda = %.3f)" % (
        CONTRAST, 4.8 * PHANTOM_SCALE * CONTRAST * np.abs(dn1).max() / cfg.lam))
    print(RULE)
    print(f"  {grid!r}; {N_ILL} illuminations, realised inclination "
          f"{th_eff.mean():.2f} +- {th_eff.std():.2f} deg; rho = {RHO:.0f}")

    # ------------------------------------------------- acquire, exp06's way
    # per-voxel accumulators (exp06) and per-ray records (for the column solver)
    numS = {(m, k): np.zeros(grid.size, np.complex128) for m in "AB" for k in CH}
    numN = {k: np.zeros(grid.size, np.complex128) for k in CH}
    J1 = {k: np.zeros(grid.size) for k in CH}
    rays = {"cols": [], "Kz": [], "Kv": [], "c": []}
    per_ray = {k: {"w": [], "gA": [], "gB": [], "gN": []} for k in CH}
    t0 = time.time()
    for j, (ki, cap) in enumerate(zip(KI, caps)):
        ixj, iyj = np.nonzero(cap["ok"])
        lin = (ixj * grid.n + iyj) * grid.nz + cap["iz"][ixj, iyj]
        phase = np.exp(-1j * cap["Kz_vox"][ixj, iyj] * cap["z_c"])
        hbase = cfg.k_m / cap["kz"][ixj, iyj]
        rays["cols"].append(ixj * grid.n + iyj)
        rays["Kz"].append(cap["Kz"][ixj, iyj])
        rays["Kv"].append(cap["Kz_vox"][ixj, iyj])
        rays["c"].append(phase)
        rb = born_rytov_ratio(Oc, cfg, grid, ki, cap)
        lgB = np.log(rb)
        gen = np.random.default_rng((SEED, j))          # exp06's per-illumination noise
        for k, dvec in CH.items():
            h = shear_factor(grid, dvec)[ixj, iyj] * hbase
            w = np.abs(h) ** 2 / sig0 ** 2
            J1[k] += np.bincount(lin, weights=w, minlength=grid.size)
            nz_ = gen.standard_normal((2, lin.size))
            noise = np.sqrt(sig0 ** 2 / 2) * (nz_[0] + 1j * nz_[1])
            gN = np.conj(h) * noise / sig0 ** 2
            bin_add(numN[k], lin, gN)
            yA = h * Otc.reshape(-1)[lin]
            yB = (sc * channel_meas(lgB, rb, grid, dvec))[ixj, iyj] * phase
            gA, gB = np.conj(h) * yA / sig0 ** 2, np.conj(h) * yB / sig0 ** 2
            bin_add(numS[("A", k)], lin, gA)
            bin_add(numS[("B", k)], lin, gB)
            for key, v in (("w", w), ("gA", gA), ("gB", gB), ("gN", gN)):
                per_ray[k][key].append(v)
        if j == 0 or (j + 1) % 60 == 0:
            print(f"    illumination {j+1}/{N_ILL}  [{time.time()-t0:.0f} s]", flush=True)
    rays = {k: np.concatenate(v) for k, v in rays.items()}
    for k in CH:
        per_ray[k] = {key: np.concatenate(v) for key, v in per_ray[k].items()}
    print(f"  {rays['cols'].size} rays; K_z rounding |K_z - K_vox|: "
          f"rms {np.sqrt(np.mean((rays['Kz']-rays['Kv'])**2)):.3f}, "
          f"max {np.abs(rays['Kz']-rays['Kv']).max():.3f} rad/um (voxel {grid.dk:.3f}); "
          f"phase error at the object's edge (|z| = 3 um): rms "
          f"{3*np.sqrt(np.mean((rays['Kz']-rays['Kv'])**2)):.2f} rad")

    # ------------------------------------------------- the column solver
    t0 = time.time()
    cw = ColumnWiener(grid, rays["cols"], rays["Kz"], rays["c"])
    T = {k: cw.data_term(per_ray[k]["w"]) for k in CH}
    bB = {k: cw.rhs(per_ray[k]["gB"]) for k in CH}
    bN = {k: cw.rhs(per_ray[k]["gN"]) for k in CH}
    print(f"  column solver built: {cw.ucols.size} columns x {grid.nz} slices "
          f"[{time.time()-t0:.0f} s]")

    def rmse(est):
        Oe = np.real(np.fft.ifftn(est * meas))
        return np.linalg.norm(Oe - Ob) / np.linalg.norm(Ob)

    def chan(m, k, t):
        num = (t * numS[(m, k)] + np.sqrt(t) * numN[k]).reshape(grid.shape)
        return hsym_c(num), hermitian_sym((t * J1[k]).reshape(grid.shape))

    def exact(ks, t, how):
        if how == "wiener":
            return cw.solve([t * T[k] for k in ks],
                            [t * bB[k] + np.sqrt(t) * bN[k] for k in ks], Sc)
        ests = [cw.solve([t * T[k]], [t * bB[k] + np.sqrt(t) * bN[k]], Sc) for k in ks]
        return sum(ests) / len(ests)

    print()
    print(RULE)
    print("1. NEAREST-VOXEL vs EXACT K_z, same data, same prior, same noise")
    print(RULE)
    print(f"  {'design':20s} {'fusion':7s} {'A ideal':>9s} {'B nearest':>10s} {'B exact':>9s}"
          f" {'exp06 A':>8s} {'exp06 B':>8s}")
    res = {}
    t0 = time.time()
    for name, ks, hows in designs:
        t = 1.0 / len(ks)
        for how in hows:
            rA = rmse(fuse([chan("A", k, t) for k in ks], Sc, how)[0])
            rB = rmse(fuse([chan("B", k, t) for k in ks], Sc, how)[0])
            rX = rmse(exact(ks, t, how))
            res[(name, how)] = (rA, rB, rX)
            eA, eB = EXP06[(name, how)]
            print(f"  {name:20s} {how:7s} {rA:9.4f} {rB:10.4f} {rX:9.4f} {eA:8.4f} {eB:8.4f}",
                  flush=True)
            assert abs(rA - eA) < 6e-5 and abs(rB - eB) < 6e-5, "exp06 not reproduced"
    print(f"  (the A and B-nearest columns reproduce exp06 to its 4 printed digits) "
          f"[{time.time()-t0:.0f} s]")

    key, pair = "{0.8,0.55} x + y", "0.8 x + 0.8 y"
    print()
    print(RULE)
    print("2. WHAT IS LEFT OF THE A -> B FACTOR, AND OF THE MARGINS")
    print(RULE)
    for i, lab in ((0, "A ideal-linear"), (1, "B, nearest voxel"), (2, "B, exact K_z")):
        f = res[(key, "wiener")][i]
        print(f"  {lab:18s} best design {f:.4f} | fusion beats averaging "
              f"{res[(key, 'naive')][i]/f:5.1f}x | best pair beats single 0.8 "
              f"{res[('0.8 x', 'wiener')][i]/f:5.1f}x | 0.8 x+y fusion beats averaging "
              f"{res[(pair, 'naive')][i]/res[(pair, 'wiener')][i]:4.1f}x")
    a, b, x = res[(key, "wiener")]
    print(f"  A -> B factor on the best design: {b/a:.2f} with nearest-voxel deposition, "
          f"{x/a:.2f} with the exact estimator,")
    print(f"  so {100*(1-np.log(x/a)/np.log(b/a)):.0f} % of it (in log terms) was the K_z rounding.")
    print(f"  total {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
