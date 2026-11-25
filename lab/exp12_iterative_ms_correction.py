"""exp12 -- Does a cheap iterative data correction restore the design margins?

DRAFT, NOT YET REVIEWED OR WRITTEN UP. Drafted by a subagent that stopped
before finishing; the `--cbs` option described below is not implemented, and
the wiki note it names does not exist yet. Nothing here is quoted in the wiki.

exp11 found that multiple scattering (multi-slice first Born, MSB, validated
against the exact CBS referee on near-axis rays) collapses the multi-shear
margins at a dn = 0.23 lambda: best design vs single 0.8 um shear 1.3x, fusion
vs averaging 1.1x, best-design RMSE 0.27. This script tests the cheapest
inverse that models multiple scattering: a fixed-point (distorted-data) loop
around the UNCHANGED linear estimator W (exp06's nearest-voxel Wiener rule),

    O_0     = W d
    e_k     = F(O_k) - L(O_k)          per illumination, channel and ray
    O_{k+1} = W (d - e_k)

with F a multiple-scattering forward model and L the linear model F is being
corrected to. Two choices of L:
  ms    L = continuum first Born (the task's suggestion). Fixed point: the
        estimate whose nonlinear excess explains the data's; the estimator's
        own gridding error (exp10) is untouched, so the ceiling is exp11's B.
  full  L = the on-voxel model h O~[voxel] that W inverts exactly (exp11's A).
        Then O_{k+1} = O_k + W (d - F(O_k)): a preconditioned Landweber step
        on the full nonlinear, gridding-free model; the ceiling is the
        gridding-free estimate (exp10's column solve), not B.
Damping alpha mixes e_k with the previous correction.

The inverse crime. Data are exp11's model C: MSB on a 38.4 um window with an
absorbing edge, measured on the whole window (lattice_samples). The loop's F
is a DIFFERENT discretisation of the same physics: `BlochMSB`, multi-slice
Born on the 12.8 um reconstruction grid itself, periodic (the estimate is
only defined there), with the illumination carried as a Bloch phase so it
need not sit on that grid's lattice; optionally with 0.2 um slices. The
PHYSICS MODEL IS SHARED (forward-only multi-slice Born), so every number is
an optimistic bound. Three controls bound how optimistic:
  - 'crime' oracle: subtract the data model's own excess of the TRUE object;
  - 'oracle': subtract BlochMSB's excess of the TRUE object (model mismatch
    only, no estimate error);
  - `--cbs`: one illumination's excess against the exact CBS solver
    (two solves on a 38.4 um window), for the truth and for the loop's final
    estimate.

Everything else is exp11's: phantom (SEED), 180 illuminations, pupil, prior,
rho, noise draws and equal-photon shares; its A, B and C columns are
reproduced and asserted before anything is corrected.

Deterministic. Data generation (3 x 180 MSB runs at 384 x 384) is cached in
lab/out/exp12_data_c*.npz and regenerated when missing.
Wiki: wiki/notes/iterative-multiple-scattering-correction.md
Run:  python exp12_iterative_ms_correction.py           contrasts 0.5, 1.0
      python exp12_iterative_ms_correction.py --quick   contrast 0.5, fewer iterations
      python exp12_iterative_ms_correction.py --cbs     the exact-referee check only
                                                        (needs a finished main run)
"""
import os
import sys
import time
import numpy as np

from shearlib import (GODT_SIM, KGrid, hermitian_sym, snap_illumination, cap_geometry,
                      born_rytov_ratio, meas_scale, cell_phantom, shell_spectrum, hsym_c,
                      fuse, MSB, CBS, ColumnWiener, lattice_samples, shear_rytov, out_path)
from exp06_multiple_scattering import (SEED, RHO, N_ILL, K_LO, PHANTOM_SCALE, PAD, EPS,
                                       shear_factor, channel_meas, bin_add)

QUICK = "--quick" in sys.argv
CONTRASTS = (0.5,) if QUICK else (0.5, 1.0)
N_ITER = 4 if QUICK else 8
BIG = 384               # exp11's data window, 38.4 um
FRAME = 0.1
Z0, NZC = 28, 72
RULE = "=" * 100
# exp06's channel order: the noise generator is drawn in this order
CH6 = {"odt": None, "0.1x": (0.1, 0.0), "0.8x": (0.8, 0.0),
       "0.8y": (0.0, 0.8), "0.55x": (0.55, 0.0), "0.55y": (0.0, 0.55)}
CH = {k: v for k, v in CH6.items() if k not in ("odt", "0.1x")}     # the ones used here
BEST = ["0.8x", "0.55x", "0.8y", "0.55y"]
PIPES = (("best", BEST, "wiener"), ("avg", BEST, "naive"), ("single", ["0.8x"], "wiener"))
# exp11 (wiki/notes/gridding-and-forward-model.md section 4): (A, B, C, Bx, Cx)
EXP11 = {0.5: {"best": (0.0171, 0.1058, 0.2686, 0.0280, 0.2802),
               "avg": (0.1642, 0.1952, 0.3080, 0.1661, 0.3024),
               "single": (0.2453, 0.2651, 0.3572, 0.2464, 0.3663)},
         1.0: {"best": (0.0172, 0.1057, 0.4917, None, 0.5353)}}


# --------------------------------------------------------------- forward models
class BlochMSB:
    """shearlib.MSB on a periodic n x n grid for an illumination that is NOT on
    that grid's lattice.

    The total field is u = exp(i k_i . r_perp) v with v periodic, so the march
    runs on v with every transverse frequency shifted by k_i. Physically: an
    infinite array of copies of the object at the grid's pitch. For k_i on the
    lattice it is shearlib.MSB exactly (checked in `selfcheck`), and its
    linear limit is periodic first Born (born_rytov_ratio) for any k_i.

    `zc` are the slice centres (uniform pitch `dz`), so the same class runs
    re-binned objects.
    """

    def __init__(self, O, dx, zc, dz, k_m, kz_floor_frac=0.1):
        import scipy.fft as sfft
        self.sfft = sfft
        self.O = np.asarray(O)
        self.n, _, self.nz = self.O.shape
        self.dx, self.dz, self.k_m = float(dx), float(dz), float(k_m)
        self.zc = np.asarray(zc, float)
        self.floor = kz_floor_frac * self.k_m
        k = 2 * np.pi * np.fft.fftfreq(self.n, self.dx)
        self.KX, self.KY = np.meshgrid(k, k, indexing="ij")
        self.active = [bool(np.any(self.O[:, :, s] != 0)) for s in range(self.nz)]

    def ratio(self, ki):
        f, fi = self.sfft.fft2, self.sfft.ifft2
        kix, kiy = float(ki[0]), float(ki[1])
        kiz = np.sqrt(self.k_m ** 2 - kix ** 2 - kiy ** 2)
        kt2 = (self.KX + kix) ** 2 + (self.KY + kiy) ** 2
        prop = kt2 < self.k_m ** 2
        kz = np.sqrt(np.maximum(self.k_m ** 2 - kt2, 0.0))
        Ph = np.where(prop, np.exp(1j * kz * self.dz / 2), 0.0)
        P = Ph * Ph
        kern = np.where(prop & (kz >= self.floor),
                        -1j * self.dz / (2 * np.maximum(kz, 1e-12)), 0.0) * Ph
        zlo = self.zc[0] - self.dz / 2
        Vk = np.zeros((self.n, self.n), complex)
        Vk[0, 0] = self.n * self.n * np.exp(1j * kiz * zlo)
        for s in range(self.nz):
            if self.active[s]:
                Vc = fi(Vk * Ph)
                Vk = Vk * P + kern * f(self.O[:, :, s] * Vc)
            else:
                Vk = Vk * P
        zex = self.zc[-1] + self.dz / 2
        return fi(Vk * np.where(prop, np.exp(-1j * kz * zex), 0.0))


def born_UR(Oxy, zc, dz, cap, ix, iy):
    """First-Rytov spectrum -i dz O^(K_perp, K_z) / (2 k_z) on the collected rays
    (ix, iy), at each ray's continuum K_z: born_rytov_ratio before its exp.
    `Oxy` is fft2 of the object over its transverse axes."""
    Kz = cap["Kz"][ix, iy]
    cols = Oxy[ix, iy, :]                              # (rays, nz)
    Ohat = (cols * np.exp(-1j * Kz[:, None] * zc[None, :])).sum(axis=1)
    return -1j * dz * Ohat / (2 * cap["kz"][ix, iy])


def rebin(O, zc, dz, b):
    """b adjacent slices -> one of thickness b dz (mean potential)."""
    if b == 1:
        return O, zc, dz
    nz = O.shape[2] // b
    return (O.reshape(O.shape[0], O.shape[1], nz, b).mean(axis=3),
            zc.reshape(nz, b).mean(axis=1), b * dz)


# ------------------------------------------------------------------ the set-up
class Setup:
    """exp11's geometry, with everything the per-illumination loops reuse."""

    def __init__(self):
        cfg = self.cfg = GODT_SIM
        grid = self.grid = KGrid(128, 0.1, 128)
        self.Kabs = grid.Kabs()
        self.k0, self.nmed = 2 * np.pi / cfg.lam, cfg.n0
        _, self.dn1 = cell_phantom(grid, np.random.default_rng(SEED), scale=PHANTOM_SCALE)
        assert not np.any(self.dn1[:, :, :Z0]) and not np.any(self.dn1[:, :, Z0 + NZC:])
        self.KI, self.th = snap_illumination(cfg, grid, N_ILL, dk=2 * np.pi / (PAD * grid.dx))
        self.caps = [cap_geometry(cfg, grid, ki) for ki in self.KI]
        self.sc = meas_scale(cfg, grid)
        self.zc = (np.arange(grid.nz) - grid.nz // 2) * grid.dz
        cover = np.zeros(grid.size, np.int32)
        self.ill = []
        for cap in self.caps:
            ix, iy = np.nonzero(cap["ok"])
            lin = (ix * grid.n + iy) * grid.nz + cap["iz"][ix, iy]
            cover += np.bincount(lin, minlength=grid.size).astype(np.int32)
            hb = cfg.k_m / cap["kz"][ix, iy]
            self.ill.append({"ix": ix, "iy": iy, "lin": lin,
                             "phase": np.exp(-1j * cap["Kz_vox"][ix, iy] * cap["z_c"]),
                             "sf": {k: shear_factor(grid, d)[ix, iy] for k, d in CH.items()},
                             "hb": hb})
        self.nray = sum(i["lin"].size for i in self.ill)
        self.supp = hermitian_sym(cover.reshape(grid.shape)) > 0
        self.ref = self.supp & (self.Kabs >= 8) & (self.Kabs <= 24)
        self.meas = self.supp & (self.Kabs >= K_LO)

    def potential(self, c):
        return self.k0 ** 2 * (self.nmed ** 2 - (self.nmed + c * self.dn1) ** 2)

    def adn(self, c):
        return 4.8 * PHANTOM_SCALE * c * np.abs(self.dn1).max() / self.cfg.lam


def forward(S, O, sig0, bins=1, lin_model="born", rays=False):
    """Binned back-projections  sum_rays conj(h) y / sig0^2  of BlochMSB's data
    for the object O (real, on the reconstruction grid) and of its linear
    counterpart, per channel.

    Returns (numF, numL[, eray]): numL is continuum first Born of the same
    (re-binned) object -- BlochMSB's exact linear limit. With rays=True also
    the ray-level excess F - L per channel (complex64, concatenated).
    """
    grid, cfg = S.grid, S.cfg
    Ob, zc, dz = rebin(O, S.zc, grid.dz, bins)
    model = BlochMSB(Ob, grid.dx, zc, dz, cfg.k_m)
    Oxy = np.fft.fft2(Ob, axes=(0, 1))
    numF = {k: np.zeros(grid.size, np.complex128) for k in CH}
    numL = {k: np.zeros(grid.size, np.complex128) for k in CH}
    er = {k: [] for k in CH}
    for ki, cap, il in zip(S.KI, S.caps, S.ill):
        r = model.ratio(ki)
        UR = born_UR(Oxy, zc, dz, cap, il["ix"], il["iy"])
        for k, dvec in CH.items():
            yF = (S.sc * np.fft.fft2(shear_rytov(r, grid.dx, dvec)))[il["ix"], il["iy"]] * il["phase"]
            yL = S.sc * il["sf"][k] * UR * il["phase"]
            hc = np.conj(il["sf"][k] * il["hb"]) / sig0 ** 2
            bin_add(numF[k], il["lin"], hc * yF)
            bin_add(numL[k], il["lin"], hc * yL)
            if rays:
                er[k].append((yF - yL).astype(np.complex64))
    if rays:
        return numF, numL, {k: np.concatenate(v) for k, v in er.items()}
    return numF, numL


# ------------------------------------------------------------------- the data
def make_data(S, c, yE=None):
    """exp11's models A, B, C for the four channels used here, binned, plus the
    data model's own linear part (c/eps) C(eps) and ray-level arrays."""
    cfg, grid = S.cfg, S.grid
    Oc = S.potential(c)
    Otc = np.fft.fftn(Oc)
    sig0 = np.sqrt(np.mean(np.abs(Otc[S.ref]) ** 2) / RHO)

    def big_O(cc):
        O = np.zeros((BIG, BIG, NZC))
        o = (BIG - grid.n) // 2
        O[o:o + grid.n, o:o + grid.n] = S.k0 ** 2 * (
            S.nmed ** 2 - (S.nmed + cc * S.dn1[:, :, Z0:Z0 + NZC]) ** 2)
        return O

    def msb_meas(model, ki, il):
        r = model.ratio(ki)
        return {k: S.sc * lattice_samples(shear_rytov(r, grid.dx, d, wrap_safe=True), grid.n,
                                          il["ix"], il["iy"]) * il["phase"] for k, d in CH.items()}

    t0 = time.time()
    if yE is None:
        m = MSB(big_O(EPS), grid.dx, grid.dz, cfg.k_m, absorb_frac=FRAME)
        yE = [msb_meas(m, ki, il) for ki, il in zip(S.KI, S.ill)]
        del m
        print(f"    data model at eps = {EPS}: {N_ILL} illuminations [{time.time()-t0:.0f} s]",
              flush=True)
    msb_c = MSB(big_O(c), grid.dx, grid.dz, cfg.k_m, absorb_frac=FRAME)
    D = {f"num{m}_{k}": np.zeros(grid.size, np.complex128) for m in "ABCLN" for k in CH}
    D.update({f"J_{k}": np.zeros(grid.size) for k in CH})
    ray = {f"{m}_{k}": [] for m in ("yC", "eC", "noise") for k in CH}
    for j, (ki, cap, il) in enumerate(zip(S.KI, S.caps, S.ill)):
        lin = il["lin"]
        rb = born_rytov_ratio(Oc, cfg, grid, ki, cap)
        lgB = np.log(rb)
        yC = msb_meas(msb_c, ki, il)
        gen = np.random.default_rng((SEED, j))               # exp06's noise draws
        for k6 in CH6:
            nz_ = gen.standard_normal((2, lin.size))
            if k6 not in CH:
                continue
            k, dvec = k6, CH[k6]
            h = il["sf"][k] * il["hb"]
            noise = np.sqrt(sig0 ** 2 / 2) * (nz_[0] + 1j * nz_[1])
            D[f"J_{k}"] += np.bincount(lin, weights=np.abs(h) ** 2 / sig0 ** 2, minlength=grid.size)
            y = {"A": h * Otc.reshape(-1)[lin],
                 "B": (S.sc * channel_meas(lgB, rb, grid, dvec))[il["ix"], il["iy"]] * il["phase"],
                 "C": yC[k], "L": (c / EPS) * yE[j][k], "N": noise}
            for m in "ABCLN":
                bin_add(D[f"num{m}_{k}"], lin, np.conj(h) * y[m] / sig0 ** 2)
            ray[f"yC_{k}"].append(y["C"].astype(np.complex64))
            ray[f"eC_{k}"].append((y["C"] - y["L"]).astype(np.complex64))
            ray[f"noise_{k}"].append(noise.astype(np.complex64))
        if (j + 1) % 60 == 0:
            print(f"    contrast {c}: illumination {j+1}/{N_ILL} [{time.time()-t0:.0f} s]", flush=True)
    D.update({key: np.concatenate(v) for key, v in ray.items()})
    D["sig0"] = sig0
    return D, yE


def load_data(S, c, cache={}):
    path = out_path(f"exp12_data_c{c}.npz")
    if not os.path.exists(path):
        print(f"  generating exp11's model-C data at contrast {c} (cached in {path})", flush=True)
        D, cache["yE"] = make_data(S, c, cache.get("yE"))
        np.savez(path, **D)
    D = dict(np.load(path))
    D["sig0"] = float(D["sig0"])
    return D


# ---------------------------------------------------------------- the estimator
class Recon:
    """exp06's nearest-voxel estimator on one contrast's data, with an optional
    per-channel correction subtracted from the back-projected data."""

    def __init__(self, S, c, D):
        self.S, self.c, self.D = S, c, D
        Otc = np.fft.fftn(S.potential(c))
        self.Ot = Otc
        self.Sc = shell_spectrum(Otc, S.Kabs, S.cfg.kmax, S.grid.dk)[0]
        self.Ob = np.real(np.fft.ifftn(Otc * S.meas))
        self.nOb = np.linalg.norm(self.Ob)

    def est(self, ks, how, model="C", corr=None):
        g, t = self.S.grid, 1.0 / len(ks)
        chans = []
        for k in ks:
            num = t * self.D[f"num{model}_{k}"] + np.sqrt(t) * self.D[f"numN_{k}"]
            if corr is not None:
                num = num - t * corr[k]
            chans.append((hsym_c(num.reshape(g.shape)),
                          hermitian_sym((t * self.D[f"J_{k}"]).reshape(g.shape))))
        return fuse(chans, self.Sc, how)[0]

    def rmse(self, est):
        Oe = np.real(np.fft.ifftn(est * self.S.meas))
        return np.linalg.norm(Oe - self.Ob) / self.nOb

    def three(self, model="C", corr=None):
        """RMSE of best design, equal-weight average and single 0.8 um shear
        on the same (corrected) data."""
        return {p: self.rmse(self.est(ks, how, model, corr)) for p, ks, how in PIPES}


def margins(r):
    return r["single"] / r["best"], r["avg"] / r["best"]


def fmt(r):
    s, a = margins(r)
    return (f"best {r['best']:.4f}  avg {r['avg']:.4f}  single {r['single']:.4f} | "
            f"vs single {s:5.2f}x  fusion vs avg {a:5.2f}x")


def excess_err(e, eC, ks=BEST):
    """||e - e_data|| / ||e_data|| over all rays of the given channels."""
    a = sum(np.sum(np.abs(e[k] - eC[k]) ** 2) for k in ks)
    b = sum(np.sum(np.abs(eC[k]) ** 2) for k in ks)
    return float(np.sqrt(a / b))


def iterate(S, R, ks, how, n_iter, mode="ms", bins=1, alpha=1.0, shared=False, log=None):
    """The fixed-point loop for one pipeline. Returns the per-iteration list of
    dicts (rmse of the pipeline, excess error at ray level, update size) and
    the final estimate (real space, full support) with its correction."""
    g, sig0 = S.grid, R.D["sig0"]
    eC = {k: R.D[f"eC_{k}"] for k in CH}
    corr = None
    est = R.est(ks, how)
    hist = [{"k": 0, "rmse": R.rmse(est), "exc": 1.0, "step": np.nan}]
    if shared:
        hist[0]["three"] = R.three()
    if log:
        log(hist[-1])
    for it in range(1, n_iter + 1):
        Ok = np.real(np.fft.ifftn(est))
        numF, numL, er = forward(S, Ok, sig0, bins=bins, rays=True)
        if mode == "ms":
            new = {k: numF[k] - numL[k] for k in CH}
        else:                                   # 'full': L = the on-voxel model W inverts
            Okt = np.fft.fftn(Ok).reshape(-1)
            new = {k: numF[k] - R.D[f"J_{k}"] * Okt for k in CH}
        corr = new if corr is None else {k: (1 - alpha) * corr[k] + alpha * new[k] for k in CH}
        est_new = R.est(ks, how, corr=corr)
        step = np.linalg.norm((est_new - est)[S.meas]) / np.linalg.norm(R.Ot[S.meas])
        est = est_new
        hist.append({"k": it, "rmse": R.rmse(est), "exc": excess_err(er, eC, ks), "step": step})
        if shared:
            hist[-1]["three"] = R.three(corr=corr)
        if log:
            log(hist[-1])
    return hist, np.real(np.fft.ifftn(est)), corr


# ------------------------------------------------------------------ self-check
def selfcheck(S):
    """BlochMSB against shearlib.MSB (on-lattice k_i) and against first Born."""
    g, cfg = S.grid, S.cfg
    O = S.potential(0.5)
    ki_on = g.dk * np.array([29.0, -6.0])                    # on the 12.8 um lattice
    a = BlochMSB(O, g.dx, S.zc, g.dz, cfg.k_m).ratio(ki_on)
    b = MSB(O, g.dx, g.dz, cfg.k_m).ratio(ki_on)
    e1 = np.abs(a - b).max()
    ki, cap, il = S.KI[7], S.caps[7], S.ill[7]               # off that lattice
    r = BlochMSB(EPS * 1e-2 * O, g.dx, S.zc, g.dz, cfg.k_m).ratio(ki)
    u1 = np.fft.fft2(np.log(r))[il["ix"], il["iy"]]
    u2 = born_UR(np.fft.fft2(EPS * 1e-2 * O, axes=(0, 1)), S.zc, g.dz, cap, il["ix"], il["iy"])
    e2 = np.linalg.norm(u1 - u2) / np.linalg.norm(u2)
    u3 = np.fft.fft2(np.log(born_rytov_ratio(EPS * 1e-2 * O, cfg, g, ki, cap)))[il["ix"], il["iy"]]
    e3 = np.linalg.norm(u3 - u2) / np.linalg.norm(u2)
    print(f"  BlochMSB = shearlib.MSB for k_i on the lattice: max |diff| {e1:.1e}")
    print(f"  BlochMSB's linear limit = continuum first Born (k_i off the lattice, contrast "
          f"{EPS*1e-2:.0e}): {e2:.1e};  born_UR = born_rytov_ratio: {e3:.1e}")
    assert e1 < 1e-9 and e2 < 2e-3 and e3 < 1e-9


# ------------------------------------------------------------------------ main
def main():
    t_all = time.time()
    S = Setup()
    print(RULE)
    print("exp12: ITERATIVE MULTIPLE-SCATTERING CORRECTION AROUND THE LINEAR ESTIMATOR")
    print(RULE)
    print(f"  {S.grid!r}; {N_ILL} illuminations at {S.th.mean():.2f} +- {S.th.std():.2f} deg; "
          f"{S.nray} rays; rho {RHO:.0f}")
    selfcheck(S)
    for c in CONTRASTS:
        run_contrast(S, c)
    print(f"  total {time.time()-t_all:.0f} s")


def run_contrast(S, c):
    D = load_data(S, c)
    R = Recon(S, c, D)
    sig0 = D["sig0"]
    eC = {k: D[f"eC_{k}"] for k in CH}
    print()
    print(RULE)
    print(f"CONTRAST {c:.2f}   (a dn/lambda = {S.adn(c):.3f})   RMSE vs band-limited truth, "
          f"nearest-voxel estimator")
    print(RULE)
    print("  1. Baselines (exp11's columns) and the ceiling")
    base = {m: R.three(m) for m in "ABCL"}
    for m, lab in (("A", "A ideal linear"), ("B", "B first Born, continuum"),
                   ("C", "C MSB data, uncorrected"), ("L", "crime oracle: C - own excess")):
        print(f"    {lab:30s} {fmt(base[m])}")
    for p, vals in EXP11.get(c, {}).items():
        for m, v in zip("ABC", vals[:3]):
            assert abs(base[m][p] - v) < 6e-5, f"exp11 {m} {p} not reproduced: {base[m][p]} vs {v}"
    print("    (A, B and C reproduce exp11's printed values)")
    sig = np.sqrt(sum(np.sum(np.abs(D[f'yC_{k}']) ** 2) for k in BEST))
    exc = np.sqrt(sum(np.sum(np.abs(eC[k]) ** 2) for k in BEST))
    print(f"    the data's nonlinear excess is {exc/sig:.3f} of the signal (all rays, four channels)")

    print("  2. Oracle: subtract the correction model's excess of the TRUE object (no iteration)")
    Otrue = S.potential(c)
    for bins in (1, 2):
        t0 = time.time()
        nF, nL, er = forward(S, Otrue, sig0, bins=bins, rays=True)
        corr = {k: nF[k] - nL[k] for k in CH}
        print(f"    BlochMSB 12.8 um periodic, dz {0.1*bins:.1f}: {fmt(R.three(corr=corr))}"
              f" | excess error {excess_err(er, eC):.3f}  [{time.time()-t0:.0f} s]", flush=True)

    print("  3. The loop, each pipeline correcting with ITS OWN estimate")
    res = {}
    for mode, alpha in (("ms", 1.0), ("full", 1.0)):
        print(f"    mode '{mode}', damping alpha = {alpha}:   [RMSE | ray-level excess error | step]")
        for p, ks, how in PIPES:
            t0 = time.time()
            hist, Ok, corr = iterate(S, R, ks, how, N_ITER, mode=mode, alpha=alpha,
                                     shared=(p == "best"))
            res[(mode, p)] = hist
            print(f"      {p:7s}" + " ".join(f"{h['rmse']:.4f}" for h in hist)
                  + f"  [{time.time()-t0:.0f} s]", flush=True)
            print(f"      {'':7s}" + " ".join(f"{h['exc']:6.3f}" for h in hist) + "  excess error")
            if p == "best":
                np.save(out_path(f"exp12_final_{mode}_c{c}.npy"), Ok.astype(np.float32))
        print(f"      margins per iteration (own estimates):")
        for i in range(N_ITER + 1):
            r = {p: res[(mode, p)][i]["rmse"] for p, _, _ in PIPES}
            print(f"        k={i}: {fmt(r)}")
        print(f"      margins per iteration (all three on the BEST design's corrected data):")
        for h in res[(mode, "best")]:
            print(f"        k={h['k']}: {fmt(h['three'])}")


if __name__ == "__main__":
    main()
