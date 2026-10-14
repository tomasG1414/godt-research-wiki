"""exp05b -- audit of exp05: precision-weighted averaging and the low-K shell.

Supports: wiki/notes/multishear-fusion-design-solution.md section 5 ([audit] rows)
and wiki/notes/multishear-fusion-derivation.md section 6.2-6.3.

Same data as exp05 (same seed, phantom, noise). For each design reports the
fused Wiener RMSE, the equal-weight ("naive") average, a per-voxel
precision-weighted average of the per-channel Wiener estimates (weights
J_j + 1/S), and the inverse-variance average of per-channel BLUEs (weights J_j),
which is algebraically the fused BLUE. Also reports the Wiener error inside the
|K| < K_LO shell that exp05 masks out.
"""
import time
import numpy as np
from shearlib import (GODT_SIM, KGrid, hermitian_sym, cap_density, azimuths, cell_phantom,
                      shell_spectrum, fuse)
import exp05_fusion_validation as e5

rng = np.random.default_rng(e5.SEED)
cfg = GODT_SIM
grid = KGrid(128, 0.1, 128)
Kabs = grid.Kabs()
phis = azimuths(e5.N_ILL)
M = hermitian_sym(cap_density(cfg, grid, e5.N_ILL, e5.N_FRAMES))
supp = M > 0
O, dn = cell_phantom(grid, rng)
Ot = np.fft.fftn(O)
S, edges, prof = shell_spectrum(Ot, Kabs, cfg.kmax, grid.dk)
ref = supp & (Kabs >= 8) & (Kabs <= 24)
sigma0 = np.sqrt(np.mean(np.abs(Ot[ref]) ** 2) / e5.RHO)
meas = supp & (Kabs >= e5.K_LO)
low = supp & (Kabs < e5.K_LO)
O_band = np.real(np.fft.ifftn(Ot * meas))
print("low shell voxels", low.sum(), " DC in supp:", bool(supp[0, 0, 0]),
      " min K_perp in low shell:", grid.KP[:, :, None].repeat(grid.nz, 2)[low].min())
print("low-shell norm frac 0.842 check:",
      np.linalg.norm(np.real(np.fft.ifftn(Ot * low))) / np.linalg.norm(O))
# fraction of that norm that is the K=0..first shell
X, Y = (1.0, 0.0), (0.0, 1.0)
def ch(m, u, t): return ((m * u[0], m * u[1]), t)
designs = [
    ("ODT", [(None, 1.0)]),
    ("0.1 x", [ch(0.1, X, 1.0)]),
    ("0.8 x", [ch(0.8, X, 1.0)]),
    ("0.8x+0.8y", [ch(0.8, X, .5), ch(0.8, Y, .5)]),
    ("{0.8,0.55}x+y", [ch(0.8, X, .25), ch(0.55, X, .25), ch(0.8, Y, .25), ch(0.55, Y, .25)]),
]
def rmse(est):
    Oe = np.real(np.fft.ifftn(est * meas))
    return np.linalg.norm(Oe - O_band) / np.linalg.norm(O_band)
for name, chans in designs:
    t0 = time.time()
    acq = e5.run_design(cfg, grid, Ot, chans, sigma0, phis, rng)
    w, J = fuse(acq, S, "wiener")
    # low-K shell: relative (power-weighted) Wiener MSE inside supp & |K|<2
    lowmse = (np.abs(w - Ot) ** 2)[low].sum() / (np.abs(Ot) ** 2)[low].sum()
    pred_low = (S / (S * J + 1))[low].sum() / (np.abs(Ot) ** 2)[low].sum()
    line = f"{name:16s} wiener {rmse(w):.4f}  lowK-shell relMSE emp {lowmse:.3f} pred {pred_low:.3f}"
    if len(chans) > 1:
        naive, _ = fuse(acq, S, "naive")
        # precision-weighted average of per-channel Wiener estimates:
        # w_j = posterior precision J_j + 1/S
        ests = [c[0] / (c[1] + 1 / S) for c in acq]
        prec = [c[1] + 1 / S for c in acq]
        pwa = sum(p * e for p, e in zip(prec, ests)) / sum(prec)
        # inverse-variance average of per-channel BLUE estimates (weights J_j)
        with np.errstate(divide="ignore", invalid="ignore"):
            bl = [np.where(c[1] > 0, c[0] / np.where(c[1] > 0, c[1], 1), 0) for c in acq]
        Jt = sum(c[1] for c in acq)
        ivb = np.where(Jt > 0, sum(c[1] * b for c, b in zip(acq, bl)) / np.where(Jt > 0, Jt, 1), 0)
        blue, _ = fuse(acq, S, "blue")
        line += (f" | naive {rmse(naive):.4f}  prec-weighted avg {rmse(pwa):.4f}"
                 f"  inv-var avg of BLUEs {rmse(ivb):.4f}  (fused BLUE {rmse(blue):.4f})")
    print(line, f"[{time.time()-t0:.0f}s]", flush=True)
