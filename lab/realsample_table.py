"""realsample_table -- real samples on the a*dn/lambda axis.

Supports wiki/notes/real-samples-on-the-scattering-axis.md. Prints the computed
column of that note's table from the inputs listed below, so every number in
the note can be regenerated. Arithmetic only; stdlib. Each input carries the
source it was read from (see the note for how each was read) or is marked
ASSUMED. Nothing here is a measurement of ours.

    python realsample_table.py            # the table
    python realsample_table.py --phantom  # also: mean contrast of the project's
                                          # exp06 phantom (needs numpy, scipy)

x = a * dn / lambda.  a: radius or half path length (um).  dn: index contrast
to the surrounding medium.  Dispersion between the wavelength of each
measurement and 0.50 / 0.53 um is ignored (stated in the note).
"""
import math
import sys

LAMS = (0.53, 0.50)          # GODT experiment (LED centre), GODT simulation


def r_sphere(v_fl):
    """Radius (um) of the sphere with volume v_fl (fL = um^3)."""
    return (3.0 * v_fl / (4.0 * math.pi)) ** (1.0 / 3.0)


# Medium indices used when a paper does not state its own (ASSUMED range):
# PBS 1.3332 and buffer medium 1.3379 (Baczewska 2021), medium 1.337 at
# 532 nm (Kim 2016).
N_MED = (1.333, 1.338)

# RBC: mean thickness from volume 93.1 fL over a disc of projected diameter 8 um
RBC_HALF_THICKNESS = 93.1 / (math.pi * 4.0 ** 2) / 2.0
# E. coli: absolute 1.384, relative 1.0387 (Balaev 2002) -> contrast to water
ECOLI_DN = 1.384 - 1.384 / 1.0387

# (label, (a_lo, a_hi), (dn_lo, dn_hi), kind)   kind: max | mean | uniform
ROWS = [
    # --- project and GODT simulations --------------------------------------
    ("exp06 phantom, contrast 1, dn_max",          (2.976, 2.976), (0.0769, 0.0769), "max"),
    ("exp06 phantom, contrast 1, cytoplasm step",  (2.976, 2.976), (0.0200, 0.0200), "mean"),
    ("GODT Shepp-Logan (n 1.33-1.335, 40 um)",     (20.0, 20.0),   (0.005, 0.005),   "max"),
    # --- GODT's experimental samples ---------------------------------------
    # Ziemczonok 2019: 30 x 25 x 11 um; polymer 1.5167 / 1.5021; oil 1.4878
    ("printed phantom, half height, high dose",    (5.5, 5.5),     (1.5167 - 1.4878,) * 2, "max"),
    ("printed phantom, half height, low dose",     (5.5, 5.5),     (1.5021 - 1.4878,) * 2, "mean"),
    ("printed phantom, lateral semi-axes",         (12.5, 15.0),   (1.5021 - 1.4878, 1.5167 - 1.4878), "mean"),
    ("printed phantom nucleolus vs nucleus",       (1.3, 1.3),     (1.5167 - 1.5021,) * 2, "uniform"),
    # Lee 2021: soma RI 1.3424-1.3990; a ASSUMED 5-10 um; medium ASSUMED
    ("neuron soma, max contrast",                  (5.0, 10.0),    (1.3990 - N_MED[1], 1.3990 - N_MED[0]), "max"),
    ("neuron soma, mean contrast ASSUMED",         (5.0, 10.0),    (0.010, 0.030),   "mean"),
    # --- other cells --------------------------------------------------------
    # Choi 2007: overall 1.36-1.40; Kim 2019: nucleolus 1.363, cytoplasm 1.349
    ("HeLa, max contrast",                         (3.0, 6.0),     (1.363 - 1.337, 1.40 - 1.337), "max"),
    ("HeLa, mean contrast",                        (3.0, 6.0),     (1.349 - 1.337, 1.375 - 1.337), "mean"),
    # Park 2008: 1.399; PBS 1.3332 (Baczewska 2021)
    ("red blood cell, half mean thickness",        (RBC_HALF_THICKNESS,) * 2, (1.399 - 1.3332,) * 2, "uniform"),
    ("red blood cell, lateral radius",             (4.0, 4.0),     (1.399 - 1.3332,) * 2, "uniform"),
    # Midtvedt 2019: n 1.402-1.420 (IQR), V 17-37 fL (IQR)
    ("yeast",                                      (r_sphere(17), r_sphere(37)), (1.402 - N_MED[1], 1.420 - N_MED[0]), "uniform"),
    # Oh 2020: V 0.5-1.4 fL; Balaev 2002 contrast
    ("E. coli, equal-volume sphere",               (r_sphere(0.5), r_sphere(1.4)), (ECOLI_DN,) * 2, "uniform"),
    # --- small, high-contrast structures (contrast to the cytoplasm) -------
    # Kim 2016: LD 1.386, cytoplasm 1.347, V 0.34-0.75 fL
    ("lipid droplet (Kim 2016)",                   (r_sphere(0.34), r_sphere(0.75)), (1.386 - 1.347,) * 2, "uniform"),
    # CP-HHT (wiki source page): LD 1.468-1.479; cytoplasm 1.349 (Kim 2019)
    ("lipid droplet (CP-HHT RI, Kim 2016 size)",   (r_sphere(0.34), r_sphere(0.75)), (1.468 - 1.349, 1.479 - 1.349), "uniform"),
    # Kim 2019: 1.363 vs 1.349; Choi 2007: up to 1.385 vs 1.355; a ASSUMED 1-2 um
    ("nucleolus vs nucleoplasm",                   (1.0, 2.0),     (0.010, 0.030),   "uniform"),
    # --- thick samples ------------------------------------------------------
    # Yasuhiko & Takeuchi 2023: diameter 140 um; mean contrast ASSUMED
    ("spheroid, 140 um, contrast ASSUMED",         (70.0, 70.0),   (0.010, 0.030),   "mean"),
    # Sun 2012 (1310 nm): grey matter 1.361-1.369, aCSF 1.342, slice 300 um
    ("acute brain slice, 300 um, grey matter",     (150.0, 150.0), (1.361 - 1.342, 1.369 - 1.342), "mean"),
]


def table():
    print(f"{'sample':46s} {'a (um)':>13s} {'dn':>15s} {'kind':>8s} "
          + " ".join(f"{'x @' + format(l, '.2f'):>13s}" for l in LAMS))
    for label, (a0, a1), (d0, d1), kind in ROWS:
        a = f"{a0:.2f}" if a0 == a1 else f"{a0:.2f}-{a1:.2f}"
        d = f"{d0:.4f}" if d0 == d1 else f"{d0:.4f}-{d1:.4f}"
        xs = []
        for lam in LAMS:
            lo, hi = a0 * d0 / lam, a1 * d1 / lam
            xs.append(f"{lo:.3f}" if abs(hi - lo) < 5e-4 else f"{lo:.3f}-{hi:.3f}")
        print(f"{label:46s} {a:>13s} {d:>15s} {kind:>8s} " + " ".join(f"{x:>13s}" for x in xs))


def phantom():
    """Mean and maximum contrast of the exp06 phantom, as exp06 builds it."""
    import numpy as np
    import shearlib as sl
    from exp06_multiple_scattering import PHANTOM_SCALE, SEED
    grid = sl.KGrid(128, 0.1, 128)
    _, dn = sl.cell_phantom(grid, np.random.default_rng(SEED), scale=PHANTOM_SCALE)
    n, nz = grid.n, grid.nz
    x = (np.arange(n) - n // 2) * grid.dx
    z = (np.arange(nz) - nz // 2) * grid.dz
    X, Y, Z = np.meshgrid(x, x, z, indexing="ij")
    a = 4.8 * PHANTOM_SCALE
    cyto = X ** 2 + Y ** 2 + Z ** 2 <= a ** 2
    dn_max, dn_mean = float(np.abs(dn).max()), float(dn[cyto].mean())
    print(f"\nexp06 phantom (scale {PHANTOM_SCALE}, contrast 1): a = {a:.3f} um, "
          f"dn_max = {dn_max:.4f}, volume-mean dn = {dn_mean:.4f}, ratio = {dn_mean / dn_max:.3f}")
    for c in (0.25, 0.5, 1.0, 2.0):
        print(f"  contrast {c:4.2f}: a dn_max/lam = {a * c * dn_max / 0.5:.3f}   "
              f"a dn_mean/lam = {a * c * dn_mean / 0.5:.3f}")


if __name__ == "__main__":
    table()
    if "--phantom" in sys.argv:
        phantom()
