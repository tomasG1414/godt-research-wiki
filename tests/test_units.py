"""Unit tests for shearlib: one property per test, so a failure names what broke.

The heavier invariants (cap density vs closed form, WPM, design solver) stay in
shearlib.selftest(), run by test_shearlib.py.
"""
import numpy as np
import pytest

import shearlib as sl


@pytest.mark.parametrize("d", [0.1, 0.33, 0.8])
def test_absH_is_the_modulus_of_H(d):
    K = np.linspace(-40, 40, 4001)
    np.testing.assert_allclose(np.abs(sl.H(K, d)), sl.absH(K, d), rtol=0, atol=1e-12)


@pytest.mark.parametrize("d, expected", [
    (0.1, []),                      # GODT simulation shear: no in-band null
    (0.33, [19.04]),                # preprint experiment
    (0.8, [7.85, 15.71, 23.56]),    # published experiment
])
def test_in_band_nulls_match_the_wiki_table(d, expected):
    # wiki/notes/finite-shear-null-experiment.md, section 2
    got = sl.nulls_in_band(d, sl.GODT_SIM.kmax)
    assert len(got) == len(expected)
    np.testing.assert_allclose(got, expected, atol=5e-3)


def test_nulls_in_band_rejects_a_nonpositive_shear():
    assert sl.nulls_in_band(0.0, 30.0) == []
    assert sl.nulls_in_band(-0.5, 30.0) == []


@pytest.mark.parametrize("rho", [1.0, 1e3, 1e9])
def test_a_null_is_lost_at_any_snr(rho):
    assert sl.wiener_rel_err(2 * np.pi / 0.8, [0.8], rho) == pytest.approx(1.0)


def test_wiener_fusion_is_the_precision_weighted_average_of_blues():
    # multishear-fusion-derivation.md section 6.2, the identity behind the "21x" claim
    rng = np.random.default_rng(0)
    shape = (6, 5, 4)
    chans = []
    for _ in range(3):
        J = rng.uniform(0.1, 5.0, shape)
        chans.append((J * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)), J))
    S = rng.uniform(0.5, 2.0, shape)
    fused, J = sl.fuse(chans, S, "wiener")
    weighted = sum(Jj * (num / Jj) for num, Jj in chans) / (J + 1.0 / S)
    np.testing.assert_allclose(fused, weighted, rtol=1e-12)


def test_blue_is_exact_without_noise_and_zero_where_nothing_samples():
    rng = np.random.default_rng(1)
    shape = (8, 8, 4)
    O = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    J = rng.uniform(0.1, 3.0, shape)
    J[0, 0, :] = 0.0                                  # a voxel column no channel sees
    est, _ = sl.fuse([(J * O, J)], np.ones(shape), "blue")
    np.testing.assert_allclose(est[J > 0], O[J > 0], rtol=1e-12)
    assert np.all(est[J == 0] == 0)


def test_equal_weight_average_dilutes_a_frequency_one_channel_cannot_see():
    shape = (4, 4, 4)
    O = np.ones(shape, complex)
    S = np.full(shape, 1e9)                           # flat prior: Wiener ~ BLUE
    seen, blind = np.full(shape, 10.0), np.zeros(shape)
    chans = [(seen * O, seen), (blind * O, blind)]
    np.testing.assert_allclose(sl.fuse(chans, S, "wiener")[0], O, rtol=1e-6)
    np.testing.assert_allclose(sl.fuse(chans, S, "naive")[0], 0.5 * O, rtol=1e-6)


def test_hermitian_sym_does_not_double_count_self_partner_voxels():
    S = sl.hermitian_sym(np.ones((4, 4, 4)))
    assert S[0, 0, 0] == 1 and S[2, 2, 2] == 1 and S[2, 0, 0] == 1   # origin / Nyquist
    assert S[1, 0, 0] == 2 and S[1, 2, 3] == 2                        # K and -K both count


def test_cap_samples_lie_inside_the_pupil_and_the_grid():
    cfg, grid, phi = sl.GODT_SIM, sl.KGrid(64, 0.1, 24), 0.3
    ix, iy, iz, kz = sl.cap_samples(cfg, grid, phi)
    kxp = grid.KX[ix, iy] + cfg.k_ill * np.cos(phi)
    kyp = grid.KY[ix, iy] + cfg.k_ill * np.sin(phi)
    assert ix.size > 0
    assert np.all(np.hypot(kxp, kyp) <= cfg.k_perp_max + 1e-9)
    assert np.all((0 <= iz) & (iz < grid.nz))
    assert np.all(kz >= cfg.kz_floor - 1e-9)
