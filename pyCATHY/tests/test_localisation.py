"""Localisation: Gaspari-Cohn conventions (B24), the matrix, and its effect on the analysis (E6)."""
import numpy as np
import pytest

from pyCATHY.DA import enkf, localisation


# ---- Gaspari-Cohn ---------------------------------------------------------------------
def test_weight_is_one_at_zero_and_vanishes_beyond_two_L():
    L = 3.0
    w = localisation.gaspari_cohn(np.array([0.0, 2 * L, 2 * L + 1e-9, 10 * L]), L)
    assert w[0] == pytest.approx(1.0)
    assert w[1] == pytest.approx(0.0, abs=1e-12)  # support is 2 L, not L
    assert w[2] == 0.0 and w[3] == 0.0


def test_weight_is_monotone_and_continuous_at_L():
    L = 2.0
    r = np.linspace(0, 2 * L, 2001)
    w = localisation.gaspari_cohn(r, L)
    assert np.all(np.diff(w) <= 1e-12)
    assert 0 < localisation.gaspari_cohn(np.array([L]), L)[0] < 1
    below, above = localisation.gaspari_cohn(np.array([L * (1 - 1e-9), L * (1 + 1e-9)]), L)
    assert abs(below - above) < 1e-6  # the two polynomial branches meet


def test_weight_depends_on_r_over_L_only():
    r = np.array([0.0, 0.7, 1.9, 3.1])
    np.testing.assert_allclose(localisation.gaspari_cohn(r, 1.0),
                               localisation.gaspari_cohn(10 * r, 10.0))


@pytest.mark.parametrize("bad", [0, -1.0])
def test_non_positive_length_scale_is_an_error(bad):
    with pytest.raises(ValueError, match="must be > 0"):
        localisation.gaspari_cohn(np.array([1.0]), bad)


# ---- matrix ---------------------------------------------------------------------------
def _loop_reference(obs, grid, L):
    """The former implementation (one Python iteration per grid point)."""
    out = np.zeros((grid.shape[0], obs.shape[0]))
    for i in range(grid.shape[0]):
        d = np.sqrt(((obs - grid[i]) ** 2).sum(axis=1))
        out[i] = localisation.gaspari_cohn(d, L)
    return out


def test_matrix_matches_the_loop_implementation_and_has_the_documented_shape():
    rng = np.random.default_rng(0)
    obs, grid = rng.uniform(0, 20, (30, 2)), rng.uniform(0, 20, (200, 2))
    loc = localisation.build_localization_matrix(obs, grid, 4.0)
    assert loc.shape == (200, 30)  # (n_grid, n_obs)
    np.testing.assert_allclose(loc, _loop_reference(obs, grid, 4.0), atol=1e-15)


def test_column_order_does_not_matter_if_obs_and_grid_agree():
    rng = np.random.default_rng(1)
    obs, grid = rng.uniform(0, 10, (8, 2)), rng.uniform(0, 10, (40, 2))
    a = localisation.build_localization_matrix(obs, grid, 3.0)
    b = localisation.build_localization_matrix(obs[:, ::-1], grid[:, ::-1], 3.0)  # [y, x] vs [x, y]
    np.testing.assert_allclose(a, b)


def test_depth_column_localises_in_depth():
    """Passing z makes the weight depth dependent (the surface obs does not see deep nodes)."""
    obs = np.array([[0.0, 0.0, 0.0]])
    grid = np.array([[0.0, 0.0, -0.1], [0.0, 0.0, -10.0]])
    loc = localisation.build_localization_matrix(obs, grid, 1.0)[:, 0]
    assert loc[0] > 0.9 and loc[1] == 0.0


def test_mismatched_coordinate_columns_raise():
    with pytest.raises(ValueError, match="coordinate columns"):
        localisation.build_localization_matrix(np.zeros((2, 2)), np.zeros((3, 3)), 1.0)


# ---- effect on the analysis (E6) -------------------------------------------------------
def _problem(seed=0, n=12, m=20, nens=15):
    rng = np.random.default_rng(seed)
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens = rng.normal(size=(n, nens))
    param = 0.4 + 0.05 * rng.normal(size=(1, nens))
    y = H @ rng.normal(size=n) + 0.3 * rng.normal(size=m)
    return y, np.eye(m) * 0.09, ens, param, H @ ens


def _localised(L, **kw):
    y, R, ens, param, pred = _problem()
    out = enkf.enkf_analysis_localized_with_inflation(
        y, R, ens, param, pred, L, Sakov=False, inflate_states=1.0, inflate_params=1.0,
        jitter_params=0.0, **kw)
    return ens, param, out[9], out[10]


def test_weights_of_one_give_the_unlocalised_analysis():
    y, R, ens, param, pred = _problem()
    plain = enkf.enkf_analysis(y, R, param, ens, pred, Sakov=False)
    _, _, ana, ana_p = _localised(np.ones((12, 20)))
    np.testing.assert_allclose(ana, plain[-2], atol=1e-10)


def test_weights_of_zero_leave_the_states_untouched_but_update_the_parameters():
    """Localisation is applied to the states only (documented), not to the parameters."""
    ens, param, ana, ana_p = _localised(np.zeros((12, 20)))
    np.testing.assert_allclose(ana, ens)
    assert not np.allclose(ana_p, param)


def test_wrong_localisation_shape_is_rejected():
    with pytest.raises(ValueError, match="Localization matrix shape"):
        _localised(np.ones((5, 20)))


def test_localised_filter_supports_perturbed_observations():
    np.random.seed(2)
    _, _, a0, _ = _localised(np.ones((12, 20)))
    np.random.seed(2)
    _, _, a1, _ = _localised(np.ones((12, 20)), perturb_obs=True)
    assert not np.allclose(a0, a1)
    np.testing.assert_allclose(a0.mean(1), a1.mean(1), atol=0.5)


def test_weights_stay_in_zero_one_even_at_the_edge_of_the_support():
    """The polynomial gave -2.8e-16 at r = 2 L (rounding)."""
    r = np.concatenate([np.linspace(0, 3, 3001), [2.0, 1.0, 1.9999999999]])
    w = localisation.gaspari_cohn(r, 1.0)
    assert w.min() >= 0.0 and w.max() <= 1.0
    assert localisation.gaspari_cohn(np.array([2.0]), 1.0)[0] == 0.0
