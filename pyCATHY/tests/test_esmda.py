"""Unit tests of the ES-MDA numerics (no CATHY run needed).

Reference: for a linear Gaussian problem ES-MDA with sum(1/alpha) = 1 must converge
(large ensemble) to the analytical Kalman posterior.
"""
import logging

import numpy as np
import pytest

from pyCATHY.DA.esmda import ESMDAConfig, compute_alphas, esmda_update


# --------------------------------------------------------------------------- alphas
@pytest.mark.parametrize("schedule", ["constant", "geometric"])
@pytest.mark.parametrize("n_iter", [1, 2, 4, 8])
def test_alphas_satisfy_constraint(schedule, n_iter):
    alphas = compute_alphas(n_iter, schedule)
    assert alphas.shape == (n_iter,)
    assert np.isclose(np.sum(1.0 / alphas), 1.0)
    assert np.all(alphas > 0)


def test_constant_alphas_equal_n_iter():
    assert np.allclose(compute_alphas(5), 5.0)


def test_geometric_alphas_decrease():
    alphas = compute_alphas(5, "geometric", geometric_ratio=0.6)
    assert np.all(np.diff(alphas) < 0)


def test_user_alphas_valid_and_invalid():
    ok = [2.0, 4.0, 4.0]  # 1/2 + 1/4 + 1/4 = 1
    assert np.allclose(compute_alphas(3, "user", ok), ok)
    with pytest.raises(ValueError, match="sum"):
        compute_alphas(3, "user", [2.0, 2.0, 2.0])
    with pytest.raises(ValueError, match="entries"):
        compute_alphas(3, "user", [1.0, 2.0])
    with pytest.raises(ValueError, match="requires"):
        compute_alphas(3, "user")
    with pytest.raises(ValueError, match="positive"):
        compute_alphas(2, "user", [-2.0, 2.0])


@pytest.mark.parametrize("bad", [0, -1, 2.5])
def test_invalid_n_iter(bad):
    with pytest.raises(ValueError, match="n_iter"):
        compute_alphas(bad)


def test_invalid_schedule_and_ratio():
    with pytest.raises(ValueError, match="schedule"):
        compute_alphas(3, "bogus")
    with pytest.raises(ValueError, match="geometric_ratio"):
        compute_alphas(3, "geometric", geometric_ratio=1.5)


# --------------------------------------------------------------------------- update
def _linear_problem(n_ens, n_param=3, n_obs=5, seed=0):
    rng = np.random.default_rng(seed)
    G = rng.normal(size=(n_obs, n_param))
    A = rng.normal(size=(n_param, n_param))
    C0 = A @ A.T / n_param + 0.5 * np.eye(n_param)
    m0 = rng.normal(size=n_param)
    m_true = rng.multivariate_normal(m0, C0)
    R = np.diag(rng.uniform(0.05, 0.2, size=n_obs))
    d = G @ m_true + rng.multivariate_normal(np.zeros(n_obs), R)
    ens = rng.multivariate_normal(m0, C0, size=n_ens).T
    return G, C0, m0, R, d, ens


def _kalman(G, C0, m0, R, d):
    K = C0 @ G.T @ np.linalg.inv(G @ C0 @ G.T + R)
    return m0 + K @ (d - G @ m0), C0 - K @ G @ C0


@pytest.mark.parametrize("schedule", ["constant", "geometric"])
def test_linear_gaussian_matches_kalman_posterior(schedule):
    G, C0, m0, R, d, ens = _linear_problem(n_ens=40000)
    mean_ref, cov_ref = _kalman(G, C0, m0, R, d)

    rng = np.random.default_rng(1)
    param = ens
    for alpha in compute_alphas(4, schedule):
        param = esmda_update(param, G @ param, d, R, alpha, rng=rng)

    assert np.allclose(param.mean(axis=1), mean_ref, atol=0.03)
    assert np.allclose(np.cov(param), cov_ref, atol=0.03)


def test_diagonal_vector_and_matrix_cov_equivalent():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=50)
    out_mat = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(3))
    out_vec = esmda_update(ens, G @ ens, d, np.diag(R), 2.0, rng=np.random.default_rng(3))
    assert np.allclose(out_mat, out_vec)


def test_full_covariance_runs():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=50)
    L = np.tril(np.ones_like(R)) * 0.01
    R_full = R + L @ L.T
    out = esmda_update(ens, G @ ens, d, R_full, 2.0, rng=np.random.default_rng(0))
    assert out.shape == ens.shape and np.all(np.isfinite(out))


def test_many_obs_uses_ensemble_space_and_agrees():
    # n_obs > n_ens triggers the Woodbury branch of enkf._solve_obs_space
    rng = np.random.default_rng(5)
    n_param, n_obs, n_ens = 2, 80, 10
    G = rng.normal(size=(n_obs, n_param))
    ens = rng.normal(size=(n_param, n_ens))
    d = rng.normal(size=n_obs)
    R = 0.1 * np.ones(n_obs)
    out = esmda_update(ens, G @ ens, d, R, 3.0, rng=np.random.default_rng(9))

    # reference with the full (m x m) solve
    rng2 = np.random.default_rng(9)
    pert = d[:, None] + np.sqrt(3.0) * np.sqrt(R)[:, None] * rng2.standard_normal((n_obs, n_ens))
    pred = G @ ens
    S = pred - pred.mean(axis=1, keepdims=True)
    dM = ens - ens.mean(axis=1, keepdims=True)
    full = S @ S.T / (n_ens - 1) + 3.0 * np.diag(R)
    ref = ens + (dM @ S.T / (n_ens - 1)) @ np.linalg.solve(full, pert - pred)
    assert np.allclose(out, ref)


def test_seed_reproducibility_and_fresh_draws():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    a = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(7))
    b = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(7))
    c = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(8))
    assert np.array_equal(a, b)
    assert not np.allclose(a, c)

    # successive calls on the same generator draw new perturbations
    rng = np.random.default_rng(7)
    first = esmda_update(ens, G @ ens, d, R, 2.0, rng=rng)
    second = esmda_update(ens, G @ ens, d, R, 2.0, rng=rng)
    assert not np.allclose(first, second)


def test_update_reduces_misfit():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=200)
    rng = np.random.default_rng(2)

    def misfit(p):
        return np.mean(np.sum((d[:, None] - G @ p) ** 2 / np.diag(R)[:, None], axis=0))

    before = misfit(ens)
    param = ens
    for alpha in compute_alphas(4):
        param = esmda_update(param, G @ param, d, R, alpha, rng=rng)
    assert misfit(param) < before


def test_localisation_zero_leaves_parameters_unchanged():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    zeros = np.zeros((ens.shape[0], d.size))
    out = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(0), localisation=zeros)
    assert np.allclose(out, ens)


def test_inflation_scales_anomalies():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    base = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(0))
    infl = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(0), inflation=1.5)
    assert np.allclose(infl.mean(axis=1), base.mean(axis=1))
    assert np.allclose(
        infl - infl.mean(axis=1, keepdims=True),
        1.5 * (base - base.mean(axis=1, keepdims=True)),
    )


def test_update_input_validation():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=10)
    pred = G @ ens
    with pytest.raises(ValueError, match="members"):
        esmda_update(ens, pred[:, :5], d, R, 2.0)
    with pytest.raises(ValueError, match="data has"):
        esmda_update(ens, pred, d[:-1], R, 2.0)
    with pytest.raises(ValueError, match="alpha"):
        esmda_update(ens, pred, d, R, 0.0)
    with pytest.raises(ValueError, match="data_cov"):
        esmda_update(ens, pred, d, np.ones(3), 2.0)
    with pytest.raises(ValueError, match="localisation"):
        esmda_update(ens, pred, d, R, 2.0, localisation=np.ones((1, 1)))
    with pytest.raises(ValueError, match="2 ensemble"):
        esmda_update(ens[:, :1], pred[:, :1], d, R, 2.0)


# --------------------------------------------------------------------------- config
def test_config_defaults_and_alphas():
    cfg = ESMDAConfig()
    assert cfg.n_iter == 4
    assert np.isclose(np.sum(1.0 / cfg.alphas), 1.0)
    assert cfg.update_initial_conditions is False


def test_config_validation():
    with pytest.raises(ValueError):
        ESMDAConfig(n_iter=0)
    with pytest.raises(ValueError):
        ESMDAConfig(alpha_schedule="user", user_alphas=[1.0, 1.0], n_iter=2)
    with pytest.raises(ValueError, match="inflation_params"):
        ESMDAConfig(inflation_params=0.0)


def test_config_empty_localisation_becomes_none():
    assert ESMDAConfig(localisation_matrix=[]).localisation_matrix is None
    assert ESMDAConfig(localisation_matrix=np.ones((2, 3))).localisation_matrix.shape == (2, 3)


# --------------------------------------------------------------------------- registry
def test_esmda_registered_and_run_analysis_contract():
    from pyCATHY.DA.analysis import ANALYSES, AnalysisConfig, run_analysis

    assert "esmda" in ANALYSES
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    state = np.random.default_rng(0).normal(size=(12, 30))
    sw = state.copy()
    cfg = AnalysisConfig(esmda_alpha=2.0, esmda_rng=np.random.default_rng(4))

    res = run_analysis(
        "esmda", d, R, ens, ["p1", "p2", "p3"], [state, sw], G @ ens, config=cfg
    )
    assert res.analysis is state  # state untouched
    assert res.analysis_param.shape == (30, 3)  # (n_members, n_param)
    assert res.extra["alpha"] == 2.0
    assert res.diagnostics is None
    # same numbers as the pure update with the same seed
    ref = esmda_update(ens, G @ ens, d, R, 2.0, rng=np.random.default_rng(4))
    assert np.allclose(res.analysis_param, ref.T)


def test_esmda_requires_parameters():
    from pyCATHY.DA.analysis import AnalysisConfig, run_analysis

    G, C0, m0, R, d, ens = _linear_problem(n_ens=10)
    state = np.zeros((4, 10))
    with pytest.raises(ValueError, match="parameters"):
        run_analysis(
            "esmda", d, R, np.zeros((0, 10)), [], [state, state], G @ ens,
            config=AnalysisConfig(),
        )


# --------------------------------------------------------------------------- driver
from pyCATHY.DA.esmda import ESMDAResult, run_esmda, weighted_misfit


def _forward_linear(G, fail=None, nan_rows=None):
    """Mock forward model g(m) = G m. ``fail[iteration]`` = members that fail then."""
    calls = {"n": 0}

    def forward(param, valid):
        n_ens = param.shape[1]
        pred = G @ param
        ok = np.ones(n_ens, dtype=bool)
        if fail and calls["n"] in fail:
            ok[list(fail[calls["n"]])] = False
        if nan_rows is not None:
            pred[nan_rows, :] = np.nan
        calls["n"] += 1
        return pred, ok

    forward.calls = calls
    return forward


def test_weighted_misfit_values():
    pred = np.array([[1.0, 3.0], [2.0, 2.0]])
    data = np.array([2.0, 2.0])
    # residuals: member0 = (1, 0), member1 = (-1, 0); var = 1 -> r^T R^-1 r = 1 each; /n_obs=2
    assert np.isclose(weighted_misfit(data, np.ones(2), pred), 0.5)
    assert np.isclose(weighted_misfit(data, np.eye(2), pred), 0.5)


def test_run_esmda_linear_matches_kalman_and_shapes():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=20000)
    mean_ref, cov_ref = _kalman(G, C0, m0, R, d)
    cfg = ESMDAConfig(n_iter=4, seed=3)
    res = run_esmda(ens, _forward_linear(G), d, R, cfg)

    assert isinstance(res, ESMDAResult)
    assert res.param.shape == ens.shape
    assert len(res.param_history) == 5  # prior + 4 updates
    assert len(res.prediction_history) == 5  # 4 runs + final
    assert res.misfit_history.shape == (5,)
    assert len(res.valid_history) == 4
    assert np.allclose(res.alphas, 4.0)
    assert np.allclose(res.param_history[0], ens)  # prior preserved, input not mutated
    assert np.allclose(res.param.mean(axis=1), mean_ref, atol=0.05)
    assert np.allclose(np.cov(res.param), cov_ref, atol=0.05)


def test_run_esmda_does_not_mutate_input_and_is_reproducible():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=40)
    before = ens.copy()
    a = run_esmda(ens, _forward_linear(G), d, R, ESMDAConfig(n_iter=3, seed=1))
    b = run_esmda(ens, _forward_linear(G), d, R, ESMDAConfig(n_iter=3, seed=1))
    c = run_esmda(ens, _forward_linear(G), d, R, ESMDAConfig(n_iter=3, seed=2))
    assert np.array_equal(ens, before)
    assert np.array_equal(a.param, b.param)
    assert not np.allclose(a.param, c.param)


def test_run_esmda_nonlinear_reduces_misfit_and_recovers_truth():
    rng = np.random.default_rng(0)
    n_ens, n_obs = 300, 6
    t = np.linspace(0.2, 2.0, n_obs)

    def g(m):  # nonlinear in the parameter, e.g. exponential decay amplitude/rate
        return m[0][None, :] * np.exp(-m[1][None, :] * t[:, None])

    m_true = np.array([2.0, 0.8])
    sigma = 0.02
    d = g(m_true[:, None])[:, 0] + sigma * rng.normal(size=n_obs)
    prior = np.vstack([rng.normal(1.5, 0.6, n_ens), rng.normal(1.0, 0.4, n_ens)])

    def forward(param, valid):
        return g(param), np.ones(param.shape[1], dtype=bool)

    res = run_esmda(prior, forward, d, np.full(n_obs, sigma**2), ESMDAConfig(n_iter=6, seed=0))
    assert res.misfit_history[-1] < res.misfit_history[0] / 5
    assert np.allclose(res.param.mean(axis=1), m_true, atol=0.1)
    assert np.all(res.param.std(axis=1) < prior.std(axis=1))


def test_run_esmda_callback_and_alphas():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=20)
    seen = []
    run_esmda(
        ens, _forward_linear(G), d, R, ESMDAConfig(n_iter=3),
        callback=lambda it, a, p, pr, v: seen.append((it, a, p.shape, pr.shape, int(v.sum()))),
    )
    assert [s[0] for s in seen] == [0, 1, 2, 3]
    assert [s[1] for s in seen] == [3.0, 3.0, 3.0, None]  # final run has no alpha
    assert all(s[2] == ens.shape and s[3] == (d.size, 20) and s[4] == 20 for s in seen)


def test_run_esmda_final_forward_false():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=20)
    fwd = _forward_linear(G)
    res = run_esmda(ens, fwd, d, R, ESMDAConfig(n_iter=3), final_forward=False)
    assert fwd.calls["n"] == 3
    assert len(res.prediction_history) == 3 and res.misfit_history.shape == (3,)


def test_run_esmda_failed_member_is_dropped_for_good():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    res = run_esmda(
        ens, _forward_linear(G, fail={1: [4, 7]}), d, R, ESMDAConfig(n_iter=3, seed=0)
    )
    assert not res.valid[4] and not res.valid[7] and res.valid.sum() == 28
    assert [int(v.sum()) for v in res.valid_history] == [30, 28, 28]
    # dropped members keep the value they had when they failed
    assert np.allclose(res.param[:, 4], res.param_history[1][:, 4])
    assert np.allclose(res.param[:, 4], res.param_history[3][:, 4])
    # their predictions are NaN afterwards, the others finite
    last = res.prediction_history[-1]
    assert np.all(np.isnan(last[:, [4, 7]])) and np.all(np.isfinite(last[:, res.valid]))


def test_run_esmda_too_few_members_raises():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=4)
    with pytest.raises(RuntimeError, match="valid member"):
        run_esmda(ens, _forward_linear(G, fail={0: [0, 1, 2]}), d, R, ESMDAConfig(n_iter=2))


def test_run_esmda_drops_nonfinite_observations():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=200)
    res = run_esmda(
        ens, _forward_linear(G, nan_rows=[1]), d, R, ESMDAConfig(n_iter=3, seed=0)
    )
    assert np.all(np.isfinite(res.param)) and np.all(np.isfinite(res.misfit_history))
    # same result as removing the observation by hand
    keep = [0, 2, 3, 4]
    ref = run_esmda(
        ens, _forward_linear(G[keep]), d[keep], R[np.ix_(keep, keep)], ESMDAConfig(n_iter=3, seed=0)
    )
    assert np.allclose(res.param, ref.param)


def test_run_esmda_all_observations_unusable_raises():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=10)
    with pytest.raises(RuntimeError, match="No usable observation"):
        run_esmda(ens, _forward_linear(G, nan_rows=list(range(d.size))), d, R, ESMDAConfig(n_iter=2))


def test_run_esmda_project_applied_after_each_update():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=40)
    lo, hi = -0.5, 0.5
    res = run_esmda(
        ens, _forward_linear(G), d, R, ESMDAConfig(n_iter=3, seed=0),
        project=lambda p: np.clip(p, lo, hi),
    )
    for p in res.param_history[1:]:
        assert p.min() >= lo and p.max() <= hi


def test_run_esmda_localisation_zero_keeps_prior():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    cfg = ESMDAConfig(n_iter=2, localisation_matrix=np.zeros((ens.shape[0], d.size)))
    res = run_esmda(ens, _forward_linear(G), d, R, cfg)
    assert np.allclose(res.param, ens)


def test_run_esmda_validation():
    G, C0, m0, R, d, ens = _linear_problem(n_ens=10)
    fwd = _forward_linear(G)
    with pytest.raises(ValueError, match="param must be"):
        run_esmda(ens[0], fwd, d, R)
    with pytest.raises(ValueError, match="data_cov"):
        run_esmda(ens, fwd, d, np.ones(2))
    with pytest.raises(ValueError, match="localisation_matrix"):
        run_esmda(ens, fwd, d, R, ESMDAConfig(localisation_matrix=np.ones((1, 1))))
    with pytest.raises(ValueError, match="min_members"):
        run_esmda(ens, fwd, d, R, min_members=1)

    bad_pred = lambda p, v: (np.zeros((3, p.shape[1])), np.ones(p.shape[1], bool))
    with pytest.raises(ValueError, match="predictions of shape"):
        run_esmda(ens, bad_pred, d, R)
    bad_ok = lambda p, v: (G @ p, np.ones(3, bool))
    with pytest.raises(ValueError, match="ok of shape"):
        run_esmda(ens, bad_ok, d, R)
    with pytest.raises(ValueError, match="project returned"):
        run_esmda(ens, fwd, d, R, project=lambda p: p[:1])


# --------------------------------------------------------------------------- progress logging
def test_run_esmda_logs_every_iteration(caplog):
    caplog.set_level(logging.INFO, logger="pyCATHY.DA.esmda")
    G, C0, m0, R, d, ens = _linear_problem(n_ens=30)
    run_esmda(ens, _forward_linear(G, fail={1: [4, 7]}), d, R, ESMDAConfig(n_iter=3, seed=0))
    msgs = [r.getMessage() for r in caplog.records if r.name == "pyCATHY.DA.esmda" and r.levelno == logging.INFO]

    runs = [m for m in msgs if "forward run done" in m]
    assert len(runs) == 3
    assert "iteration 1/3" in runs[0] and "alpha=3" in runs[0] and "30/30 members valid" in runs[0]
    assert "28/30 members valid" in runs[1]                         # two members dropped at the 2nd run
    ups = [m for m in msgs if "updated parameters" in m]
    assert len(ups) == 3 and all("p0 " in u and "+/-" in u for u in ups)
    final = [m for m in msgs if m.startswith("ES-MDA final run")]
    assert len(final) == 1 and "28/30 members valid" in final[0] and "misfit" in final[0]
    # chronological: run 1, update 1, run 2, update 2, run 3, update 3, final
    order = [msgs.index(x) for x in (runs[0], ups[0], runs[1], ups[1], runs[2], ups[2], final[0])]
    assert order == sorted(order)


def test_run_esmda_without_final_forward_logs_no_final_line(caplog):
    caplog.set_level(logging.INFO, logger="pyCATHY.DA.esmda")
    G, C0, m0, R, d, ens = _linear_problem(n_ens=20)
    run_esmda(ens, _forward_linear(G), d, R, ESMDAConfig(n_iter=2), final_forward=False)
    assert not [r for r in caplog.records if r.getMessage().startswith("ES-MDA final run")]
