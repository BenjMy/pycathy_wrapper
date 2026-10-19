"""Unit tests of the EnKF analysis step (no CATHY run needed).

Synthetic linear problem with many strongly correlated observations (like a
whole ERT survey of a quasi-homogeneous medium). The reference is the exact
ensemble Kalman gain K = P_xo (S S^T/(N-1) + R)^-1.
"""
import importlib.util
import pathlib

import numpy as np
import pytest

try:
    from pyCATHY.DA import enkf
except ImportError:  # heavy geo dependencies missing: load the module alone
    _path = pathlib.Path(__file__).resolve().parents[1] / "DA" / "enkf.py"
    _spec = importlib.util.spec_from_file_location("enkf_standalone", _path)
    enkf = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(enkf)


def _problem(n_state=50, n_obs=2000, n_ens=30, sigma=1.0, seed=0, n_param=1):
    rng = np.random.default_rng(seed)
    # correlated rows: every observation is a mix of a few state modes
    H = rng.normal(size=(n_obs, 3)) @ rng.normal(size=(3, n_state))
    H /= np.sqrt(n_state)
    x_true = rng.normal(size=(n_state, 1))
    ens = rng.normal(size=(n_state, n_ens))
    # parameter that the data cannot see (like ZROOT for an ERT survey)
    param = 0.4 + 0.005 * rng.normal(size=(n_param, n_ens))
    predict_obs = H @ ens
    y = H @ x_true
    data = y + sigma * rng.normal(size=(n_obs, n_ens))  # perturbed observations
    R = sigma**2 * np.eye(n_obs)
    return data, R, param, ens, predict_obs


def _run(sakov, **kw):
    data, R, param, ens, pred = _problem(**kw)
    out = enkf.enkf_analysis(data, R, param, ens, pred, Sakov=sakov)
    return out[-2], out[-1], param, ens  # analysis, analysis_param (N x n_param)


def test_sakov_equals_full_gain():
    """The ensemble-space form must give the same analysis as the full form."""
    ana_s, par_s, _, _ = _run(True)
    ana_f, par_f, _, _ = _run(False)
    np.testing.assert_allclose(ana_s, ana_f, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(par_s, par_f, rtol=1e-6, atol=1e-8)


@pytest.mark.parametrize("n_obs", [10, 200, 2000])
def test_sakov_gain_does_not_scale_with_n_obs(n_obs):
    """The update of an unobservable parameter stays within a few prior sd."""
    _, par_s, param, _ = _run(True, n_obs=n_obs)
    prior_sd = param.std()
    shift = np.abs(par_s.T - param).max()
    assert shift < 5 * prior_sd, f"param shift {shift:.3g} >> prior sd {prior_sd:.3g}"


def test_sakov_state_update_is_bounded():
    """Analysis states stay of the order of the prior (no blow-up)."""
    ana_s, _, _, ens = _run(True)
    assert np.abs(ana_s).max() < 10 * np.abs(ens).max()


def test_sakov_single_observation_matches_full():
    """With one observation the two forms were always equivalent up to H P H^T."""
    ana_s, _, _, _ = _run(True, n_obs=1, sigma=10.0)
    ana_f, _, _, _ = _run(False, n_obs=1, sigma=10.0)
    np.testing.assert_allclose(ana_s, ana_f, rtol=1e-6, atol=1e-8)


def test_solve_obs_space_both_branches_agree():
    """m < N (direct) and m > N (Woodbury) give the same solution as the dense solve."""
    rng = np.random.default_rng(3)
    for m, n_ens in [(8, 30), (200, 15)]:
        S = rng.normal(size=(m, n_ens))
        R = np.diag(rng.uniform(0.5, 2.0, m))
        rhs = rng.normal(size=(m, n_ens))
        ref = np.linalg.solve(S @ S.T / (n_ens - 1) + R, rhs)
        np.testing.assert_allclose(enkf._solve_obs_space(S, R, rhs, n_ens), ref, rtol=1e-8, atol=1e-10)
        Rfull = R + 0.1 * np.ones((m, m))  # non-diagonal R
        ref = np.linalg.solve(S @ S.T / (n_ens - 1) + Rfull, rhs)
        np.testing.assert_allclose(enkf._solve_obs_space(S, Rfull, rhs, n_ens), ref, rtol=1e-8, atol=1e-10)


def test_perturbed_obs_restores_posterior_spread():
    """Without perturbed observations the analysis spread is short by K R K^T;
    with perturb_obs=True it matches the exact (Kalman) posterior."""
    rng = np.random.default_rng(0)
    n, m, N = 5, 40, 4000
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    R = np.eye(m)
    ens = rng.normal(size=(n, N))
    y = H @ rng.normal(size=n) + rng.normal(size=m)
    P = np.cov(ens)
    K = P @ H.T @ np.linalg.inv(H @ P @ H.T + R)
    exact = np.trace((np.eye(n) - K @ H) @ P)
    no_param = np.zeros((0, N))

    np.random.seed(1)
    plain = enkf.enkf_analysis(y, R, no_param, ens, H @ ens, Sakov=True)[-2]
    np.random.seed(1)
    pert = enkf.enkf_analysis(y, R, no_param, ens, H @ ens, Sakov=True, perturb_obs=True)[-2]

    assert np.trace(np.cov(plain)) < 0.3 * exact          # collapsed spread (the old behaviour)
    assert abs(np.trace(np.cov(pert)) / exact - 1) < 0.1  # statistically consistent
    # the mean update is the same in both forms
    np.testing.assert_allclose(pert.mean(1), plain.mean(1), atol=0.05)


def test_perturb_obs_is_off_by_default_and_reproducible():
    a1, _, _, _ = _run(True)
    a2, _, _, _ = _run(True)
    np.testing.assert_allclose(a1, a2)  # default: deterministic


@pytest.mark.parametrize("name", ["enkf_analysis_inflation", "enkf_analysis_inflation_multiparm"])
def test_inflation_variants_reduce_to_plain_enkf_when_alpha_is_one(name):
    """alpha = 1 must give the plain analysis (the ensemble mean used to be scrambled)."""
    rng = np.random.default_rng(1)
    n, m, N = 20, 60, 12
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens = rng.normal(size=(n, N))
    param = 0.4 + 0.01 * rng.normal(size=(1, N))
    y = H @ rng.normal(size=n) + rng.normal(size=m)
    R = np.eye(m)
    ref = enkf.enkf_analysis(y, R, param, ens, H @ ens, Sakov=False)
    out = getattr(enkf, name)(y, R, param, ens, H @ ens, Sakov=False, alpha=1.0)
    np.testing.assert_allclose(out[-2], ref[-2], atol=1e-10)


@pytest.mark.parametrize(
    "da_type",
    ["enkf_Evensen2009", "enkf_Evensen2009_Sakov", "enkf_analysis_localized_with_inflation"],
)
def test_perturb_obs_reaches_every_filter_that_supports_it(da_type):
    """Through the registered adapters (the DA_type entry points): the option is not ignored."""
    from pyCATHY.DA.cathy_DA import ANALYSES, AnalysisConfig

    rng = np.random.default_rng(0)
    n, m, nens = 10, 15, 20
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens = rng.normal(size=(n, nens))
    param = 0.4 + 0.05 * rng.normal(size=(1, nens))
    y = H @ rng.normal(size=n)
    R = np.eye(m) * 0.1

    def run(perturb_obs):
        np.random.seed(7)
        cfg = AnalysisConfig(perturb_obs=perturb_obs)
        if da_type == "enkf_analysis_localized_with_inflation":
            cfg.localisation_matrix = np.ones((n, m))
        return ANALYSES[da_type](y, R, param, ens, H @ ens, cfg).analysis

    off1, off2, on = run(False), run(False), run(True)
    np.testing.assert_allclose(off1, off2)       # default: deterministic
    assert not np.allclose(off1, on)             # the option changes the analysis
