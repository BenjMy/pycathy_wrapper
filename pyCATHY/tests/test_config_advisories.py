"""Statistical properties of the analysis step, pinned numerically, and the advisories derived from them."""
import numpy as np
import pytest

from pyCATHY.DA import enkf, esmda
from pyCATHY.DA.analysis import (INFLATION_CONVENTION, AnalysisConfig, DAConfig, config_advisories, run_analysis)
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.issues import IssueLog


@pytest.mark.parametrize("m,N,diag", [(5, 20, True), (60, 20, True), (60, 20, False), (300, 16, True)])
def test_ensemble_space_solve_equals_the_direct_solve(m, N, diag):
    rng = np.random.default_rng(0)
    S = rng.normal(size=(m, N)); S -= S.mean(1, keepdims=True)
    A = rng.normal(size=(m, m))
    R = np.diag(rng.uniform(.5, 2, m)) if diag else A @ A.T / m + .5 * np.eye(m)
    rhs = rng.normal(size=(m, 3))
    ref = np.linalg.solve(S @ S.T / (N - 1) + R, rhs)
    np.testing.assert_allclose(enkf._solve_obs_space(S, R, rhs, N), ref, atol=1e-10)


def _linear_gaussian(perturb_obs, N=2000):
    rng = np.random.default_rng(3)
    P = np.diag([1., 2., .5, 1.5]); R = np.eye(4) * .5
    Pa = P - P @ np.linalg.inv(P + R) @ P                       # Kalman posterior covariance
    x = np.linalg.cholesky(P) @ rng.normal(size=(4, N))
    np.random.seed(0)
    out = enkf.enkf_analysis(np.zeros(4), R, np.zeros((0, N)), x, x.copy(), Sakov=True, perturb_obs=perturb_obs)
    return np.var(out[9], axis=1, ddof=1) / np.diag(Pa)


def test_posterior_variance_needs_perturbed_observations():
    """Burgers et al. (1998): without perturbed observations the K R K^T term of the posterior covariance is missing."""
    assert np.all(_linear_gaussian(False) < 0.6)                 # 20-50 % of the Kalman variance
    assert np.all(np.abs(_linear_gaussian(True) - 1) < 0.1)


def test_inflation_conventions_differ_between_filters():
    rng = np.random.default_rng(1)
    N = 400
    x = rng.normal(size=(3, N)); y = rng.normal(size=(2, 3)) @ x
    huge = np.eye(2) * 1e12                                      # no update: the spread change is the inflation alone
    sd = lambda a: a.std(1, ddof=1).mean()
    plain = enkf.enkf_analysis(np.zeros(2), huge, np.zeros((0, N)), x, y, alpha=1.21, Sakov=True)[9]
    assert sd(plain) / sd(x) == pytest.approx(1.1, abs=1e-3)     # variance factor: std x sqrt(1.21)
    smoother = esmda.esmda_update(x, y, np.zeros(2), np.array([1e12, 1e12]), 1.0, rng=rng, inflation=1.21)
    assert sd(smoother) / sd(x) == pytest.approx(1.21, abs=1e-3)  # std factor
    assert "variance" in INFLATION_CONVENTION["enkf_Evensen2009"] and "std factor" in INFLATION_CONVENTION["esmda"]


@pytest.mark.parametrize("da_type", ["enkf_dual_analysis", "enkf_dual_analysis_Sakov"])
def test_dual_filter_ignores_inflation_and_says_so(da_type):
    rng = np.random.default_rng(0)
    state, param, pred = rng.normal(size=(5, 30)), rng.normal(size=(2, 30)), rng.normal(size=(3, 30))
    res = []
    for infl in (1.0, 3.0):
        np.random.seed(0)
        cfg = AnalysisConfig(inflation_states=infl, sakov=True)
        res.append(run_analysis(da_type, np.zeros(3), np.eye(3), param.copy(), ["a", "b"], [state.copy()] * 2, pred.copy(),
                                config=cfg).analysis_param)
    np.testing.assert_allclose(res[0], res[1])                   # the option really has no effect ...
    msgs = " ".join(m for _, m in config_advisories(da_type, AnalysisConfig(inflation_states=3.0)))
    assert "inflation=3 has no effect" in msgs                   # ... and the user is told


def test_advisories_by_filter():
    adv = lambda t, **kw: " | ".join(m for _, m in config_advisories(t, AnalysisConfig(**kw)))
    assert adv("enkf_Evensen2009") == ""                                            # perturb_obs=True is the default (S19)
    assert "without perturbed observations" in adv("enkf_Evensen2009", perturb_obs=False)
    assert adv("esmda") == "" and adv("esmda", inflation_params=1.2) == ""
    assert "never perturbs" in adv("enkf_analysis_inflation")
    assert "no effect" in adv("enkf_analysis_inflation", perturb_obs=True)         # said once, in the filter's own advisory
    assert "no effect with pf" in adv("pf", inflation_states=1.5) and adv("pf") == ""


def test_cycle_config_logs_the_advisories_once_per_distinct_message(tmp_path):
    sim = DA.__new__(DA)
    sim.config, sim.inflation, sim.localisationMatrix, sim.count_DA_cycle = DAConfig(), 1.0, None, 0
    sim._run_da_type = "enkf_Evensen2009"
    sim.config = DAConfig(analysis=AnalysisConfig(perturb_obs=False))        # the advisory only fires on an explicit False
    sim.issue_log = IssueLog(str(tmp_path / "i.log"))
    for c in range(3):
        sim.count_DA_cycle = c
        sim._analysis_config_of_cycle()
    s = sim.issue_log.summary()
    assert s["n_distinct"] == 1 and s["n_events"] == 3 and s["by_kind"] == {"config": 3}     # one message, seen 3 cycles
    assert "perturbed" in open(str(tmp_path / "i.log")).read()


def test_manifest_records_the_inflation_convention(tmp_path):
    sim = DA.__new__(DA)
    sim.config, sim._run_da_type, sim._run_params_updated = DAConfig(), "esmda", ["Ks"]
    assert "compounding" in sim._run_options_manifest()["inflation_convention"]


@pytest.mark.parametrize("n_iter", [1, 2, 4, 8])
def test_esmda_matches_one_kalman_update_in_the_linear_gaussian_case(n_iter):
    """Emerick & Reynolds (2013): sum(1/alpha) = 1 makes the n_iter updates equal one update with the true R."""
    rng = np.random.default_rng(5)
    n, m, N = 3, 4, 20000
    G = rng.normal(size=(m, n)); P = np.diag([1., .5, 2.]); R = np.eye(m) * .3
    m0 = np.array([0., 1., -1.]); d = G @ np.array([.5, .2, 1.])
    K = P @ G.T @ np.linalg.inv(G @ P @ G.T + R)
    mean_k, Pk = m0 + K @ (d - G @ m0), (np.eye(n) - K @ G) @ P
    x = m0[:, None] + np.linalg.cholesky(P) @ rng.normal(size=(n, N))
    for a in esmda.compute_alphas(n_iter):
        x = esmda.esmda_update(x, G @ x, d, np.diag(R), a, rng=rng)
    assert np.abs(x.mean(1) - mean_k).max() < 0.03
    np.testing.assert_allclose(np.var(x, 1, ddof=1) / np.diag(Pk), 1.0, atol=0.05)


def test_particle_filter_survives_a_degenerate_likelihood_and_resamples_exactly():
    from pyCATHY.DA import pf

    rng = np.random.default_rng(0)
    w = rng.dirichlet(np.ones(50)); np.random.seed(1)
    cnt = np.bincount(pf.systematic_resample(w, 50), minlength=50)
    assert np.all((cnt >= np.floor(50 * w)) & (cnt <= np.ceil(50 * w)))      # systematic: floor/ceil of N w
    N, m = 40, 2000                                                          # 2000 data, misfit of ~60 sigma
    out = pf.particle_filter_analysis(np.full(m, 6.0), np.eye(m) * .01, rng.normal(size=(2, N)), rng.normal(size=(30, N)),
                                      rng.normal(size=(m, N)), use_enkf_update=False, resample_threshold=0.5,
                                      jitter_std_param=0.0, jitter_std_state=0.0, temper_neff_frac=0.5, liu_west_a=0.0)
    w = np.asarray(out["weights"])
    assert np.all(np.isfinite(w)) and w.sum() == pytest.approx(1.0)
    assert out["n_eff_raw"] < 1.5 and out["n_eff"] == pytest.approx(0.5 * N)  # exact likelihood degenerate, tempered N_eff ok
    assert np.all(np.isfinite(out["Analysis"])) and np.all(np.isfinite(out["Analysisparam"]))


# ---- realised prior / posterior of the updated parameters
def _entry(prior, finals=(), **kw):
    e = dict(ini_perturbation=np.asarray(prior, float), nominal=1.0, sampling_type="normal", per_type=None,
             clip_min=None, clip_max=None)
    e.update(kw)
    for k, f in enumerate(finals, 1):
        e[f"update_nb{k}"] = np.asarray(f, float)
    return e


def test_parameter_summary_reports_the_realised_prior_not_the_requested_one():
    from pyCATHY.DA.param_tables import parameter_summary

    np.random.seed(0)
    from pyCATHY.DA import perturbate as P
    s = P.sampling_dist("lognormal", 1.0, 1.0, 20000)                     # nominal 1.0 is the MEDIAN
    t = P.sampling_dist_trunc(0.2, 0.8, 20000, loc=0.5, scale=0.3)        # requested sd 0.3
    out = parameter_summary({"Ks0": _entry(s, sampling_type="lognormal"), "n0": _entry(t, clip_min=0.2, clip_max=0.8)})
    assert out["Ks0"]["prior_median"] == pytest.approx(1.0, rel=0.03)
    assert out["Ks0"]["prior_mean_over_nominal"] == pytest.approx(np.exp(0.5), rel=0.05)   # +65 %
    assert out["n0"]["prior_sd"] < 0.18 and out["n0"]["n_updates"] == 0


def test_parameter_flags_catch_collapse_and_bound_pile_up():
    from pyCATHY.DA.param_tables import parameter_flags, parameter_summary

    rng = np.random.default_rng(0)
    prior = rng.normal(1, .3, 40)
    ok = _entry(prior, [rng.normal(1, .2, 40)], clip_min=0.0, clip_max=2.0)
    collapsed = _entry(prior, [rng.normal(1, .2, 40), rng.normal(1.1, .005, 40)])
    stuck = _entry(prior, [np.where(rng.random(40) < .6, 0.0, rng.normal(1, .2, 40))], clip_min=0.0, clip_max=2.0)
    summ = parameter_summary({"ok": ok, "collapsed": collapsed, "stuck": stuck, "scalar": 3.0})
    assert "scalar" not in summ and summ["stuck"]["frac_at_bound_final"] >= 0.5
    flags = dict(parameter_flags(summ))
    assert "ok" not in flags and "collapsed" in flags["collapsed"] and "bound decides" in flags["stuck"]
    assert parameter_flags(summ, collapse=2.0, at_bound=2.0) == []


def test_manifest_has_the_parameter_summary_and_logs_the_flags(tmp_path):
    import json

    sim = DA.__new__(DA)
    sim.project_name, sim.workdir = "prj", str(tmp_path)
    sim.config, sim.count_DA_cycle, sim.NENS = DAConfig(), 2, 40
    sim.dict_parm_pert = {"Ks0": _entry(np.random.default_rng(0).normal(1, .3, 40),
                                        [np.full(40, 1.05) + 1e-4 * np.arange(40)])}
    sim.issue_log = IssueLog(str(tmp_path / "i.log"))
    path = sim.write_run_summary()
    data = json.load(open(path))
    assert data["parameters"]["Ks0"]["contraction"] > 0.95 and data["n_issues"] == 1
    assert "Ks0: ensemble spread contracted" in open(str(tmp_path / "i.log")).read()


# ---- localisation
def test_gaspari_cohn_is_a_valid_compact_correlation_function():
    from pyCATHY.DA.localisation import build_localization_matrix, gaspari_cohn

    r = np.linspace(0, 2.5, 2501)
    w = gaspari_cohn(r, 1.0)
    assert w[0] == 1.0 and gaspari_cohn(np.array([1.0]), 1.0)[0] == pytest.approx(5 / 24)   # published value at r = L
    assert np.all(w[r >= 2.0] == 0) and np.all(np.diff(w) <= 1e-12)                       # support 2L, monotone
    assert abs(gaspari_cohn(np.array([1 - 1e-9]), 1.)[0] - gaspari_cohn(np.array([1 + 1e-9]), 1.)[0]) < 1e-8
    x = np.linspace(0, 10, 200)[:, None]
    assert np.linalg.eigvalsh(build_localization_matrix(x, x, 1.0)).min() > -1e-10         # positive semi-definite
    with pytest.raises(ValueError):
        gaspari_cohn(r, 0.0)


def test_localized_filter_reduces_to_the_plain_filter_and_respects_the_mask():
    rng = np.random.default_rng(0)
    N, ns, m = 30, 12, 5
    state, param = rng.normal(size=(ns, N)), rng.normal(size=(2, N))
    pred, data, R = state[:m] + .1 * rng.normal(size=(m, N)), rng.normal(size=m), np.eye(m) * .2

    def run(t, loc=None):
        np.random.seed(0)
        r = run_analysis(t, data, R, param.copy(), ["a", "b"], [state.copy()] * 2, pred.copy(),
                         config=AnalysisConfig(sakov=True, localisation_matrix=loc))
        return np.asarray(r.analysis), np.asarray(r.analysis_param)

    a0, p0 = run("enkf_Evensen2009_Sakov")
    a1, p1 = run("enkf_analysis_localized_with_inflation", np.ones((ns, m)))
    np.testing.assert_allclose(a1, a0, atol=1e-12); np.testing.assert_allclose(p1, p0, atol=1e-12)   # taper of ones
    az, _ = run("enkf_analysis_localized_with_inflation", np.zeros((ns, m)))
    np.testing.assert_array_equal(az, state)                                                         # taper of zeros
    half = np.ones((ns, m)); half[6:] = 0
    ah, _ = run("enkf_analysis_localized_with_inflation", half)
    assert np.abs(ah[:6] - state[:6]).max() > 1e-6 and np.array_equal(ah[6:], state[6:])           # masked rows untouched


def test_log10_covariance_propagation_matches_monte_carlo():
    from pyCATHY.DA.cathy_DA import _log10_cov

    rng = np.random.default_rng(0)
    d = np.array([20., 100., 400.]); sig = .05 * d
    y = d[:, None] + sig[:, None] * rng.standard_normal((3, 200000))
    np.testing.assert_allclose(np.diag(_log10_cov(d, np.diag(sig ** 2))), np.var(np.log10(y + 1e-3), axis=1), rtol=0.02)


# ---- round S19: perturb_obs=True is the default ------------------------------------------------------------------
def test_perturb_obs_is_on_by_default_and_the_default_run_is_silent(tmp_path):
    assert AnalysisConfig().perturb_obs is True
    sim = DA.__new__(DA)
    sim.config, sim.inflation, sim.localisationMatrix, sim.count_DA_cycle = DAConfig(), 1.0, None, 0
    sim._run_da_type = "enkf_Evensen2009"
    sim.issue_log = IssueLog(str(tmp_path / "i.log"))
    assert sim._analysis_config_of_cycle().perturb_obs is True               # the config built when analysis=None
    assert sim.issue_log.summary()["n_events"] == 0


def test_default_config_gives_the_kalman_posterior_variance():
    """Through the registered adapter with a default AnalysisConfig (not the low-level function)."""
    from pyCATHY.DA.cathy_DA import ANALYSES
    rng = np.random.default_rng(3)
    P = np.diag([1., 2., .5, 1.5]); R = np.eye(4) * .5; N = 2000
    Pa = P - P @ np.linalg.inv(P + R) @ P
    x = np.linalg.cholesky(P) @ rng.normal(size=(4, N))
    np.random.seed(0)
    res = ANALYSES["enkf_Evensen2009"](np.zeros(4), R, np.zeros((0, N)), x, x.copy(), AnalysisConfig())
    ratio = np.var(res.analysis, axis=1, ddof=1) / np.diag(Pa)
    assert np.all(np.abs(ratio - 1) < 0.1)


@pytest.mark.parametrize("da_type", ["enkf_analysis_inflation", "enkf_dual_analysis", "pf"])
def test_default_change_does_not_touch_filters_that_do_not_perturb(da_type):
    """Bit-identical results for perturb_obs True / False: the new default cannot change them."""
    from pyCATHY.DA.cathy_DA import ANALYSES
    rng = np.random.default_rng(0)
    n, m, N = 6, 5, 24
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens, param = rng.normal(size=(n, N)), 0.4 + 0.05 * rng.normal(size=(1, N))
    y = H @ rng.normal(size=n); R = np.eye(m) * .1
    out = []
    for flag in (True, False):
        np.random.seed(5)
        out.append(ANALYSES[da_type](y, R, param.copy(), ens.copy(), H @ ens, AnalysisConfig(perturb_obs=flag)).analysis)
    np.testing.assert_array_equal(out[0], out[1])


def test_manifest_records_the_effective_perturb_obs():
    sim = DA.__new__(DA)
    sim.config, sim._run_params_updated = DAConfig(), []
    sim._run_da_type = "enkf_Evensen2009"
    o = sim._run_options_manifest()
    assert o["perturb_obs"] is True and o["sakov"] is True and o["perturb_obs_applied"] is True
    sim._run_da_type = "enkf_analysis_inflation"
    o = sim._run_options_manifest()
    assert o["perturb_obs"] is True and o["perturb_obs_applied"] is False      # configured, but this filter never perturbs


# ---- round S19c: opt-in centring of the observation noise --------------------------------------------------------
def test_centred_noise_has_zero_member_mean_and_same_random_stream():
    R = np.diag([.5, 1., 2.])
    np.random.seed(0); raw = enkf._sample_obs_noise(R, 16)
    np.random.seed(0); cen = enkf._sample_obs_noise(R, 16, center=True)
    np.testing.assert_allclose(cen.mean(1), 0, atol=1e-14)
    np.testing.assert_allclose(cen, raw - raw.mean(1, keepdims=True))
    np.random.seed(0); again = enkf._sample_obs_noise(R, 16)
    np.testing.assert_array_equal(raw, again)                      # default untouched


def _analysis(center, perturb=True, N=40, seed=0):
    rng = np.random.default_rng(1)
    P = np.diag([1., 2., .5, 1.5]); R = np.eye(4) * .5
    x = np.linalg.cholesky(P) @ rng.normal(size=(4, N))
    y = rng.normal(size=4)
    np.random.seed(seed)
    out = enkf.enkf_analysis(y, R, np.zeros((0, N)), x, x.copy(), Sakov=True,
                             perturb_obs=perturb, perturb_obs_center=center)
    return out[9], x, y, R


def test_centred_noise_gives_the_analysis_mean_of_the_unperturbed_data():
    cen, *_ = _analysis(True)
    plain, *_ = _analysis(False, perturb=False)
    raw, *_ = _analysis(False)
    np.testing.assert_allclose(cen.mean(1), plain.mean(1), atol=1e-10)      # exact, not statistical
    assert np.abs(raw.mean(1) - plain.mean(1)).max() > 1e-3                 # the sampling error that was there


def test_centred_noise_keeps_the_kalman_posterior_variance():
    rng = np.random.default_rng(3)
    P = np.diag([1., 2., .5, 1.5]); R = np.eye(4) * .5; N = 2000
    Pa = P - P @ np.linalg.inv(P + R) @ P
    x = np.linalg.cholesky(P) @ rng.normal(size=(4, N))
    np.random.seed(0)
    out = enkf.enkf_analysis(np.zeros(4), R, np.zeros((0, N)), x, x.copy(), Sakov=True, perturb_obs=True,
                             perturb_obs_center=True)[9]
    assert np.all(np.abs(np.var(out, axis=1, ddof=1) / np.diag(Pa) - 1) < 0.1)


def test_default_is_not_centred_and_reaches_the_adapters():
    from pyCATHY.DA.cathy_DA import ANALYSES
    assert AnalysisConfig().perturb_obs_center is False
    rng = np.random.default_rng(0)
    n, m, N = 6, 5, 24
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens, param = rng.normal(size=(n, N)), 0.4 + 0.05 * rng.normal(size=(1, N))
    y = H @ rng.normal(size=n); R = np.eye(m) * .1
    res = {}
    for key, cfg in (("def", AnalysisConfig()), ("off", AnalysisConfig(perturb_obs_center=False)),
                     ("on", AnalysisConfig(perturb_obs_center=True))):
        np.random.seed(5)
        res[key] = ANALYSES["enkf_Evensen2009"](y, R, param.copy(), ens.copy(), H @ ens, cfg).analysis
    np.testing.assert_array_equal(res["def"], res["off"])
    assert not np.allclose(res["def"], res["on"])
    legacy = AnalysisConfig.from_legacy_kwargs({"perturb_obs_center": True})
    assert legacy.perturb_obs_center is True


def test_localized_filter_forwards_the_centring():
    from pyCATHY.DA.cathy_DA import ANALYSES
    rng = np.random.default_rng(0)
    n, m, N = 6, 5, 24
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens, param = rng.normal(size=(n, N)), 0.4 + 0.05 * rng.normal(size=(1, N))
    y = H @ rng.normal(size=n); R = np.eye(m) * .1
    out = {}
    for c in (False, True):
        np.random.seed(5)
        out[c] = ANALYSES["enkf_analysis_localized_with_inflation"](
            y, R, param.copy(), ens.copy(), H @ ens, AnalysisConfig(perturb_obs_center=c, localisation_matrix=np.ones((n, m)))).analysis
    assert not np.allclose(out[False], out[True])
