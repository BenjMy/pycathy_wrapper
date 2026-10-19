import numpy as np
import pandas as pd

from pyCATHY.ERT.ERT_checks import quadrupole_quality_mask
from pyCATHY.DA.cathy_DA import DA


def _df():
    rng = np.random.default_rng(0)
    n = 100
    k = 10 ** rng.normal(1.5, 0.1, n)
    rhoa = 15 + rng.normal(0, 0.5, n)
    k[-5:] = 3e5          # k_too_large
    rhoa[0] = -1.0        # non positive
    rhoa[1] = np.nan
    valid = np.ones(n)
    valid[2] = 0
    return pd.DataFrame({"rhoa": rhoa, "k": k, "valid": valid})


def test_mask_flags_each_problem():
    keep, reasons = quadrupole_quality_mask(_df())
    assert keep.size == 100
    assert not keep[[0, 1, 2]].any() and not keep[-5:].any()
    assert keep[3:-5].all()
    assert {"rhoa_invalid", "flagged_invalid", "k_too_large"} <= set(reasons)


def test_mask_without_k_still_filters_rhoa():
    keep, _ = quadrupole_quality_mask(pd.DataFrame({"rhoa": [1.0, -2.0, np.inf]}))
    assert keep.tolist() == [True, False, False]


def test_signal_below_noise_optional():
    df = pd.DataFrame({"rhoa": [15.0, 15.0], "k": [10.0, 1e4]})
    keep, _ = quadrupole_quality_mask(df)
    assert keep.all()
    keep, r = quadrupole_quality_mask(df, noise_abs=1e-2)  # R=1.5e-3 V < 10*1e-2
    assert not keep[1] and "signal_below_noise" in r


def test_da_mask_layout_and_switch():
    df = _df()
    da = type("X", (), {})()
    da.count_DA_cycle = 0
    da.dict_obs = {0: {"ERT0": {"data": df}}}
    da.mask_bad_ert = True
    da.ert_check_reports = []
    keep = DA._ert_quality_mask(da, ["ERT0"], 100, np.ones(100, dtype=bool))
    assert keep.sum() == 92
    assert DA._ert_quality_mask(da, ["ERT0"], 50) is None  # size mismatch: nothing dropped
    da.mask_bad_ert = False
    assert DA._ert_quality_mask(da, ["ERT0"], 100) is None


def test_exclusions_are_flagged_not_deleted():
    """The observed data stay complete; the excluded indices are recorded."""
    df = _df()
    before = df.copy()
    da = type("X", (), {})()
    da.count_DA_cycle = 3
    da.dict_obs = {t: {"ERT0": {"data": df}} for t in (0, 100, 200, 300)}
    da.mask_bad_ert = True
    da.ert_check_reports = []
    keep = DA._ert_quality_mask(da, ["ERT0"], 100, None)
    pd.testing.assert_frame_equal(df, before)  # nothing removed or rewritten
    rec = da.ert_check_reports[0]
    assert rec["stage"] == "mask" and rec["cycle"] == 3
    assert rec["n_total"] == 100 and rec["n_excluded"] == 8
    assert rec["excluded_idx"] == np.flatnonzero(~keep).tolist()
    assert "ok" not in rec  # not counted as a failed check


# ---- the covariance-localisation matrix has one column per observation -------------------
def test_localisation_matrix_follows_the_quadrupole_mask(caplog):
    """With ERT, a matrix built for all the data failed in the filter once one quadrupole was
    excluded (ValueError: Localization matrix shape (10, 20) does not match (10, 18))."""
    import logging
    from types import SimpleNamespace

    from pyCATHY.DA.cathy_DA import ANALYSES, AnalysisConfig

    rng = np.random.default_rng(0)
    n, m, nens = 10, 20, 15
    H = rng.normal(size=(m, n)) / np.sqrt(n)
    ens = rng.normal(size=(n, nens))
    param = 0.4 + 0.05 * rng.normal(size=(1, nens))
    y, R, pred = H @ rng.normal(size=n), np.eye(m) * 0.1, H @ ens
    keep = np.ones(m, bool)
    keep[[3, 7]] = False
    full = np.arange(n * m, dtype=float).reshape(n, m) / (n * m)

    cfg = AnalysisConfig(localisation_matrix=full)
    cfg2 = DA._mask_localisation_matrix(SimpleNamespace(), cfg, keep)
    assert cfg2.localisation_matrix.shape == (n, 18)
    np.testing.assert_array_equal(cfg2.localisation_matrix, full[:, keep])
    ANALYSES["enkf_analysis_localized_with_inflation"](
        y[keep], R[np.ix_(keep, keep)], param, ens, pred[keep], cfg2)  # runs

    # already built for the kept observations: unchanged
    kept_only = full[:, keep]
    same = DA._mask_localisation_matrix(SimpleNamespace(), AnalysisConfig(localisation_matrix=kept_only), keep)
    assert same.localisation_matrix is kept_only or np.array_equal(same.localisation_matrix, kept_only)

    # neither size: a warning, matrix untouched
    with caplog.at_level(logging.WARNING):
        odd = DA._mask_localisation_matrix(SimpleNamespace(), AnalysisConfig(localisation_matrix=np.ones((n, 5))), keep)
    assert odd.localisation_matrix.shape == (n, 5) and "left unchanged" in caplog.text

    # no matrix: nothing to do
    assert DA._mask_localisation_matrix(SimpleNamespace(), AnalysisConfig(), keep).localisation_matrix is None
