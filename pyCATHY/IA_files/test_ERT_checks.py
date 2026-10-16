"""Assertion-based smoke tests of ERT_checks (no pygimli / pyCATHY needed)."""
import json, pickle, warnings
import numpy as np
import pandas as pd
import ERT_checks as C

rng = np.random.default_rng(0)
N = 200
a = np.arange(N); b = a + 1; m = a + 2; n = a + 3
k = 10 ** rng.normal(1.5, 0.2, N)
rho = 50.0
rhoa = rho * (1 + 0.02 * rng.standard_normal(N))
good = pd.DataFrame(dict(a=a, b=b, m=m, n=n, k=k, rhoa=rhoa, valid=np.ones(N)))

def codes(rep): return {i["code"] for i in rep.issues}

# clean data -> no issue
assert C.check_geometric_factor(k).ok
assert C.check_rhoa(rhoa, np.full(10, rho)).ok
assert C.check_forward_data(good, res0=np.full(10, rho), noise_abs=1e-4).ok
assert C.check_scheme(a, b, m, n, n_elec=N + 3).ok
assert C.check_observations(good.assign(r=good.rhoa / good.k, rec_err=1.0)).ok

# singular quadrupoles + negative / exploding rhoa (the symptom of the log)
k_bad = k.copy(); k_bad[[3, 50]] = [1e-9, 1e9]
rep = C.check_geometric_factor(k_bad)
assert {"k_outlier", "k_too_large"} <= codes(rep) and {3, 50} <= set(rep.issues[0]["idx"])
rhoa_bad = rhoa.copy(); rhoa_bad[[7, 8]] = [-2.3e6, 2.4e7]
rep = C.check_rhoa(rhoa_bad, np.full(10, rho))
assert {"rhoa_nonpositive", "rhoa_implausible"} <= codes(rep)
bad = good.assign(k=k_bad, rhoa=rhoa_bad, valid=np.r_[0, np.ones(N - 1)])
assert {"k_outlier", "rhoa_nonpositive", "valid_flag"} <= codes(C.check_forward_data(bad, res0=[rho]))

# noise-dominated voltage
tiny = good.assign(rhoa=good.k * 1e-6)
assert "signal_below_noise" in codes(C.check_forward_data(tiny, noise_abs=1e-4))

# scheme: -1 shift, shared electrode, out of range, duplicate
rep = C.check_scheme([0, 1, 5, 5], [-1, 2, 6, 6], [2, 1, 7, 7], [3, 4, 99, 99], n_elec=50)
assert {"scheme_negative_index", "scheme_shared_electrode",
        "scheme_index_out_of_range", "scheme_duplicates"} <= codes(rep)

# forward vs observations: off-by-one and k mismatch
shifted = good.assign(a=good.a - 1, k=good.k * 1.5)
assert {"fwd_obs_quadrupoles", "fwd_obs_k"} <= codes(C.check_forward_vs_obs(shifted, good))
assert C.check_forward_vs_obs(good, good).ok
assert "fwd_obs_size" in codes(C.check_forward_vs_obs(good.iloc[:10], good))

# observations
o = good.assign(r=good.rhoa / good.k * 2, rec_err=np.r_[0.0, np.ones(N - 1)])
assert {"obs_rhoa_k_r", "obs_error_invalid"} <= codes(C.check_observations(o))
assert "obs_no_rhoa" in codes(C.check_observations(pd.DataFrame(dict(a=[1]))))

# data covariance units
assert "cov_invalid" in codes(C.check_data_cov([0.0, 1.0], [1.0, 1.0]))
assert "cov_implausible" in codes(C.check_data_cov(np.full(5, 1e6), np.full(5, 50.0)))

# ensemble: duplicate member, regime change, nan, nonpositive, collapse
ens = np.column_stack([rhoa * (1 + 0.1 * rng.standard_normal()) * (1 + 0.05 * rng.standard_normal(N))
                       for _ in range(8)])
assert C.check_ensemble_predictions(ens, rhoa).ok
ens2 = ens.copy(); ens2[:, 5] = ens2[:, 2]; ens2[:, 7] *= 50; ens2[3, 0] = -1.0; ens2[4, 1] = np.nan
rep = C.check_ensemble_predictions(ens2, rhoa, ens_ids=range(1, 9))
assert {"member_duplicated", "member_regime_change", "pred_nonpositive", "pred_nan"} <= codes(rep)
assert "ensemble_collapse" in codes(C.check_ensemble_predictions(np.tile(rhoa[:, None], 4)))
assert "systematic_misfit" in codes(C.check_ensemble_predictions(ens * 1e4, rhoa))
rank1 = np.outer(rhoa, 1 + 0.1 * rng.standard_normal(6))
assert "rank1_ensemble" in codes(C.check_ensemble_predictions(rank1, rhoa))

# ensemble of forward DataFrames: geometry / k must be identical
members = [good.assign(rhoa=r) for r in ens.T]
assert C.check_ensemble_forward(members, good, ens_ids=range(1, 9)).ok
members[3] = members[3].assign(k=members[3].k * 2)
members[4] = members[4].assign(a=members[4].a + 1)
assert {"ensemble_k_differs", "ensemble_quadrupoles_differ"} <= codes(
    C.check_ensemble_forward(members, good))

# innovation and the log / covariance units mismatch
R = np.full(N, 25.0)
assert C.check_innovation(rhoa, ens, R).ok
assert "innovation_large" in codes(C.check_innovation(rhoa, ens * 1e4, R))
assert "R_not_log_transformed" in codes(C.check_innovation(rhoa, ens, R, use_log=True))
full = C.check_assimilation_inputs(np.r_[-1.0, rhoa[1:]], ens, np.diag(R), use_log=True)
assert {"obs_invalid", "R_not_log_transformed"} <= codes(full)

# Archie inputs
assert C.check_archie_inputs([0.3, 0.5], [0.4, 0.4], {"rFluid_Archie": [17.0]}).ok
assert {"sat_out_of_range", "porosity_out_of_range", "archie_param_invalid"} <= codes(
    C.check_archie_inputs([0.0, 1.5], [0.0, 0.4], {"m_Archie": [-1.7]}))

# warn(): logger.warning, one record per issue; Python warnings only if enabled
import logging
class _H(logging.Handler):
    def __init__(self): super().__init__(); self.recs=[]
    def emit(self, r): self.recs.append(r)
h=_H(); C.logger.addHandler(h); C.logger.setLevel(logging.DEBUG)
C.check_geometric_factor(k_bad, name="member 3").warn()
assert len(h.recs)==2 and all(r.levelno==logging.WARNING for r in h.recs)
assert "member 3" in h.recs[0].getMessage()
h.recs.clear(); C.check_geometric_factor(k).warn(); assert not h.recs
C.EMIT_PYTHON_WARNINGS=True
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always"); C.check_geometric_factor(k_bad).warn()
assert len(w)==2 and all(issubclass(x.category, C.ERTCheckWarning) for x in w)
C.EMIT_PYTHON_WARNINGS=False; C.logger.removeHandler(h)
assert C.issues_summary([]) == "OK"

# bookkeeping: attrs survive pickling (multiprocessing), JSON export, failure count
df = good.copy(); df.attrs["ert_check"] = C.check_geometric_factor(k_bad).to_dict()
assert pickle.loads(pickle.dumps(df)).attrs["ert_check"]["ok"] is False
reports = [{"stage": "forward", "cycle": 0, "members": {1: df.attrs["ert_check"], 2: C.check_geometric_factor(k).to_dict()}}]
assert C.count_failed(reports) == 1
C.save_reports(reports, "/tmp/_ert_checks_test.json")
assert json.load(open("/tmp/_ert_checks_test.json"))[0]["members"]["1"]["ok"] is False

# degenerate input must not crash
for bad_in in (None, [], [np.nan, np.nan]):
    C.check_geometric_factor(bad_in); C.check_rhoa(bad_in)
C.check_ensemble_predictions(np.ones((5, 1)))
print("all ERT_checks smoke tests passed")
