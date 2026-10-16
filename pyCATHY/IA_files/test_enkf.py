"""enkf.py: refactored (logger) vs original outputs, helper vs dense, collapse warning."""
import sys, types, logging, importlib.util, io, contextlib
import numpy as np
# stub rich (not installed everywhere)
rich = types.ModuleType("rich"); rich.print = lambda *a, **k: None
rc = types.ModuleType("rich.console"); rp = types.ModuleType("rich.progress"); rp.track = lambda x, *a, **k: x
rich.console = rc; sys.modules.update({"rich": rich, "rich.console": rc, "rich.progress": rp})

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m
new = load("enkf.py", "enkf_new"); old = load("orig_repo/enkf.py", "enkf_old")

rng = np.random.default_rng(0)
n_state, n_obs, n_ens, n_par = 60, 7, 12, 2
ens = rng.normal(size=(n_state, n_ens)); par = rng.normal(size=(n_par, n_ens))
pred = rng.normal(size=(n_obs, n_ens)); data = rng.normal(size=n_obs)
R = np.diag(rng.uniform(.5, 1, n_obs))
for sakov in (False, True):
    with contextlib.redirect_stdout(io.StringIO()):
        r_old = old.enkf_analysis(data, R, par, ens, pred, Sakov=sakov)
    r_new = new.enkf_analysis(data, R, par, ens, pred, Sakov=sakov)
    for a, b in zip(r_old, r_new):
        assert np.allclose(a, b)
print("enkf_analysis: identical outputs")

# helper == dense
for P in (rng.normal(size=(30, 9)), rng.normal(size=(5, 4))):
    C = P @ P.T / (P.shape[1] - 1)
    d, o = new._cov_diag_offdiag_means(P)
    assert np.isclose(d, np.mean(np.diag(C))) and np.isclose(o, np.mean(C - np.diag(np.diag(C))))
with contextlib.redirect_stdout(io.StringIO()):
    s_old = old.analyze_covariances(R, pred - pred.mean(1, keepdims=True), ens - ens.mean(1, keepdims=True))
s_new = new.analyze_covariances(R, pred - pred.mean(1, keepdims=True), ens - ens.mean(1, keepdims=True))
for k in s_old: assert np.allclose(s_old[k], s_new[k]), k
print("analyze_covariances: identical summary, no dense n_state^2 matrix")

# warnings go to the logger: no collapse warning for a spread ensemble, one for identical members
class H(logging.Handler):
    def __init__(s): super().__init__(); s.m = []
    def emit(s, r): s.m.append((r.levelno, r.getMessage()))
h = H(); new.logger.addHandler(h); new.logger.setLevel(logging.INFO)
new.enkf_analysis(data, R, par, ens, pred)
assert not any("collapse" in m for _, m in h.m), h.m
pred_c = np.tile(pred[:, :1], (1, n_ens)) + 0.0
try: new.enkf_analysis(data, R, par, ens, pred_c)
except np.linalg.LinAlgError: pass
assert any("collapse" in m for _, m in h.m), h.m
assert any("N_ens=12" in m for _, m in h.m)
print("enkf logging OK")
