"""Round-trip of DA parameter transforms + back-transform of analysed params."""
import ast, logging, numpy as np
src = open("cathy_DA.py").read()
cls = [n for n in ast.parse(src).body if isinstance(n, ast.ClassDef) and n.name == "DA"][0]
meths = [n for n in cls.body if isinstance(n, ast.FunctionDef)
         and n.name in ("transform_parameters", "_back_transform_analysis_param")]
ns = {"np": np, "logger": logging.getLogger("t")}
mod = ast.Module([ast.ClassDef("T", [], [], meths, [])], [])
ast.fix_missing_locations(mod)
exec(compile(mod, "s", "exec"), ns)
T = ns["T"]()
T.dict_parm_pert = {
    "ks": {"transf_type": "log"},
    "poro": {"transf_type": "log-ratio", "transf_bounds": {"A": 0.1, "B": 0.6}},
    "h": {"transf_type": "hyperbolic"},
    "plain": {"transf_type": None},
}
rng = np.random.default_rng(1)
vals = {"ks": rng.uniform(1e-6, 1e-3, 8), "poro": rng.uniform(.15, .55, 8),
        "h": rng.normal(0, 3, 8), "plain": rng.normal(size=8)}
names = ["St. var.", "ks", "poro", "h", "plain"]
fw = np.vstack([T.transform_parameters([n], param=vals[n])[0] for n in names[1:]])  # (n_param, n_mem)
assert np.isfinite(fw).all()
back = T._back_transform_analysis_param(names, fw.T)       # (n_mem, n_param)
for j, n in enumerate(names[1:]):
    assert np.allclose(back[:, j], vals[n], rtol=1e-9), n
# wrong shape -> returned untouched, no crash
assert T._back_transform_analysis_param(names, np.zeros((8, 2))).shape == (8, 2)
assert T._back_transform_analysis_param(names, np.zeros((8, 0))).size == 0
print("transform round-trips OK")
