"""_log10_cov: first-order propagation, finite for non-positive observations."""
import ast, numpy as np
src = open("cathy_DA.py").read()
ns = {"np": np}
for n in ast.parse(src).body:
    if (isinstance(n, ast.FunctionDef) and n.name == "_log10_cov") or \
       (isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "_LOG_EPS"):
        exec(compile(ast.Module([n], []), "s", "exec"), ns)
f = ns["_log10_cov"]
d = np.array([50., 100., 200.]); R = np.diag((0.1 * d) ** 2)                 # 10 % error
assert np.allclose(np.sqrt(np.diag(f(d, R))), 0.1 / np.log(10), rtol=1e-3)
Rn = f(np.array([-5., 0., 10.]), np.eye(3))
assert np.isfinite(Rn).all() and np.all(np.diag(Rn) > 0)
print("_log10_cov OK")
