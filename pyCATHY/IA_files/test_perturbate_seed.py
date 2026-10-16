"""sampling_dist: different parameters must not share the same random draws."""
import ast, numpy as np
src = open("perturbate.py").read()
fn = [n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "sampling_dist"][0]
ns = {"np": np}; exec(compile(ast.Module([fn], []), "s", "exec"), ns); f = ns["sampling_dist"]
a = f("normal", 0.0, 1.0, 64); b = f("normal", 0.0, 1.0, 64)
assert not np.allclose(a, b) and abs(np.corrcoef(a, b)[0, 1]) < 0.5
ens = np.vstack([np.log(f("lognormal", m, 0.3, 64)) for m in (1e-5, 0.4, 17.0, 2.0)])
assert np.linalg.matrix_rank(ens - ens.mean(1, keepdims=True)) == 4        # not rank-1
np.random.seed(5); x = f("normal", 0, 1, 8); np.random.seed(5); y = f("normal", 0, 1, 8)
assert np.allclose(x, y)                                                     # seed set once by the user still reproduces
print("sampling_dist: independent draws per parameter, reproducible via a user seed")
