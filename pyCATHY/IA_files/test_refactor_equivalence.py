"""Original vs refactored Archie_rho_DA (same seed) + _captured_stdout behaviour.
Run: python test_refactor_equivalence.py [path/to/original/petro_Archie.py]"""
import ast, sys, logging, io, contextlib, types
import numpy as np

ORIG = sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/uploads/petro_Archie.py"

def extract(path, name, ns):
    src = open(path).read()
    node = [n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == name][0]
    exec(compile(ast.Module([node], []), path, "exec"), ns)
    return ns[name]

class _Quiet:                       # stands in for rich.console.Console
    def __init__(self, *a, **k): pass
    def print(self, *a, **k): pass
    rule = print
rich = types.SimpleNamespace(console=types.SimpleNamespace(Console=_Quiet))
ns_old = {"np": np, "rich": rich, "print": lambda *a, **k: None}
ns_new = {"np": np, "logger": logging.getLogger("t")}
old = extract(ORIG, "Archie_rho_DA", ns_old)
new = extract("petro_Archie.py", "Archie_rho_DA", ns_new)

rng = np.random.default_rng(0)
sat = [rng.uniform(.2, .9, 50), rng.uniform(.2, .9, 50)]
por = [np.full(50, .4), np.full(50, .35)]
kw = dict(rFluid_Archie=[17., 20.], a_Archie=[1., 1.], m_Archie=[1.7, 2.], n_Archie=[2., 2.])
sat2 = [rng.uniform(.2, .9, (50, 4)), rng.uniform(.2, .9, (50, 4))]
por2 = [np.full((50, 4), .4), np.full((50, 4), .35)]
for pert in ([0, 0],):   # noise-free: identical to the original (the noise formula was changed on purpose)
    for sat_, por_ in ((sat, por), (sat2, por2)):   # 1-D and 2-D sigma
        p = pert
        np.random.seed(3); r_old = old(sat=sat_, porosity=por_, pert_sigma_Archie=p, **kw)
        np.random.seed(3); r_new = new(sat=sat_, porosity=por_, pert_sigma_Archie=p, **kw)
        for a, b in zip(np.atleast_1d(r_old) if not isinstance(r_old, tuple) else r_old,
                        np.atleast_1d(r_new) if not isinstance(r_new, tuple) else r_new):
            assert np.allclose(np.asarray(a), np.asarray(b), equal_nan=True), p
print("Archie_rho_DA: identical outputs (no noise)")

# noise: relative std = pert, positive conductivity, independent per node and per time row
for sat_, por_ in ((sat, por), (sat2, por2)):
    base = new(sat=sat_, porosity=por_, pert_sigma_Archie=[0, 0], **kw)
    noisy = np.concatenate([np.ravel(x) for x in [new(sat=sat_, porosity=por_, pert_sigma_Archie=[0.05, 0.05], **kw)]])
    b0 = np.ravel(base); rel = (1 / noisy) / (1 / b0) - 1                  # relative change of sigma
    assert 0.03 < np.std(rel) < 0.08 and np.all(noisy > 0), np.std(rel)
    assert not np.allclose(rel, rel[0])
print("Archie_rho_DA: noise std ~ pert * sigma")

# _captured_stdout
src = open("simulate_ERT.py").read()
tree = ast.parse(src)
fn = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_captured_stdout"][0]
ns = {"contextlib": contextlib, "io": io, "logging": logging, "sys": sys,
      "os": __import__("os"), "logger": logging.getLogger("t2")}
exec(compile(ast.Module([fn], []), "s", "exec"), ns)
cs = ns["_captured_stdout"]
orig = sys.stdout
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    with cs(True): print("inside-verbose")
    with cs(False): print("inside-muted")
    try:
        with cs(False): print("before-exc"); raise ValueError("x")
    except ValueError: pass
    print("after")
out = buf.getvalue()
assert "inside-muted" not in out and "before-exc" not in out and "after" in out, out
assert sys.stdout is orig
print("_captured_stdout: captures, restores after exception")
