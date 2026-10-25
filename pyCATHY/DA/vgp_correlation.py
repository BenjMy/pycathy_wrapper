"""Correlated van Genuchten parameters (Botto, Belluco & Camporese 2018, HESS 22, 4251-4266).

``perturbate`` draws every van Genuchten parameter (VGP) independently. Real VGP are correlated, and the ensemble
filters only know the correlation that the prior ensemble contains. Following Botto et al. (2018):

1. each VGP is mapped to a (nearly) Gaussian variable ``y = T(x)`` by a Johnson transformation
   (``log`` -> LN, ``log-ratio`` with the bounds A, B -> SB, ``hyperbolic`` -> SU, ``None`` -> identity), the same
   transformations as ``transf_type`` of the scenario;
2. the transformed variables are correlated through a correlation matrix ``R`` (3 x 3 for SATCELL_VG, n_VG, thetar_VG):
   ``Sigma = D R D`` with ``D = diag(s)``, ``s_i`` the standard deviation of ``y_i``;
3. the members are ``y = u + sigma_VG * L z``, ``u = T(nominal)``, ``L`` the Cholesky factor of ``Sigma`` and
   ``z ~ N(0, I)`` (Botto: ``y = u + T' z`` with T the upper Cholesky factor), then back-transformed, ``x = T^-1(y)``,
   and clipped to the bounds of the scenario.

Two entry points:

* ``apply_botto_vgp``: Botto et al. (2018) exactly: the upper-triangular Cholesky factors ``T`` of the covariance of the
  transformed (theta_r, alpha, n) of Table 2 (``BOTTO_2018_T2``: sandy loam, silt loam; derived from Carsel & Parrish 1988)
  and ``y = u + sigma_VG T' z`` with sigma_VG = 0.1 (VGP updated) or 0.25 (not updated). The spread comes from ``T``, not
  from the scenario; ``u`` is the transformed nominal of the scenario; the scenario bounds are the Johnson A, B bounds and
  the clip.  CATHY's ``SATCELL_VG`` is 1/alpha (m), so alpha = 1/SATCELL_VG is what is transformed.
* ``apply_vgp_correlation``: generic, a correlation matrix of your own (``R``) with the scenario standard deviations.

Usage (after the independent draw, before the cache / the run)::

    var_per, report = apply_vgp_correlation(var_per, R)             # one matrix for every layer/material group
    var_per, report = apply_vgp_correlation(var_per, {0: R_A, 1: R_B, 2: R_C})   # one per group

The marginal sd of every parameter is the one of the scenario (``sd`` / ``per_sigma``, mapped to the transformed space
by a finite difference); ``sigma_vg`` rescales the whole spread (Botto: 0.1 when the VGP are updated, 0.25 when not).
"""
import logging
import re
import warnings

import numpy as np

logger = logging.getLogger(__name__)

VGP_NAMES = ("SATCELL_VG", "n_VG", "thetar_VG")

# Botto, Belluco & Camporese (2018), HESS 22, 4251, Table 2: factored covariance matrices (upper-triangular T, Sigma = T'T)
# of the transformed (theta_r, alpha, n), rows/columns in this order.
BOTTO_ORDER = ("thetar_VG", "alpha", "n_VG")
BOTTO_2018_T2 = {
    "sandy_loam": np.array([[0.538, 0.017, -0.194], [0.0, 0.014, 0.019], [0.0, 0.0, 0.108]]),
    "silt_loam": np.array([[0.522, 0.030, -0.170], [0.0, 0.082, 0.234], [0.0, 0.0, 0.158]]),
}
# Johnson transformation of each variable (Botto: SB (log-ratio) for all, except n of the sandy loam and alpha of the silt
# loam, which are LN): scenario ``transf_type`` names.
BOTTO_2018_TRANSFORMS = {
    "sandy_loam": ("log-ratio", "log-ratio", "log"),
    "silt_loam": ("log-ratio", "log", "log-ratio"),
}
SIGMA_VG_UPDATED, SIGMA_VG_NOT_UPDATED = 0.1, 0.25      # Botto: VGP updated / not updated by the filter
_KIND = {None: "None", "None": "None", "log": "LN", "log-ratio": "SB", "hyperbolic": "SU"}


def forward(x, kind, A=None, B=None):
    """Johnson transform of the scenario ``transf_type`` ('log', 'log-ratio', 'hyperbolic', None)."""
    k = _KIND[kind]
    x = np.asarray(x, dtype=float)
    if k == "None":
        return x
    if k == "LN":
        return np.log(x)
    if k == "SB":
        return np.log((x - A) / (B - x))
    return np.arcsinh(x)


def inverse(y, kind, A=None, B=None):
    k = _KIND[kind]
    y = np.asarray(y, dtype=float)
    if k == "None":
        return y
    if k == "LN":
        return np.exp(y)
    if k == "SB":
        u = np.exp(y)
        return (A + B * u) / (1.0 + u)
    return np.sinh(y)


def check_correlation(R, name="R"):
    """Validated copy of a correlation matrix: square, symmetric, unit diagonal, |r| <= 1, positive definite.
    A matrix that is only positive semi-definite / slightly indefinite is repaired (eigenvalues floored) with a warning."""
    R = np.array(R, dtype=float)
    if R.ndim != 2 or R.shape[0] != R.shape[1]:
        raise ValueError(f"{name}: a square matrix is required, got shape {R.shape}")
    if not np.allclose(R, R.T, atol=1e-8):
        raise ValueError(f"{name} is not symmetric")
    if not np.allclose(np.diag(R), 1.0, atol=1e-8):
        raise ValueError(f"{name}: the diagonal must be 1 (a CORRELATION matrix)")
    if np.any(np.abs(R) > 1.0 + 1e-12):
        raise ValueError(f"{name}: |r| must be <= 1")
    w, V = np.linalg.eigh(R)
    if w.min() < 1e-10:
        warnings.warn(f"{name} is not positive definite (smallest eigenvalue {w.min():.3g}): repaired", RuntimeWarning)
        R = V @ np.diag(np.maximum(w, 1e-8)) @ V.T
        d = np.sqrt(np.diag(R))
        R = R / np.outer(d, d)
    return R


def _scale(entry, kind, A, B):
    """Standard deviation of the transformed variable: half the spread of T(mean +- sd) (exact for 'None', ~sd/mean for log)."""
    m, sd = float(entry["mean"]), float(entry["sd"])
    if kind == "log" and entry.get("sampling_type") == "lognormal":
        return max(sd, 1e-12)                 # lognormal sampling: sd is already the std of ln(x)
    lo, hi = m - sd, m + sd
    if kind in ("log", "log-ratio"):
        eps = 1e-9 * max(abs(m), 1.0)
        lo = max(lo, (A if A is not None else 0.0) + eps)
        if kind == "log-ratio":
            hi = min(hi, B - eps)
    s = 0.5 * abs(float(forward(hi, kind, A, B)) - float(forward(lo, kind, A, B)))
    return max(s, 1e-12)


def _group_id(name, base):
    m = re.fullmatch(re.escape(base) + r"(\d+)", name)
    return int(m.group(1)) if m else None


def apply_vgp_correlation(var_per, corr, sigma_vg=1.0, names=VGP_NAMES, rng=None):
    """Replace the independent draws of the VGP entries of ``var_per`` (``perturbate_parm`` output, keys such as
    ``SATCELL_VG0``) by correlated ones. ``corr``: matrix (all groups) or ``{group: matrix}``. Updates ``ini_perturbation``
    and ``sampling`` in place; returns ``(var_per, report)`` with the target and the sample correlation per group.
    Uses ``np.random`` (reproducible with ``np.random.seed``) unless ``rng`` (a ``np.random.Generator``) is given."""
    k = len(names)
    groups = sorted({g for n in names for key in var_per if (g := _group_id(key, n)) is not None})
    report = {}
    for g in groups:
        keys = [f"{n}{g}" for n in names]
        if any(key not in var_per for key in keys):
            logger.warning("group %s: not all of %s are perturbed: left independent", g, keys)
            continue
        R = corr.get(g, corr.get("default")) if isinstance(corr, dict) else corr
        if R is None:
            continue
        R = check_correlation(R, f"corr[{g}]")
        if R.shape != (k, k):
            raise ValueError(f"corr[{g}] must be {k} x {k} (order {names}), got {R.shape}")
        ents = [var_per[key] for key in keys]
        kinds, As, Bs = [], [], []
        for e in ents:
            tt = e.get("transf_type")
            kind = tt[0] if isinstance(tt, (list, tuple)) else tt
            kinds.append(kind)
            tb = e.get("transf_bounds") or {}
            As.append(tb.get("A", e.get("clip_min")))
            Bs.append(tb.get("B", e.get("clip_max")))
        n_ens = int(ents[0]["ensemble_size"])
        u = np.array([float(forward(e["mean"], kd, a, b)) for e, kd, a, b in zip(ents, kinds, As, Bs)])
        s = np.array([_scale(e, kd, a, b) for e, kd, a, b in zip(ents, kinds, As, Bs)])
        L = np.linalg.cholesky(np.outer(s, s) * R)
        z = (rng.standard_normal((k, n_ens)) if rng is not None else np.random.standard_normal((k, n_ens)))
        y = u[:, None] + float(sigma_vg) * (L @ z)
        x = np.vstack([inverse(y[i], kinds[i], As[i], Bs[i]) for i in range(k)])
        for i, e in enumerate(ents):
            lo, hi = e.get("clip_min"), e.get("clip_max")
            xi = np.clip(x[i], -np.inf if lo is None else lo, np.inf if hi is None else hi)
            e["ini_perturbation"] = xi
            e["sampling"] = xi.copy()
            e["correlated_with"] = [kk for kk in keys if kk != keys[i]]
        yy = np.vstack([forward(e["ini_perturbation"], kd, a, b) for e, kd, a, b in zip(ents, kinds, As, Bs)])
        report[g] = dict(names=keys, target=R, sample=np.corrcoef(yy) if n_ens > 2 else None, n=n_ens, sigma_vg=sigma_vg)
        logger.info("group %s: correlated %s (sigma_VG=%s, N=%d)", g, keys, sigma_vg, n_ens)
    return var_per, report


CARSEL_CODE = {"silt_loam": "SIL", "silt": "SI", "clay": "C", "sand": "S", "clay_loam": "CL"}   # petro.Carsel_Parrish_1988 rows


def carsel_parrish_bounds(soil_class):
    """Johnson bounds (A, B) and transform type of (theta_r, alpha, n) of a soil class from ``pyCATHY.petro``
    (Carsel & Parrish 1988): ``{'thetar_VG': (A, B, kind), 'alpha': (A, B, kind), 'n_VG': (A, B, kind)}`` or ``None`` when the
    class is not tabulated there (e.g. sandy loam). alpha is in 1/cm in the table."""
    code = CARSEL_CODE.get(soil_class)
    if code is None:
        return None
    from pyCATHY import petro
    _, tr, _ = petro.Carsel_Parrish_1988(code)
    kmap = {"SB": "log-ratio", "LN": "log"}
    out = {}
    for var, key in (("thetar_VG", "theta_r"), ("alpha", "alpha"), ("n_VG", "N")):
        r = tr.loc[key]
        out[var] = (float(r["A"]), float(r["B"]), kmap[r["transf_type"]])
    return out


def botto_covariance(soil_class):
    """Covariance (Sigma = T'T) and correlation of the transformed (theta_r, alpha, n) of a Botto (2018) soil class."""
    T = BOTTO_2018_T2[soil_class]
    cov = T.T @ T
    sd = np.sqrt(np.diag(cov))
    return cov, cov / np.outer(sd, sd)


def apply_botto_vgp(var_per, soil_class, sigma_vg, rng=None, groups=None, bounds="carsel"):
    """Draw the VGP of ``var_per`` as in Botto et al. (2018): ``y = u + sigma_vg T' z`` in the transformed space.

    ``soil_class``: 'sandy_loam' | 'silt_loam' (all groups) or ``{group: class}`` (group = layer / material index, e.g.
    ``{0: 'sandy_loam', 1: 'silt_loam', 2: 'silt_loam'}``). Entries ``thetar_VG<g>``, ``SATCELL_VG<g>``, ``n_VG<g>`` are
    replaced; ``SATCELL_VG`` is 1/alpha. Bounds (A, B of SB and the final clip) are the scenario ``clip_min`` / ``clip_max``.
    ``bounds='carsel'``: A, B of theta_r and n from ``petro.Carsel_Parrish_1988`` when the class is tabulated there (silt loam;
    alpha in 1/cm -> 1/m), else, or with ``'scenario'``, the scenario clips. Returns ``(var_per, report)``."""
    names = ("thetar_VG", "SATCELL_VG", "n_VG")
    found = sorted({g for n in names for key in var_per if (g := _group_id(key, n)) is not None})
    report = {}
    for g in found:
        cls = soil_class.get(g, soil_class.get("default")) if isinstance(soil_class, dict) else soil_class
        if cls is None or (groups is not None and g not in groups):
            continue
        keys = [f"{n}{g}" for n in names]
        if any(k not in var_per for k in keys):
            logger.warning("group %s: not all of %s are perturbed: left independent", g, keys)
            continue
        T = BOTTO_2018_T2[cls]
        kinds = BOTTO_2018_TRANSFORMS[cls]                    # in BOTTO_ORDER: theta_r, alpha, n
        e_thr, e_sat, e_n = (var_per[k] for k in keys)
        n_ens = int(e_thr["ensemble_size"])
        # physical bounds and nominal in Botto's variables (alpha = 1 / SATCELL_VG)
        lo_s, hi_s = e_sat.get("clip_min"), e_sat.get("clip_max")
        a_lo = 1.0 / hi_s if hi_s else None
        a_hi = 1.0 / lo_s if lo_s else None
        nom = (float(e_thr["mean"]), 1.0 / float(e_sat["mean"]), float(e_n["mean"]))
        A = (e_thr.get("clip_min"), a_lo, e_n.get("clip_min"))
        B = (e_thr.get("clip_max"), a_hi, e_n.get("clip_max"))
        cp = carsel_parrish_bounds(cls) if bounds == "carsel" else None
        if cp is not None:
            if tuple(c[2] for c in cp.values()) != tuple(kinds):
                logger.warning("group %s: petro transforms %s differ from Botto %s", g, [c[2] for c in cp.values()], kinds)
            A = (cp["thetar_VG"][0], cp["alpha"][0] * 100.0 or a_lo, cp["n_VG"][0])
            B = (cp["thetar_VG"][1], cp["alpha"][1] * 100.0, cp["n_VG"][1])
        elif bounds == "carsel":
            logger.warning("group %s: %s not in petro.Carsel_Parrish_1988: scenario bounds used", g, cls)
        u = np.array([float(forward(nom[i], kinds[i], A[i], B[i])) for i in range(3)])
        z = rng.standard_normal((3, n_ens)) if rng is not None else np.random.standard_normal((3, n_ens))
        y = u[:, None] + float(sigma_vg) * (T.T @ z)
        x = [inverse(y[i], kinds[i], A[i], B[i]) for i in range(3)]
        thr, alpha, n = x
        sat = 1.0 / alpha
        for e, v in ((e_thr, thr), (e_sat, sat), (e_n, n)):
            lo, hi = e.get("clip_min"), e.get("clip_max")
            v = np.clip(v, -np.inf if lo is None else lo, np.inf if hi is None else hi)
            e["ini_perturbation"] = v
            e["sampling"] = v.copy()
            e["correlated_with"] = [k for k in keys if k != e["type_parm"]]
        yy = np.vstack([forward(a, kinds[i], A[i], B[i]) for i, a in enumerate((e_thr["ini_perturbation"],
                                                                            1.0 / e_sat["ini_perturbation"],
                                                                            e_n["ini_perturbation"]))])
        report[g] = dict(names=keys, soil_class=cls, target=botto_covariance(cls)[1],
                         sample=np.corrcoef(yy) if n_ens > 2 else None, n=n_ens, sigma_vg=sigma_vg,
                         order="(theta_r, alpha, n) in Botto's transformed space")
        logger.info("group %s: Botto 2018 %s, sigma_VG=%s, N=%d", g, cls, sigma_vg, n_ens)
    return var_per, report
