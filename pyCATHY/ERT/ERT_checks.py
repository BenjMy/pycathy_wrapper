"""Sanity checks for ERT observations, forward-modelled data and ensemble predictions.

Every ``check_*`` function returns an :class:`ERTCheckReport`; nothing is raised.
Call ``report.warn()`` to log one warning per failed check (logger of this
module; ``EMIT_PYTHON_WARNINGS = True`` adds ``ERTCheckWarning``), and
``report.to_dict()`` to get a JSON-friendly record for the pipeline report.

Checked, by stage:

* scheme / observations : electrode indices, geometric factor ``k``, ``rhoa``,
  ``rhoa == k * r`` consistency, ``valid`` flags, error model (``rec_err``)
* single forward run    : ``k`` outliers, non-finite / non-positive / implausible
  ``rhoa`` relative to the input resistivity model, signal below noise floor
* forward vs observed   : same quadrupoles (indices, off-by-one), same ``k``
* ensemble of forwards  : same geometry in all members, duplicated members,
  member regime changes, systematic misfit, ensemble collapse / rank-1 spread
* analysis inputs       : data covariance, normalised innovation (and the
  log-transform / covariance units mismatch)
"""

from __future__ import annotations

import json
import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

__all__ = [
    "ERTCheckWarning",
    "ERTCheckConfig",
    "ERTCheckReport",
    "check_scheme",
    "check_geometric_factor",
    "check_rhoa",
    "check_forward_data",
    "check_forward_vs_obs",
    "check_observations",
    "check_data_cov",
    "check_archie_inputs",
    "check_ensemble_forward",
    "check_ensemble_predictions",
    "check_innovation",
    "check_assimilation_inputs",
    "save_reports",
    "count_failed",
    "issues_summary",
]


logger = logging.getLogger(__name__)

# Failed checks are logged (logger.warning). Set to True to also emit them as
# Python `ERTCheckWarning` warnings (e.g. to turn them into errors in tests).
EMIT_PYTHON_WARNINGS = False


class ERTCheckWarning(UserWarning):
    """A sanity check on ERT data failed."""


def issues_summary(issues) -> str:
    """'OK' or 'N issue(s): code1, code2, ...' for a list of issue dicts."""
    if not issues:
        return "OK"
    codes = sorted({it["code"] for it in issues})
    return f"{len(issues)} issue(s): " + ", ".join(codes)


@dataclass
class ERTCheckConfig:
    """Thresholds of the checks (defaults are deliberately permissive)."""

    k_mad_factor: float = 8.0  # |log10|k| - median| > factor * spread -> outlier
    k_min_spread_decades: float = 0.5  # floor of the spread (decades) in that test
    k_abs_max: float = 1e5  # hard bound on |k| [m]
    k_ratio_tol: float = 0.05  # |k_fwd / k_obs - 1| tolerance
    rhoa_factor: float = 10.0  # plausible rhoa in [rho_min/f, rho_max*f]
    max_invalid_frac: float = 0.05  # share of valid==0 above which we warn
    snr_min: float = 10.0  # |R| should exceed snr_min * noise_abs
    rel_err_range: tuple = (1e-4, 1.0)  # plausible sqrt(var)/|data|
    dynamic_range_max: float = 1e4  # max/min of positive observed rhoa
    member_regime_factor: float = 5.0  # member max|Hx| vs median over members
    misfit_decades: float = 2.0  # systematic misfit: |log10(Hx/obs)| > this ...
    misfit_member_frac: float = 0.9  # ... in at least this share of members
    collapse_rel_spread: float = 1e-6  # median ensemble sd / |mean| below this
    rank1_var_frac: float = 0.99  # first singular mode share of variance
    nis_max: float = 3.0  # rms normalised innovation above this -> warn
    log_eps: float = 1e-3  # offset used by the log transform of the DA
    max_listed_idx: int = 20


DEFAULT_CONFIG = ERTCheckConfig()


@dataclass
class ERTCheckReport:
    """Outcome of a group of checks."""

    name: str
    issues: List[Dict[str, Any]] = field(default_factory=list)
    stats: Dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.issues

    def add(
        self,
        code: str,
        message: str,
        mask: Optional[np.ndarray] = None,
        count: Optional[int] = None,
        max_idx: int = DEFAULT_CONFIG.max_listed_idx,
    ) -> None:
        idx: List[int] = []
        if mask is not None:
            mask = np.asarray(mask, dtype=bool)
            count = int(mask.sum())
            idx = np.flatnonzero(mask)[:max_idx].tolist()
        self.issues.append(
            {"code": code, "message": message, "count": count, "idx": idx}
        )

    def merge(self, other: "ERTCheckReport") -> "ERTCheckReport":
        self.issues.extend(other.issues)
        self.stats.update(other.stats)
        return self

    def summary(self) -> str:
        return issues_summary(self.issues)

    def warn(self, stacklevel: int = 2) -> "ERTCheckReport":
        """Log every failed check (logger.warning); also a Python warning if enabled."""
        for it in self.issues:
            extra = ""
            if it["count"]:
                extra = f" (n={it['count']}"
                extra += f", first idx={it['idx']})" if it["idx"] else ")"
            msg = f"[{self.name}] {it['code']}: {it['message']}{extra}"
            logger.warning(msg)
            if EMIT_PYTHON_WARNINGS:
                warnings.warn(msg, ERTCheckWarning, stacklevel=stacklevel)
        return self

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "ok": self.ok,
            "issues": self.issues,
            "stats": self.stats,
        }


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _get(container: Any, *keys: str) -> Optional[np.ndarray]:
    """First available key of a pygimli container / DataFrame / dict as a float vector."""
    for key in keys:
        try:
            if hasattr(container, "haveData") and not container.haveData(key):
                continue
            return np.asarray(container[key], dtype=float).ravel()
        except Exception:  # missing key, wrong container: try the next one
            continue
    return None


def _stat(x: np.ndarray, fn) -> Optional[float]:
    x = x[np.isfinite(x)]
    return float(fn(x)) if x.size else None


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


# ---------------------------------------------------------------------------
# scheme and geometric factor
# ---------------------------------------------------------------------------


def check_scheme(a, b, m, n, n_elec: Optional[int] = None, name: str = "scheme"):
    """Electrode indices of the quadrupoles (0-based, -1 = electrode at infinity)."""
    rep = ERTCheckReport(name)
    arrs = [np.asarray(v, dtype=float).ravel() for v in (a, b, m, n)]
    if len({x.size for x in arrs}) != 1:
        rep.add("scheme_length", "a, b, m, n have different lengths")
        return rep
    allv = np.vstack(arrs)
    rep.stats["n_quadrupoles"] = int(allv.shape[1])
    if allv.shape[1] == 0:
        rep.add("scheme_empty", "no quadrupole in the scheme")
        return rep

    fin = np.isfinite(allv).all(axis=0)
    if (~fin).any():
        rep.add("scheme_nan", "non-finite electrode index", mask=~fin)
    ia, ib, im, in_ = np.rint(np.where(np.isfinite(allv), allv, 0)).astype(int)
    iv = np.vstack([ia, ib, im, in_])

    neg = fin & (iv < 0).any(axis=0)
    if neg.any():
        rep.add(
            "scheme_negative_index",
            "negative electrode index (-1 is legitimate only for pole configurations; "
            "otherwise a 0/1-based indexing shift, e.g. the `-1` applied to .txt sequences)",
            mask=neg,
        )
    if n_elec is not None:
        oor = fin & (iv >= n_elec).any(axis=0)
        if oor.any():
            rep.add(
                "scheme_index_out_of_range",
                f"electrode index >= number of electrodes ({n_elec})",
                mask=oor,
            )
    same_ab = fin & (ia == ib) & (ia >= 0)
    same_mn = fin & (im == in_) & (im >= 0)
    if same_ab.any():
        rep.add("scheme_a_equals_b", "current electrodes A and B coincide", mask=same_ab)
    if same_mn.any():
        rep.add("scheme_m_equals_n", "potential electrodes M and N coincide", mask=same_mn)
    shared = fin & (
        ((ia == im) & (ia >= 0))
        | ((ia == in_) & (ia >= 0))
        | ((ib == im) & (ib >= 0))
        | ((ib == in_) & (ib >= 0))
    )
    if shared.any():
        rep.add(
            "scheme_shared_electrode",
            "a current electrode is also a potential electrode (k diverges or is ill-defined)",
            mask=shared,
        )
    n_dup = allv.shape[1] - np.unique(iv.T, axis=0).shape[0]
    if n_dup > 0:
        rep.add("scheme_duplicates", "duplicated quadrupoles", count=int(n_dup))
    return rep


def check_geometric_factor(k, cfg: Optional[ERTCheckConfig] = None, name: str = "k"):
    """Geometric factor: finite, non-zero, no gross outlier (decades from the median)."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    k = None if k is None else np.asarray(k, dtype=float).ravel()
    if k is None or k.size == 0:
        rep.add("k_missing", "geometric factor not available")
        return rep

    fin = np.isfinite(k)
    if (~fin).any():
        rep.add("k_nonfinite", "NaN/inf geometric factor", mask=~fin)
    zero = fin & (k == 0)
    if zero.any():
        rep.add(
            "k_zero",
            "k == 0 (geometric factors not computed: rhoa is then 0 or meaningless)",
            mask=zero,
        )
    good = fin & (k != 0)
    if good.sum() >= 3:
        lk = np.log10(np.abs(k[good]))
        med = np.median(lk)
        spread = max(1.4826 * np.median(np.abs(lk - med)), cfg.k_min_spread_decades)
        out = np.zeros(k.size, dtype=bool)
        out[good] = np.abs(lk - med) > cfg.k_mad_factor * spread
        if out.any():
            rep.add(
                "k_outlier",
                f"|k| more than {cfg.k_mad_factor * spread:.1f} decades from the median "
                f"(median |k| = {10 ** med:.3g}): near-singular quadrupoles blow up rhoa",
                mask=out,
            )
    big = fin & (np.abs(k) > cfg.k_abs_max)
    if big.any():
        rep.add("k_too_large", f"|k| > {cfg.k_abs_max:g} m", mask=big)
    rep.stats.update(
        k_min=_stat(k, np.min),
        k_max=_stat(k, np.max),
        k_median_abs=_stat(np.abs(k), np.median),
        k_frac_negative=float(np.mean(k[fin] < 0)) if fin.any() else None,
    )
    return rep


# ---------------------------------------------------------------------------
# apparent resistivity
# ---------------------------------------------------------------------------


def check_rhoa(
    rhoa,
    res0=None,
    cfg: Optional[ERTCheckConfig] = None,
    name: str = "rhoa",
):
    """rhoa finite, positive and, if the input model ``res0`` is known, plausible."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    rhoa = None if rhoa is None else np.asarray(rhoa, dtype=float).ravel()
    if rhoa is None or rhoa.size == 0:
        rep.add("rhoa_missing", "apparent resistivity not available")
        return rep

    fin = np.isfinite(rhoa)
    if (~fin).any():
        rep.add("rhoa_nonfinite", "NaN/inf apparent resistivity", mask=~fin)
    nonpos = fin & (rhoa <= 0)
    if nonpos.any():
        rep.add(
            "rhoa_nonpositive",
            "rhoa <= 0 (they are later forced to 1e-3, which creates new outliers "
            "and breaks log transforms)",
            mask=nonpos,
        )
    if res0 is not None:
        r0 = np.asarray(res0, dtype=float).ravel()
        r0 = r0[np.isfinite(r0) & (r0 > 0)]
        if r0.size:
            lo, hi = r0.min() / cfg.rhoa_factor, r0.max() * cfg.rhoa_factor
            implausible = fin & ((np.abs(rhoa) > hi) | ((rhoa > 0) & (rhoa < lo)))
            if implausible.any():
                rep.add(
                    "rhoa_implausible",
                    f"|rhoa| outside [{lo:.3g}, {hi:.3g}] Ohm.m although the model "
                    f"spans {r0.min():.3g}-{r0.max():.3g} Ohm.m "
                    "(near-zero potentials / noise-dominated or singular quadrupoles)",
                    mask=implausible,
                )
            rep.stats.update(res0_min=float(r0.min()), res0_max=float(r0.max()))
    rep.stats.update(
        rhoa_min=_stat(rhoa, np.min),
        rhoa_max=_stat(rhoa, np.max),
        rhoa_median=_stat(rhoa, np.median),
    )
    return rep


def check_forward_data(
    het,
    res0=None,
    noise_abs: Optional[float] = None,
    cfg: Optional[ERTCheckConfig] = None,
    name: str = "fwd ERT",
):
    """All checks on the output of one forward run (pygimli container or DataFrame).

    ``noise_abs`` is the absolute voltage noise of the simulation [V]; when given,
    quadrupoles whose resistance ``R = rhoa / k`` is below ``snr_min * noise_abs``
    are reported (their rhoa is dominated by noise, then multiplied by ``k``).
    """
    cfg = cfg or DEFAULT_CONFIG
    k = _get(het, "k")
    rhoa = _get(het, "rhoa")
    rep = check_geometric_factor(k, cfg, name=name)
    rep.merge(check_rhoa(rhoa, res0, cfg, name=name))

    valid = _get(het, "valid")
    if valid is not None:
        invalid = valid == 0
        if invalid.mean() > cfg.max_invalid_frac:
            rep.add(
                "valid_fraction",
                f"more than {100 * cfg.max_invalid_frac:.0f}% of the data flagged invalid",
                mask=invalid,
            )
        elif invalid.any():
            rep.add("valid_flag", "data flagged invalid by the simulator", mask=invalid)

    if noise_abs is not None and k is not None and rhoa is not None and k.size == rhoa.size:
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(k != 0, rhoa / k, np.nan)
        weak = np.isfinite(r) & (np.abs(r) < cfg.snr_min * noise_abs)
        rep.stats["noise_abs_assumed"] = float(noise_abs)
        if weak.any():
            rep.add(
                "signal_below_noise",
                f"|R| = |rhoa/k| below {cfg.snr_min:g} x assumed absolute noise "
                f"({noise_abs:g} V): rhoa is noise-dominated",
                mask=weak,
            )
    return rep


def check_forward_vs_obs(
    pred, obs, cfg: Optional[ERTCheckConfig] = None, name: str = "fwd vs obs"
):
    """Forward data must use the same quadrupoles (and k) as the observations."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    if obs is None or pred is None:
        return rep
    pa = [_get(pred, key) for key in "abmn"]
    oa = [_get(obs, key) for key in "abmn"]
    if any(v is None for v in pa + oa):
        return rep  # observations without electrode indices (e.g. resipy csv)
    if pa[0].size != oa[0].size:
        rep.add(
            "fwd_obs_size",
            f"{pa[0].size} forward data vs {oa[0].size} observations: "
            "the prediction cannot be matched to the data",
        )
        return rep
    differ = np.zeros(pa[0].size, dtype=bool)
    shifts = set()
    for p, o in zip(pa, oa):
        differ |= p != o
        shifts.update(np.unique((p - o)[p != o]).tolist())
    if differ.any():
        hint = (
            f" (constant index shift {sorted(shifts)}: 0/1-based mismatch?)"
            if len(shifts) == 1
            else ""
        )
        rep.add(
            "fwd_obs_quadrupoles",
            "forward and observed quadrupoles differ" + hint,
            mask=differ,
        )
    kp, ko = _get(pred, "k"), _get(obs, "k", "K")
    if kp is not None and ko is not None and kp.size == ko.size:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = kp / ko
        ok = np.isfinite(ratio)
        bad = ~ok | (np.abs(ratio - 1) > cfg.k_ratio_tol)
        if bad.any():
            rep.add(
                "fwd_obs_k",
                f"forward k differs from observed k by more than {100 * cfg.k_ratio_tol:.0f}% "
                "(different topography / numerical vs analytical k?)",
                mask=bad,
            )
        rep.stats["k_ratio_median"] = _stat(ratio, np.median)
    return rep


# ---------------------------------------------------------------------------
# observations and error model
# ---------------------------------------------------------------------------


def quadrupole_quality_mask(
    df,
    cfg: Optional[ERTCheckConfig] = None,
    res0=None,
    noise_abs: Optional[float] = None,
):
    """Boolean mask (True = keep) of the quadrupoles usable for assimilation.

    Applies the same criteria as the checks, but returns the *full* mask instead of
    the first few indices: non-finite or non-positive rhoa, ``valid == 0``, and,
    when ``k`` is available, ``k`` non-finite / zero / far from the median (decades)
    / larger than ``k_abs_max`` (these near-singular quadrupoles blow up rhoa).
    ``res0`` (optional model resistivities) also drops implausible rhoa, and
    ``noise_abs`` (V, optional) drops ``|rhoa/k| < snr_min * noise_abs``.
    Returns ``(keep, reasons)`` with ``reasons`` a dict code -> number of data.
    """
    cfg = cfg or DEFAULT_CONFIG
    rhoa = _get(df, "rhoa")
    if rhoa is None:
        return np.array([], dtype=bool), {}
    rhoa = np.asarray(rhoa, dtype=float).ravel()
    bad = {}
    bad["rhoa_invalid"] = ~np.isfinite(rhoa) | (rhoa <= 0)

    valid = _get(df, "valid")
    if valid is not None and np.size(valid) == rhoa.size:
        bad["flagged_invalid"] = np.asarray(valid).ravel() == 0

    k = _get(df, "k", "K")
    if k is not None and np.size(k) == rhoa.size:
        k = np.asarray(k, dtype=float).ravel()
        fin = np.isfinite(k)
        bad["k_invalid"] = ~fin | (k == 0)
        good = fin & (k != 0)
        out = np.zeros(k.size, dtype=bool)
        if good.sum() >= 3:
            lk = np.log10(np.abs(k[good]))
            med = np.median(lk)
            spread = max(1.4826 * np.median(np.abs(lk - med)), cfg.k_min_spread_decades)
            out[good] = np.abs(lk - med) > cfg.k_mad_factor * spread
        bad["k_outlier"] = out
        bad["k_too_large"] = fin & (np.abs(k) > cfg.k_abs_max)
        if noise_abs is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(good, rhoa / np.where(k == 0, np.nan, k), np.nan)
            bad["signal_below_noise"] = np.isfinite(r) & (np.abs(r) < cfg.snr_min * noise_abs)

    if res0 is not None:
        r0 = np.asarray(res0, dtype=float).ravel()
        r0 = r0[np.isfinite(r0) & (r0 > 0)]
        if r0.size:
            lo, hi = r0.min() / cfg.rhoa_factor, r0.max() * cfg.rhoa_factor
            bad["rhoa_implausible"] = np.isfinite(rhoa) & ((np.abs(rhoa) > hi) | ((rhoa > 0) & (rhoa < lo)))

    drop = np.zeros(rhoa.size, dtype=bool)
    for m in bad.values():
        drop |= m
    reasons = {code: int(m.sum()) for code, m in bad.items() if m.any()}
    return ~drop, reasons


def check_observations(df, cfg: Optional[ERTCheckConfig] = None, name: str = "ERT obs"):
    """Observed ERT data container: rhoa, k, r, error model, validity, scheme."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    rhoa = _get(df, "rhoa")
    k = _get(df, "k", "K")
    r = _get(df, "r", "resist")
    if rhoa is None and k is not None and r is not None:
        rhoa = k * r
    if rhoa is None:
        rep.add("obs_no_rhoa", "no apparent resistivity (rhoa, or k and r) in the observations")
        return rep

    rep.merge(check_rhoa(rhoa, None, cfg, name=name))
    pos = rhoa[np.isfinite(rhoa) & (rhoa > 0)]
    if pos.size and pos.max() / pos.min() > cfg.dynamic_range_max:
        rep.add(
            "obs_dynamic_range",
            f"rhoa spans {pos.max() / pos.min():.1e} (> {cfg.dynamic_range_max:.0e}): "
            "outliers or unfiltered data",
        )
    if k is not None:
        rep.merge(check_geometric_factor(k, cfg, name=name))
        if r is not None and k.size == r.size == rhoa.size:
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.abs(rhoa - k * r) / np.abs(rhoa)
            inc = np.isfinite(rel) & (rel > 1e-2)
            if inc.any():
                rep.add("obs_rhoa_k_r", "rhoa != k * r (inconsistent columns)", mask=inc)

    err = _get(df, "rec_err", "err")
    if err is not None:
        bad = ~np.isfinite(err) | (err <= 0)
        if bad.any():
            rep.add(
                "obs_error_invalid",
                "reciprocal/measurement error is NaN or <= 0 "
                "(zero variance makes the data covariance singular)",
                mask=bad,
            )
        rep.stats["err_median"] = _stat(np.abs(err), np.median)

    valid = _get(df, "valid")
    if valid is not None and (valid == 0).any():
        rep.add("obs_flagged_invalid", "observations flagged invalid still present", mask=valid == 0)

    idx = [_get(df, key) for key in "abmn"]
    if all(v is not None for v in idx):
        rep.merge(check_scheme(*idx, name=name))
    return rep


def check_data_cov(
    var_diag, values=None, cfg: Optional[ERTCheckConfig] = None, name: str = "data cov"
):
    """Diagonal of the data covariance R: positive and plausible relative to the data."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    var = np.asarray(var_diag, dtype=float).ravel()
    bad = ~np.isfinite(var) | (var <= 0)
    if bad.any():
        rep.add("cov_invalid", "data variance is NaN or <= 0", mask=bad)
    if values is not None:
        vals = np.asarray(values, dtype=float).ravel()
        if vals.size == var.size:
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.sqrt(var) / np.abs(vals)
            med = _stat(rel, np.median)
            rep.stats["rel_err_median"] = med
            lo, hi = cfg.rel_err_range
            if med is not None and not (lo <= med <= hi):
                rep.add(
                    "cov_implausible",
                    f"median sqrt(variance)/|data| = {med:.3g}, outside [{lo:g}, {hi:g}] "
                    "(percent vs fraction vs variance units?)",
                )
    return rep


def check_archie_inputs(
    sat, porosity, params: Optional[Dict[str, Any]] = None, name: str = "Archie"
):
    """Saturation, porosity and Archie parameters entering the state -> resistivity map."""
    rep = ERTCheckReport(name)
    sat = np.asarray(sat, dtype=float).ravel()
    poro = np.asarray(porosity, dtype=float).ravel()
    if (~np.isfinite(sat)).any():
        rep.add("sat_nonfinite", "NaN/inf saturation", mask=~np.isfinite(sat))
    out = np.isfinite(sat) & ((sat <= 0) | (sat > 1 + 1e-6))
    if out.any():
        rep.add("sat_out_of_range", "saturation outside (0, 1]", mask=out)
    badp = ~np.isfinite(poro) | (poro <= 0) | (poro > 1)
    if badp.any():
        rep.add("porosity_out_of_range", "porosity outside (0, 1]", mask=badp)
    for key, value in (params or {}).items():
        v = np.asarray(value, dtype=float).ravel()
        if key != "pert_sigma_Archie" and ((~np.isfinite(v)) | (v <= 0)).any():
            rep.add("archie_param_invalid", f"{key} must be finite and > 0, got {v.tolist()[:5]}")
    return rep


# ---------------------------------------------------------------------------
# ensemble
# ---------------------------------------------------------------------------


def _ids(ens_ids, n):
    return list(ens_ids) if ens_ids is not None and len(ens_ids) == n else list(range(n))


def check_ensemble_predictions(
    Hx,
    obs=None,
    ens_ids: Optional[Sequence[int]] = None,
    cfg: Optional[ERTCheckConfig] = None,
    name: str = "ERT ensemble",
):
    """Predicted observations ``Hx`` (n_obs x n_members): validity, duplicates, spread."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    Hx = np.asarray(Hx, dtype=float)
    if Hx.ndim == 1:
        Hx = Hx[:, None]
    n_obs, n_ens = Hx.shape
    ids = _ids(ens_ids, n_ens)
    rep.stats.update(n_obs=int(n_obs), n_members=int(n_ens))

    nan_cols = np.isnan(Hx).any(axis=0)
    if nan_cols.any():
        rep.add("pred_nan", f"NaN predictions in members {[ids[i] for i in np.flatnonzero(nan_cols)]}",
                mask=np.isnan(Hx).any(axis=1))
    nonpos = np.nan_to_num(Hx, nan=1.0) <= 0
    if nonpos.any():
        per = nonpos.sum(axis=0)
        rep.add(
            "pred_nonpositive",
            f"non-positive predictions (max {int(per.max())} per member; members "
            f"{[ids[i] for i in np.flatnonzero(per)][:10]})",
            mask=nonpos.any(axis=1),
        )

    seen: Dict[bytes, int] = {}
    dups = []
    for j in range(n_ens):
        key = np.round(Hx[:, j], 6).tobytes()
        if key in seen:
            dups.append((ids[seen[key]], ids[j]))
        else:
            seen[key] = j
    if dups:
        rep.add(
            "member_duplicated",
            f"identical predictions in members {dups[:10]} (copied or mis-indexed member: "
            "the effective ensemble size is reduced)",
            count=len(dups),
        )

    mx = np.nanmax(np.abs(Hx), axis=0)
    if n_ens >= 3 and np.isfinite(mx).all() and np.median(mx) > 0:
        med = np.median(mx)
        reg = (mx > cfg.member_regime_factor * med) | (mx < med / cfg.member_regime_factor)
        if reg.any():
            rep.add(
                "member_regime_change",
                f"max |Hx| of members {[ids[i] for i in np.flatnonzero(reg)]} differs by more "
                f"than x{cfg.member_regime_factor:g} from the ensemble median ({med:.3g})",
                count=int(reg.sum()),
            )

    if obs is not None:
        o = np.asarray(obs, dtype=float).ravel()
        if o.size == n_obs:
            with np.errstate(divide="ignore", invalid="ignore"):
                lr = np.log10(np.clip(np.abs(Hx), 1e-300, None)) - np.log10(
                    np.clip(np.abs(o), 1e-300, None)
                )[:, None]
            frac = np.nanmean(np.abs(lr) > cfg.misfit_decades, axis=1)
            sysm = frac >= cfg.misfit_member_frac
            rep.stats["median_abs_log10_misfit"] = _stat(np.abs(lr), np.median)
            if sysm.any():
                rep.add(
                    "systematic_misfit",
                    f"prediction off by > {cfg.misfit_decades:g} decades from the observation in "
                    f">= {100 * cfg.misfit_member_frac:.0f}% of the members",
                    mask=sysm,
                )
        else:
            rep.add("obs_size", f"{o.size} observations but {n_obs} predicted values per member")

    if n_ens >= 3:
        mean = np.nanmean(Hx, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.nanstd(Hx, axis=1, ddof=1) / (np.abs(mean) + 1e-12)
        med_rel = _stat(rel, np.median)
        rep.stats["median_rel_spread"] = med_rel
        if med_rel is not None and med_rel < cfg.collapse_rel_spread:
            rep.add("ensemble_collapse", f"median ensemble spread/|mean| = {med_rel:.1e}: no spread")
        L = np.log10(np.clip(np.nan_to_num(np.abs(Hx), nan=1.0), 1e-3, None))
        A = L - L.mean(axis=1, keepdims=True)
        s = np.linalg.svd(A, compute_uv=False)
        tot = float((s**2).sum())
        if tot > 0:
            frac1 = float(s[0] ** 2 / tot)
            rep.stats["first_mode_var_frac"] = frac1
            if frac1 > cfg.rank1_var_frac:
                rep.add(
                    "rank1_ensemble",
                    f"{100 * frac1:.1f}% of the predicted-data variance is in one mode: "
                    "the update is essentially a uniform shift",
                )
    return rep


def check_ensemble_forward(
    preds: Sequence[Any],
    obs_df=None,
    ens_ids: Optional[Sequence[int]] = None,
    cfg: Optional[ERTCheckConfig] = None,
    name: str = "ERT ensemble",
):
    """Ensemble of forward results (one DataFrame per member): geometry + predictions."""
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    if len(preds) == 0:
        rep.add("ensemble_empty", "no forward result")
        return rep
    ids = _ids(ens_ids, len(preds))

    ref = {key: _get(preds[0], key) for key in ("a", "b", "m", "n", "k")}
    geo_bad, k_bad = [], []
    for j, p in enumerate(preds[1:], start=1):
        for key in "abmn":
            cur = _get(p, key)
            if cur is None or ref[key] is None or not np.array_equal(cur, ref[key]):
                geo_bad.append(ids[j])
                break
        kj = _get(p, "k")
        if kj is None or ref["k"] is None or kj.shape != ref["k"].shape or not np.allclose(
            kj, ref["k"], rtol=1e-6, equal_nan=True
        ):
            k_bad.append(ids[j])
    if geo_bad:
        rep.add("ensemble_quadrupoles_differ", f"members {geo_bad} use other quadrupoles than member {ids[0]}")
    if k_bad:
        rep.add("ensemble_k_differs", f"geometric factor of members {k_bad} differs from member {ids[0]}")

    rhoas = [_get(p, "rhoa") for p in preds]
    if all(r is not None for r in rhoas) and len({r.size for r in rhoas}) == 1:
        obs = _get(obs_df, "rhoa") if obs_df is not None else None
        rep.merge(check_ensemble_predictions(np.column_stack(rhoas), obs, ids, cfg, name))
    else:
        rep.add("ensemble_rhoa_shape", "members have missing or different-sized rhoa")
    return rep


def check_innovation(
    data,
    Hx,
    R_diag,
    use_log: bool = False,
    cfg: Optional[ERTCheckConfig] = None,
    name: str = "innovation",
):
    """Normalised innovation ``(d - mean(Hx))^2 / (spread^2 + R)`` should be O(1).

    With ``use_log`` the data and predictions are transformed as the DA does
    (``log10(x + log_eps)``) and R is checked for the units mismatch (R is not
    transformed together with the data).
    """
    cfg = cfg or DEFAULT_CONFIG
    rep = ERTCheckReport(name)
    d = np.asarray(data, dtype=float).ravel()
    H = np.asarray(Hx, dtype=float)
    R = np.asarray(R_diag, dtype=float).ravel()
    if H.ndim == 1:
        H = H[:, None]
    if not (d.size == H.shape[0] == R.size):
        rep.add("innovation_shape", f"data {d.size}, Hx {H.shape}, R {R.size}: sizes differ")
        return rep
    if use_log:
        with np.errstate(invalid="ignore", divide="ignore"):
            d = np.log10(np.clip(d, 0.0, None) + cfg.log_eps)
            H = np.log10(np.clip(H, 0.0, None) + cfg.log_eps)
        med_R = _stat(R, np.median)
        if med_R is not None and med_R > 1.0:
            rep.add(
                "R_not_log_transformed",
                f"data are assimilated in log10 but the median data variance is {med_R:.3g} "
                "(linear units): R is far too large and the data are ignored",
            )
    innov = d - np.nanmean(H, axis=1)
    spread_var = np.nanvar(H, axis=1, ddof=1) if H.shape[1] > 1 else np.zeros(d.size)
    with np.errstate(divide="ignore", invalid="ignore"):
        nis = innov**2 / (spread_var + R)
    rms = _stat(nis, lambda x: np.sqrt(np.mean(x)))
    rep.stats.update(
        innovation_mean=_stat(innov, np.mean),
        innovation_std=_stat(innov, np.std),
        rms_normalised_innovation=rms,
    )
    if rms is not None and rms > cfg.nis_max:
        rep.add(
            "innovation_large",
            f"rms normalised innovation = {rms:.3g} (expected ~1, limit {cfg.nis_max:g}): "
            "model bias, outliers in the predictions or under-estimated errors",
        )
    return rep


def check_assimilation_inputs(
    data,
    Hx,
    data_cov,
    use_log: bool = False,
    ens_ids: Optional[Sequence[int]] = None,
    cfg: Optional[ERTCheckConfig] = None,
    name: str = "ERT analysis",
):
    """Observation vector, predictions and data covariance as they enter the filter."""
    cfg = cfg or DEFAULT_CONFIG
    d = np.asarray(data, dtype=float).ravel()
    cov = np.asarray(data_cov, dtype=float)
    R = np.diag(cov) if cov.ndim == 2 else cov.ravel()

    rep = ERTCheckReport(name)
    bad = ~np.isfinite(d) | (d <= 0)
    if bad.any():
        rep.add("obs_invalid", "non-finite or non-positive observations enter the analysis", mask=bad)
    rep.merge(check_ensemble_predictions(Hx, d, ens_ids, cfg, name))
    rep.merge(check_data_cov(R, d, cfg, name))
    rep.merge(check_innovation(d, Hx, R, use_log, cfg, name))
    return rep


# ---------------------------------------------------------------------------
# bookkeeping
# ---------------------------------------------------------------------------


def count_failed(obj: Any) -> int:
    """Number of report dicts (``to_dict`` output, possibly nested) that are not ok."""
    if isinstance(obj, dict):
        n = 1 if obj.get("ok") is False else 0
        return n + sum(count_failed(v) for k, v in obj.items() if k != "issues")
    if isinstance(obj, (list, tuple)):
        return sum(count_failed(v) for v in obj)
    return 0


def save_reports(reports: Sequence[Dict[str, Any]], path: str) -> None:
    """Write the accumulated reports to a JSON file."""
    with open(path, "w") as fh:
        json.dump(list(reports), fh, indent=2, default=_json_default)
