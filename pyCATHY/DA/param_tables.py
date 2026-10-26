"""Helpers on the soil / Feddes parameter tables and on the validity of updated parameters.
Moved out of ``cathy_DA`` without change; every name is still importable from
``pyCATHY.DA.cathy_DA``."""
from __future__ import annotations

from typing import Any, Callable, List, Optional, Tuple
import logging
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Names of the perturbed variables, by family (the digits of the zone / layer
# are stripped first: "Ks0" -> "Ks", "ZROOT1" -> "ZROOT").
_VG_PARM_NAMES = ("n_VG", "thetar_VG", "alpha_VG", "SATCELL_VG")


_ARCHIE_PARM_NAMES = (
    "porosity",
    "rFluid_Archie",
    "a_Archie",
    "m_Archie",
    "n_Archie",
    "pert_sigma_Archie",
)


_FEDDES_PARM_NAMES = ("PCANA", "PCREF", "PCWLT", "ZROOT", "PZ", "OMGC")


# Columns of the soil table (`SPP_map`) set by each perturbed soil parameter.
_SOIL_PARM_COLUMNS = {
    "Ks": ["PERMX", "PERMY", "PERMZ"],
    "porosity": ["POROS"],
    "n_VG": ["VGNCELL"],
    "thetar_VG": ["VGRMCCELL"],
    "SATCELL_VG": ["VGPSATCELL"],
}


def _assign_soil_value(
    df: pd.DataFrame,
    pert_control_name: str,
    index: Optional[int],
    update_cols: List[str],
    value: Any,
    veg_map: Optional[np.ndarray] = None,
    zone: Optional[np.ndarray] = None,
    layer_groups: Optional[List[List[int]]] = None,
) -> None:
    """Set ``value`` in the columns ``update_cols`` of the soil table ``df`` (in place).

    ``df`` is indexed by (zone, layer). Where the value goes depends on how the
    perturbation was built: ``"zone"`` (zone number ``index + 1``), ``"layers"``
    (layer ``index``), ``"root_map"`` (the zones of vegetation ``index + 1``,
    from ``veg_map`` and ``zone``), ``"material"`` (every layer of ``layer_groups[index]``, in
    ALL zones), ``"zone_material"`` (``layer_groups[index] = [zone, [layers]]``: those layers of
    that zone ONLY, so a compacted zone keeps its own parameters); anything else, or no
    ``index``: every row.
    """
    if pert_control_name == "zone" and index is not None:
        df.loc[(index + 1, slice(None)), update_cols] = value
    elif pert_control_name == "layers" and index is not None:
        df.loc[(slice(None), index), update_cols] = value
    elif pert_control_name == "material" and index is not None:
        if layer_groups is None:
            raise ValueError("the 'material' control needs layer_groups (scenario 'material_layers')")
        df.loc[(slice(None), list(layer_groups[index])), update_cols] = value
    elif pert_control_name == "zone_material" and index is not None:
        if layer_groups is None:
            raise ValueError("the 'zone_material' control needs layer_groups (scenario 'zone_material_layers')")
        zone_id, layers = layer_groups[index]
        df.loc[(int(zone_id), list(layers)), update_cols] = value
    elif pert_control_name == "root_map" and index is not None:
        zones_in_vegi = zone[index + 1 == veg_map]
        df.loc[zones_in_vegi, update_cols] = value
    else:
        df.loc[:, update_cols] = value

    # Ks "c-lite": keep the anisotropy of the base table (PERMZ = Ks * base kz/kx of each row)
    ratio = df.attrs.get("kz_ratio")
    if ratio is not None and "PERMZ" in update_cols and "PERMX" in update_cols:
        mask = pd.DataFrame(False, index=df.index, columns=["_m"])
        _assign_soil_value(mask, pert_control_name, index, ["_m"], True, veg_map, zone, layer_groups)
        sel = mask["_m"].to_numpy()
        df.loc[sel, "PERMZ"] = np.asarray(df.loc[sel, "PERMX"], dtype=float) * np.asarray(
            ratio, dtype=float
        )[sel]


def _reflect_into_bounds(values, lo=None, hi=None) -> np.ndarray:
    """Fold ``values`` back into ``[lo, hi]`` by reflection (NaN untouched).

    Unlike ``np.clip`` it does not send every out-of-bounds member to the same value, so the
    ensemble keeps its spread. With one bound, the value is mirrored about it.
    """
    v = np.array(values, dtype=float)
    with np.errstate(invalid="ignore"):
        if lo is not None and hi is not None:
            w = hi - lo
            if w <= 0:
                return np.full_like(v, lo)
            r = np.mod(v - lo, 2 * w)
            return lo + np.where(r > w, 2 * w - r, r)
        if lo is not None:
            return np.where(v < lo, 2 * lo - v, v)
        if hi is not None:
            return np.where(v > hi, 2 * hi - v, v)
    return v


def _has_name_token(key: str, token: str) -> bool:
    """True if `token` appears in the (casefolded) parameter name `key` at the start
    of the name or after a non-alphanumeric character.

    Keeps the suffixes added per zone / layer ("Ks1", "porosity_zone2") but no longer
    matches inside another word: "a_Archie" used to match "pert_sigma_Archie".
    """
    return re.search(r"(?<![a-z0-9])" + re.escape(token.casefold()), key) is not None


def _nonphysical_values(
    name: str, values: np.ndarray, root_depth_limit: Callable[[], float]
) -> Optional[Tuple[np.ndarray, str]]:
    """Test the updated values of one parameter against its physical bounds.

    Returns ``(mask, label)``, ``mask`` being True for the members whose value
    is not physical, or ``None`` when no bound is known for the parameter. The
    comparisons are strict (a value on a bound is accepted) and NaN is not
    flagged. ``root_depth_limit`` is only called for a root depth.
    """
    key = name.casefold()
    if _has_name_token(key, "rFluid"):
        return values < 0, "negative"
    if _has_name_token(key, "porosity"):  # a volume fraction
        return (values < 0) | (values > 1), "out of range"
    if _has_name_token(key, "a_Archie"):
        return (values < 0) | (values > 3), "out of range"
    if _has_name_token(key, "Ks"):
        return values < 0, "negative Ks"
    if _has_name_token(key, "PCREF"):  # a pressure head: must be negative
        return values > 0, "positive PCREF"
    if _has_name_token(key, "Zroot"):
        # a root depth is a length: negative values are as unfeasible as too deep ones
        return (values < 0) | (values > root_depth_limit()), "unfeasible root depth"
    return None


def _set_feddes_value(feddes, column, row_position, value):
    """Set one value of a Feddes table in place, by column name and row position.

    ``feddes[column].iloc[row] = value`` is a chained assignment: with pandas >= 3
    (Copy-on-Write always on) it modifies a temporary and the table is left unchanged,
    so a parameter update was silently lost.
    """
    feddes.iloc[row_position, feddes.columns.get_loc(column)] = value


def _valid_zone_values(zone_map):
    """Sorted zone values of a localisation map, nodata (-9999) excluded."""
    values = np.unique(zone_map)
    return values[values != -9999]


def _check_positive_pressure_heads(ensemble_psi, label="forecast"):
    """Warn when members have a positive pressure head; returns ``{member: share of nodes}``.

    ``label`` says which ensemble was checked: the forecast (before the analysis) or the
    analysed one, whose positive heads are the effect of the update itself. Members are
    numbered from 1 and NaN columns (invalid members) are ignored.
    """
    psi = np.asarray(ensemble_psi, dtype=float)
    with np.errstate(invalid="ignore"):
        positive = psi > 0
    if not positive.any():
        return {}
    share = positive.mean(axis=0)
    members = {int(i) + 1: float(share[i]) for i in np.flatnonzero(share > 0)}
    logger.warning(
        "positive pressure heads in the %s ensemble: members %s (max %.3g, share of nodes %s)",
        label,
        sorted(members),
        float(np.nanmax(psi)),
        {m: f"{v:.0%}" for m, v in members.items()},
    )
    return members


def _update_key(cycle_nb: int) -> str:
    """Key of `dict_parm_pert[parameter]` holding the ensemble of a DA cycle.

    Cycle 0 (and below) is the initial perturbation, cycle n > 0 the n-th update.
    """
    return "ini_perturbation" if cycle_nb <= 0 else f"update_nb{cycle_nb}"


# ---------------------------------------------------------------------------
# Realised prior / posterior of the updated parameters (what the ensemble really is, not what was asked for)
# ---------------------------------------------------------------------------
def _ensembles_of(entry: dict) -> List[np.ndarray]:
    """Prior then every analysed ensemble of one ``dict_parm_pert`` entry, each flattened to (n_member,)."""
    rows = [np.ravel(np.asarray(entry["ini_perturbation"], dtype=float))]
    keys = sorted(int(k[len("update_nb"):]) for k in entry if str(k).startswith("update_nb") and str(k)[9:].isdigit())
    rows += [np.ravel(np.asarray(entry[f"update_nb{k}"], dtype=float)) for k in keys]
    return rows


def parameter_summary(dict_parm_pert: dict, bound_tol: float = 1e-6) -> dict:
    """Per parameter: sampling law, bounds, REALISED prior mean / median / sd, final mean / sd, contraction
    ``1 - sd_final / sd_prior`` and the fraction of the final ensemble sitting on a bound.

    Why: the requested prior is not the realised one. ``lognormal`` takes the nominal value as the MEDIAN unless ``lognormal_center="mean"`` is set (the mean is
    then higher by exp(sd^2 / 2): +65 % at sd = 1), and a truncated normal has a smaller sd than requested (0.30 -> 0.16 with
    bounds at +/- 0.67 sd). Entries without an ensemble are skipped.
    """
    out = {}
    for name, entry in (dict_parm_pert or {}).items():
        if not isinstance(entry, dict) or "ini_perturbation" not in entry:
            continue
        ens = _ensembles_of(entry)
        first, last = ens[0], ens[-1]
        nominal = entry.get("nominal")
        nominal = float(nominal) if isinstance(nominal, (int, float, np.integer, np.floating)) else None
        sd0, sd1 = float(np.nanstd(first, ddof=1)), float(np.nanstd(last, ddof=1))
        lo, hi = entry.get("clip_min"), entry.get("clip_max")
        scale = max(abs(float(np.nanmean(last))), abs(float(np.nanstd(last))), 1e-300)
        at_bound = np.zeros(last.shape, bool)
        for b in (lo, hi):
            if b is not None:
                at_bound |= np.abs(last - float(b)) <= bound_tol * max(abs(float(b)), scale)
        out[name] = dict(
            sampling_type=entry.get("sampling_type"), per_type=entry.get("per_type"), nominal=nominal,
            clip_min=lo, clip_max=hi, n_members=int(np.sum(~np.isnan(first))), n_updates=len(ens) - 1,
            prior_mean=float(np.nanmean(first)), prior_median=float(np.nanmedian(first)), prior_sd=sd0,
            prior_mean_over_nominal=(float(np.nanmean(first)) / nominal) if nominal else None,
            final_mean=float(np.nanmean(last)), final_sd=sd1,
            contraction=(1.0 - sd1 / sd0) if sd0 > 0 else None,
            frac_at_bound_final=float(np.mean(at_bound[~np.isnan(last)])) if (lo is not None or hi is not None) else None,
        )
    return out


def parameter_flags(summary: dict, collapse: float = 0.95, at_bound: float = 0.25) -> List[Tuple[str, str]]:
    """``[(parameter, message)]`` for the situations that make an updated parameter ensemble useless:
    spread contracted by ``collapse`` or more (no further information can be assimilated), or ``at_bound`` of the members
    piled on a bound (the bound, not the data, decides the value)."""
    out = []
    for name, s in summary.items():
        if s["n_updates"] and s["contraction"] is not None and s["contraction"] >= collapse:
            out.append((name, f"ensemble spread contracted by {100 * s['contraction']:.0f} % "
                              f"(sd {s['prior_sd']:.3g} -> {s['final_sd']:.3g}): collapsed, check inflation / perturbed obs"))
        if s["frac_at_bound_final"] is not None and s["frac_at_bound_final"] >= at_bound:
            out.append((name, f"{100 * s['frac_at_bound_final']:.0f} % of the final members sit on a bound "
                              f"[{s['clip_min']}, {s['clip_max']}]: the bound decides, not the data"))
    return out
