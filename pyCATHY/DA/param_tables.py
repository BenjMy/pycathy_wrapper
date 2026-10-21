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
    from ``veg_map`` and ``zone``), ``"material"`` (every layer of ``layer_groups[index]``);
    anything else, or no ``index``: every row.
    """
    if pert_control_name == "zone" and index is not None:
        df.loc[(index + 1, slice(None)), update_cols] = value
    elif pert_control_name == "layers" and index is not None:
        df.loc[(slice(None), index), update_cols] = value
    elif pert_control_name == "material" and index is not None:
        if layer_groups is None:
            raise ValueError("the 'material' control needs layer_groups (scenario 'material_layers')")
        df.loc[(slice(None), list(layer_groups[index])), update_cols] = value
    elif pert_control_name == "root_map" and index is not None:
        zones_in_vegi = zone[index + 1 == veg_map]
        df.loc[zones_in_vegi, update_cols] = value
    else:
        df.loc[:, update_cols] = value


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
