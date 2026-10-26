""" Managing Data Assimilation process pertubation. 
    Take a scenario dictionnary as input describing which and how parameters are perturbated
    Check consistency of the distribution
    Prepare for DA class
"""

from contextlib import contextmanager
import pandas as pd
import numpy as np
from scipy.stats import stats, qmc, norm
from scipy.stats import truncnorm
import matplotlib.pyplot as plt

from pyCATHY import cathy_utils as utils_CT
from pyCATHY.plotters import cathy_plots as plt_CT
import re 
import os
import logging

logger = logging.getLogger(__name__)

#%%
def check_distribution(parm2check):
    if parm2check["type_parm"] == "porosity":
        if parm2check["nominal"] < 0:
            raise ValueError

    if parm2check["type_parm"] == "thetar_VG":
        if parm2check["nominal"] > 1:
            logger.warning("nominal thetar_VG is above 1: value too high for a residual water content")

    # References for Archie paarmeters range
    # https://www.sciencedirect.com/science/article/pii/S1018363918306123
    # a : TYPE, optional
    #     Tortuosity factor. The default is [1.0].
    # m : TYPE, optional
    #     Cementation exponent. The default is [2.0]. (usually in the range 1.3 -- 2.5 for sandstones)
    # n : TYPE, optional
    #     Saturation exponent. The default is [2.0].


LAYER_LIKE = ("layers", "material", "zone_material")  # controls sampled as one value per list entry


def _require_nlayers(type_parm, nlayers):
    """The 'layers' control samples one value per layer: fail clearly without the count."""
    if nlayers is None:
        raise ValueError(
            f"{type_parm}: the 'layers' control needs the number of layers: pass nlayers=... "
            "to perturbate_parm (entries built by perturbate() already carry it)"
        )


def _apply_scenario_bounds(list_pert, scenario):
    """Store the scenario ``per_bounds`` on every entry, whatever branch built it.

    The ``zone`` and ``root_map`` branches created their entries with no bounds at all, so
    ``per_bounds`` was dropped silently (neither the prior nor the clipping of the analysed
    values saw it). Defaults an entry already carries (e.g. ZROOT >= 0) are kept for the
    side the scenario leaves open.
    """
    names = list(scenario.get("per_name", []))
    for entry in list_pert:
        root = re.match(r"^(.*?)(\d*)$", entry["type_parm"])
        name, suffix = root.group(1), root.group(2)
        if name not in names:
            continue
        i = names.index(name)
        nominal = scenario["per_nom"][i]
        listed = isinstance(nominal, list)
        entry["clip_min"], entry["clip_max"] = check4bounds(
            scenario, i, entry.get("clip_min"), entry.get("clip_max"),
            het_size=len(nominal) if listed else 1,
            het_nb=int(suffix) if (listed and suffix != "") else 0,
        )
        # Same for the parameter transformation (log, log-ratio, hyperbolic): the hand-built
        # entries (scalar Archie, ...) never copied the scenario ``transf_type``, so the
        # analysis updated them in linear space. An entry that already carries one is kept.
        transf = scenario.get("transf_type")
        if transf is not None and i < len(transf):
            spec = transf[i]
            if isinstance(spec, (list, tuple)):
                j = int(suffix) if suffix != "" else 0
                spec = spec[j] if j < len(spec) else None
            if isinstance(spec, str) and spec.strip().lower() == "none":
                spec = None
            if entry.get("transf_type") in (None, "None"):
                entry["transf_type"] = spec


def _n_classes(raster):
    """Number of zones (or vegetation classes) of a map; the nodata value -9999 is not one.

    Counting it made a map with nodata ask for one scenario entry too many, while the
    zone numbers written to the soil table (1..K) never include it.
    """
    values = np.unique(np.asarray(raster))
    return int(np.count_nonzero(values != -9999))


def _is_number(x):
    return isinstance(x, (int, float, np.integer, np.floating)) and not isinstance(x, bool)


def _is_unset(b):
    return b is None or (isinstance(b, str) and b.strip().lower() == "none")


def parse_bounds(bounds, where="per_bounds"):
    """Parse one bounds specification into ``(min, max)``; a side set to None is open.

    Accepted: ``None`` / ``"None"`` (no bounds), ``{min = a, max = b}`` (either key
    may be omitted) and ``[a, b]``. Anything else raises ``ValueError`` naming
    ``where``: a malformed bound used to be ignored silently (bare ``except``),
    leaving the parameter unbounded.
    """
    if _is_unset(bounds):
        return None, None
    if isinstance(bounds, dict):
        unknown = set(bounds) - {"min", "max"}
        if unknown:
            raise ValueError(
                f"{where}: unknown bound key(s) {sorted(unknown)}; use 'min' and 'max'."
            )
        lo, hi = bounds.get("min"), bounds.get("max")
    elif isinstance(bounds, (list, tuple)) and len(bounds) == 2:
        lo, hi = bounds
    else:
        raise ValueError(
            f"{where}: cannot read bounds {bounds!r}; expected 'None', "
            "{min = a, max = b} or [a, b]."
        )
    for side, v in (("min", lo), ("max", hi)):
        if v is not None and not _is_number(v):
            raise ValueError(f"{where}: {side} must be a number, got {v!r}.")
    if lo is not None and hi is not None and lo >= hi:
        raise ValueError(f"{where}: min ({lo}) must be lower than max ({hi}).")
    return lo, hi


def check4bounds(scenario, index, clip_min, clip_max, het_size=1, het_nb=None):
    """Bounds of parameter ``index`` of the scenario (``per_bounds`` is optional).

    ``clip_min`` / ``clip_max`` are the defaults, kept for a side the scenario leaves
    open. With ``het_size > 1`` (several layers, zones or classes) ``per_bounds[index]``
    is either one bounds spec shared by all of them, or a list with one spec per entry
    (``het_nb`` selects it).
    """
    per_bounds = scenario.get("per_bounds") if hasattr(scenario, "get") else None
    if per_bounds is None:
        return clip_min, clip_max
    name = scenario["per_name"][index] if "per_name" in scenario else index
    if index >= len(per_bounds):
        logger.warning("per_bounds has no entry for %s: parameter left unbounded.", name)
        return clip_min, clip_max

    spec = per_bounds[index]
    is_pair = isinstance(spec, (list, tuple)) and len(spec) == 2 and all(
        _is_unset(v) or _is_number(v) for v in spec
    )
    if (
        het_size > 1
        and isinstance(spec, (list, tuple))
        and not is_pair
        and not _is_unset(spec)
    ):
        if len(spec) != het_size:
            raise ValueError(
                f"per_bounds for {name}: {len(spec)} entries for {het_size} layers/zones."
            )
        spec = spec[het_nb]
    lo, hi = parse_bounds(spec, where=f"per_bounds for {name}")
    return (clip_min if lo is None else lo), (clip_max if hi is None else hi)


def create_dict_entry(
                         name,scenario,
                         index,nstri, 
                         scenario_nom, scenario_mean, scenario_sd, 
                         clip_min,clip_max,
                         NENS,
                         per_type,
                         pert_control_name=None,
                         transf_type=None
                      ):
    
    # per_type = scenario["per_type"][index]
    # if type(nstri) == int: 
    #     per_type = scenario["per_type"][index][nstri]
        
    # Sampling law. Historically always "normal" (the scenario's ``sampling_type`` was ignored here); it is now honoured
    # for Ks only (log-normal prior: ``per_sigma`` = std of ln(Ks)), so that existing scenario files keep their behaviour.
    sampling = "normal"
    if name == "Ks":
        _samp = scenario.get("sampling_type")
        if isinstance(_samp, (list, tuple)) and index < len(_samp) and _samp[index] in ("normal", "lognormal", "uniform"):
            sampling = _samp[index]

    dict_parm_entry = {
        "type_parm": name + str(nstri),
        "nominal": scenario_nom,  # nominal value
        "mean": scenario_mean,
        "sd": scenario_sd,
        "units": "",  # units
        "sampling_type": sampling,
        "ensemble_size": NENS,  # size of the ensemble
        "per_type": per_type,
        "transf_type": transf_type,
        "savefig": name + str(nstri) + ".png",
        "surf_zones_param": nstri,
        "clip_min": clip_min,
        "clip_max": clip_max,
        "pert_control_name": pert_control_name,
    }
    
    return dict_parm_entry

def process_perturbation_block(
    scenario, 
    NENS,
    index,
    nlayers=None, 
    zone_raster=None,
    veg_raster=None,
):
    """
    Process perturbation for a single parameter (given by index) from scenario dict.

    Returns a list of perturbation dict entries (does NOT append to external list).

    Parameters:
    - scenario: dict with keys 'per_name', 'per_nom', 'per_mean', 'per_sigma', 'per_type', optional 'per_bounds' and 'pert_control_name'.
    - NENS: int, number of ensemble members.
    - index: int, index of the parameter in scenario dict lists to process.
    - nlayers: int, required for 'layers' mode.
    - zone_raster: np.ndarray, required for 'zone' mode.
    - veg_raster: np.ndarray, required for 'root_map' mode.
    """
    entries = []
    pname = scenario["per_name"][index]
    pert_control_name = scenario.get('pert_control_name', [None]*len(scenario["per_name"]))[index]
    scenario_nom = scenario["per_nom"][index]
    scenario_mean = scenario["per_mean"][index]
    scenario_sd = scenario["per_sigma"][index]
    per_type = scenario["per_type"][index]
    transf_type = scenario.get("transf_type", [None] * len(scenario["per_name"]))[index]

    # Bounds (per_bounds, optional): one spec for the whole parameter, or one spec per layer / zone /
    # class / material (check4bounds selects it with het_size / het_nb). Same parser as everywhere else:
    # {min, max} and [min, max] are both accepted.

    # Handle perturbation by control type
    if isinstance(scenario_nom, list):
        if pert_control_name in LAYER_LIKE:
            if pert_control_name == 'layers' and nlayers is None:
                raise ValueError("nlayers must be provided for 'layers' perturbation control.")
            n_items = nlayers if pert_control_name == 'layers' else len(scenario_nom)
            for layer_idx in range(n_items):
                clip_min, clip_max = check4bounds(scenario, index, None, None,
                                                  het_size=n_items, het_nb=layer_idx)
                entry = create_dict_entry(
                    pname, scenario, index, layer_idx,
                    scenario_nom[layer_idx],
                    scenario_mean[layer_idx],
                    scenario_sd[layer_idx],
                    clip_min, clip_max, NENS,
                    per_type=per_type[layer_idx],
                    pert_control_name=pert_control_name,
                    transf_type = transf_type
                )
                check_distribution(entry)
                entries.append(entry)

        elif pert_control_name == 'zone':
            if zone_raster is None:
                raise ValueError("zone_raster must be provided for 'zone' perturbation control.")
            for zone_idx in range(_n_classes(zone_raster)):
                clip_min, clip_max = check4bounds(scenario, index, None, None,
                                                  het_size=_n_classes(zone_raster), het_nb=zone_idx)
                z_nom = scenario_nom[zone_idx]
                z_mean = scenario_mean[zone_idx]
                z_sd = scenario_sd[zone_idx]
                z_type = per_type[zone_idx]

                entry = create_dict_entry(
                    pname, scenario, index, zone_idx,
                    z_nom, z_mean, z_sd,
                    clip_min, clip_max, NENS,
                    per_type=z_type,
                    pert_control_name=pert_control_name,
                    transf_type = transf_type
                )
                check_distribution(entry)
                entries.append(entry)

        elif pert_control_name == 'root_map':
            if veg_raster is None:
                raise ValueError("veg_raster must be provided for 'root_map' perturbation control.")
            for veg_idx in range(_n_classes(veg_raster)):
                clip_min, clip_max = check4bounds(scenario, index, None, None,
                                                  het_size=_n_classes(veg_raster), het_nb=veg_idx)
                r_nom = scenario_nom[veg_idx]
                r_mean = scenario_mean[veg_idx]
                r_sd = scenario_sd[veg_idx]
                r_type = per_type[veg_idx] if isinstance(per_type, list) else per_type

                entry = create_dict_entry(
                    pname, scenario, index, veg_idx,
                    r_nom, r_mean, r_sd,
                    clip_min, clip_max, NENS,
                    per_type=r_type,
                    pert_control_name=pert_control_name,
                    transf_type = transf_type
                )
                check_distribution(entry)
                entries.append(entry)

        else:
            raise ValueError(
                f"Unknown list-type perturbation control '{pert_control_name}' for parameter '{pname}'. "
                "Expected 'layers', 'zone', or 'root_map'."
            )

    else:
        # Scalar perturbation
        clip_min, clip_max = check4bounds(scenario, index, None, None)
        entry = create_dict_entry(
            pname, scenario, index, '',
            scenario_nom, scenario_mean, scenario_sd,
            clip_min, clip_max, NENS,
            per_type=per_type,
            pert_control_name=pert_control_name,
            transf_type = transf_type
        )
        check_distribution(entry)
        entries.append(entry)

    return entries



def perturbate(simu_DA, scenario, NENS,
               **kwargs):
    """Build the list of perturbation entries described by a scenario.

    Each entry is a dictionary (``type_parm``, ``nominal``, ``mean``, ``sd``,
    ``clip_min``, ``clip_max``, ``sampling_type``, ``per_type``, ``transf_type``,
    ``pert_control_name``, ...) later sampled by ``perturbate_parm``.

    Parameters
    ----------
    simu_DA : DA
        Simulation, read for the soil table, the zone map and the vegetation map.
    scenario : dict
        One scenario of ``scenario.toml``: lists aligned with ``per_name``
        (``per_nom``, ``per_mean``, ``per_sigma``, ``per_type``, ``sampling_type``,
        ``transf_type``, ``per_bounds``, ``pert_control_name``) and ``listUpdateParm``.
    NENS : int
        Ensemble size.
    nzones, nlayers : int, optional (keyword)
        Override the number of zones / layers read from the soil table.
    validate_scenario : bool, optional (keyword, default True)
        Check the scenario before building and the names of ``listUpdateParm`` after
        (``pyCATHY.DA.scenario``). False skips both checks.

    Notes
    -----
    ``pert_control_name`` selects what an entry applies to; the entries of a list
    are written in order, so overlapping entries override each other (see
    ``scenario.describe_ensemble``):

    - ``None``: the whole domain (key ``<name>``, or ``<name>0`` for Feddes parameters
      such as ``ZROOT``);
    - ``"layers"``: one entry per layer, keys ``<name>0``, ``<name>1`` ... (layer index);
    - ``"zone"``: one entry per zone, key index ``k`` is zone number ``k + 1``;
    - ``"root_map"``: one entry per vegetation class, key index ``k`` is class ``k + 1``.

    Zone and class numbers must run from 1 to K without gaps.
    ``per_bounds`` is ``None``, ``{min = a, max = b}`` or ``[a, b]`` (see
    ``parse_bounds``); bounds are used when sampling the prior and again to clip the
    analysed values.
    """
    
    df_SPP, df_FP = simu_DA.read_inputs('soil')
    zone_raster, zone_header = simu_DA.read_inputs('zone')
    veg_raster, veg_header = simu_DA.read_inputs('root_map')
        
    if 'nzones' in kwargs: 
        nzones = kwargs.pop('nzones')
    else:
        nzones = len(df_SPP.index.get_level_values(0).unique())
    if 'nlayers' in kwargs: 
        nlayers = kwargs.pop('nlayers')
    else:
        nlayers = len(df_SPP.index.get_level_values(1).unique())
        # nlayers_PERMX = len(df_SPP['PERMX'].unique())

    # fail early, with every problem listed, instead of silently ignoring a misspelled
    # parameter or a bad bound (validate_scenario=False skips the check)
    check_scenario = kwargs.pop("validate_scenario", True)
    if check_scenario:
        from pyCATHY.DA.scenario import validate_scenario as _validate

        try:
            n_feddes_rows = len(simu_DA.read_inputs("soil")[1].index.unique())
        except Exception:  # no Feddes table readable: ZROOT / PCREF layout is then not checked
            n_feddes_rows = None
        _, scenario_warnings = _validate(
            scenario, nlayers=nlayers, zone_raster=zone_raster, veg_raster=veg_raster,
            n_feddes_rows=n_feddes_rows,
        )
        for message in scenario_warnings:
            logger.warning("scenario: %s", message)

    #%%
    list_pert = []

    #%% Initial and boundary conditions parameters
    # ------------------------------------------------------------------------
    
    if "WTPOSITION" in scenario["per_name"]:

        index = scenario["per_name"].index("WTPOSITION")

        clip_min = 0
        clip_max = simu_DA.dem_parameters['base']
        
        clip_min, clip_max = check4bounds(
            scenario,
            index,
            clip_min,
            clip_max,
        )
            
        WTPOSITION = {
            "type_parm": "WTPOSITION",
            "nominal": scenario["per_nom"][index],  # nominal value
            "mean": scenario["per_mean"][index],
            "sd": scenario["per_sigma"][index],
            "units": "Water table position $(m)$",  # units
            "sampling_type": "normal",
            "ensemble_size": NENS,  # size of the ensemble
            "per_type": scenario["per_type"][index],
            "savefig": "WTPOSITION.png",
            "clip_min": clip_min,
            "clip_max": clip_max,
        }
        list_pert.append(WTPOSITION)
        
        
    if "ic" in scenario["per_name"]:

        index = scenario["per_name"].index("ic")

        pert_control_name = None
        if 'pert_control_name' in scenario:
            pert_control_name = scenario['pert_control_name'][index]
        
        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        
        
        clip_min = None
        clip_max = None
        if type(scenario["per_nom"][index]) is list:    
            if pert_control_name=='layers':
                for nstri in range(nlayers):
                    
                    if nlayers > 1:
                        scenario_nom = scenario["per_nom"][index][nstri]
                        scenario_mean = scenario["per_mean"][index][nstri]
                        scenario_sd = scenario["per_sigma"][index][nstri]
                        scenario_per_type = scenario["per_type"][index][nstri]
                        transf_type = scenario["transf_type"][index][nstri]
        
                    ic = create_dict_entry(
                                'ic',scenario,index,nstri, 
                                scenario_nom, scenario_mean, scenario_sd, 
                                clip_min,clip_max,
                                NENS,
                                per_type=scenario_per_type,
                                transf_type = transf_type,
                                pert_control_name=pert_control_name
                                )
                    
                    list_pert.append(ic)
                    
            elif pert_control_name=='zone':
                for zonei in range(_n_classes(zone_raster)):
                    # Adjust bounds and retrieve parameters for the current zone
                    clip_min, clip_max = check4bounds(
                    scenario, index, None, None, het_size=_n_classes(zone_raster), het_nb=zonei
                )
                    zone_nom, zone_mean, zone_sd = scenario_nom[zonei], scenario_mean[zonei], scenario_sd[zonei]
                    scenario_per_type = scenario["per_type"][index][zonei]

                    ic = create_dict_entry(
                        'ic', scenario, index, zonei, zone_nom, zone_mean, zone_sd, 
                        clip_min, clip_max, NENS,
                        per_type=scenario_per_type,
                        pert_control_name=pert_control_name,
                        transf_type = scenario["transf_type"][index]
                    )
                    check_distribution(ic)
                    list_pert.append(ic)        
            else:
                raise ValueError(f"List-type perturbation provided, but unknown control name: '{pert_control_name}'. "
                         f"Expected 'layers' or 'zone'.")
        else: 
            # bounds are optional; e.g. a pressure head ic must stay negative
            # for an unsaturated initial state (IPOND=0)
            clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
            ic = create_dict_entry(
                        'ic',scenario,index,'', 
                        scenario_nom, scenario_mean, scenario_sd, 
                        clip_min,clip_max,
                        NENS,
                        per_type=scenario["per_type"][index],
                        transf_type = scenario["transf_type"][index]
                        )
            list_pert.append(ic)

    if "atmbc" in scenario["per_name"]:
        index = scenario["per_name"].index("atmbc")

        atmbc_nom = scenario["per_nom"]  # 1e-4
        atmbc_mean = scenario["per_mean"]  # sampling mean
        atmbc_sd = scenario["per_sigma"]  # 1.53
        

        df = simu_DA.atmbc['atmbc_df']
        # Get unique 'time' values
        unique_times = df['time'].unique()

        # type_parm
        # extract value of each node of the mesh for all the atmbc times
        for value_idx in range(len(df) // len(unique_times)):
            # print(value_idx)
            values = []
            for time in unique_times:
                # print(time)
                value = df.loc[df['time'] == time, 'value'].iloc[value_idx]
                values.append(value)
            # print(values)
            atmbc = {
                "type_parm": "atmbc" + str(value_idx),
                "nominal": atmbc_nom[index], 
                "units": "$mm.h^{-1}$",
                "mean": atmbc_nom[index],
                "sd": atmbc_sd[index],
                "atmbc_units": "atmospheric forcing $mm.h^{-1}$",
                "sampling_type": 'normal',
                "ensemble_size": NENS,
                "per_type": scenario["per_type"][index],
                "time_variable": True,
                "data2perturbate": {'time': simu_DA.atmbc['time'], 'VALUE': values},
                "time_decorrelation_len": scenario["time_decorrelation_len"][index],
                "savefig": "atmbc.png"
            }
            
            list_pert.append(atmbc)

    #%% Petro-Pedophysical parameters
    # ------------------------------------------------------------------------
    
    if "SATCELL_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("SATCELL_VG")
        SATCELL_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        list_pert.extend(SATCELL_VG)


    if "thetar_VG" in scenario["per_name"]:
        
        index = scenario["per_name"].index("thetar_VG")
        thetar_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        list_pert.extend(thetar_VG)
        
    if "alpha_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("alpha_VG")          
        alpha_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        list_pert.extend(alpha_VG)

    if "n_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("n_VG")
        n_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        
        list_pert.extend(n_VG)


    # if "VGP" in scenario["per_name"]:
    #     # - 'PERMX' (NSTR, nstriONE): saturated hydraulic conductivity - xx
    #     # - 'PERMY' (NSTR, nstriONE): saturated hydraulic conductivity - yy
    #     # - 'PERMZ' (NSTR, nstriONE): saturated hydraulic conductivity - zz
    #     # - 'ELSTOR' (NSTR, nstriONE): specific storage
    #     # - 'POROS'  (NSTR, nstriONE): porosity (moisture content at saturation) = \thetaS

    #     # retention curves parameters VGN, VGRMC, and VGPSAT
    #     # - 'VGNCELL' (NSTR, nstriONE): van Genuchten curve exponent  = n
    #     # - 'VGRMCCELL' (NSTR, nstriONE): residual moisture content = \thetaR
    #     # - 'VGPSATCELL' (NSTR, nstriONE): van Genuchten curve exponent -->
    #     #                               VGPSAT == -1/alpha (with alpha expressed in [L-1]);
    #     # ['ks','ss','phi','thetar','alpha','n'])
    #     # ['PERMX','ELSTOR','POROS','VGRMCCELL','VGPSATCELL','VGNCELL'])

    #     # need to account for transformation of the parameters see Boto

    #     VGP_units = ["", "", "", "", "", ""]

    #     for i, p in enumerate(["ks", "ss", "phi", "thetar", "alpha", "r"]):
    #         index = scenario["per_name"].index("VGP")

    #         if p not in scenario["per_nom"][index]:
    #             pass
    #         else:
    #             p_VGP = {
    #                 "type_parm": p,
    #                 "nominal": scenario["per_nom"][index][p],  # nominal value
    #                 "mean": scenario["per_mean"][index][p],
    #                 "sd": scenario["per_sigma"][index][p],
    #                 "units": VGP_units[i],  # units
    #                 "sampling_type": "normal",
    #                 "ensemble_size": NENS,  # size of the ensemble
    #                 "per_type": scenario["per_type"][index],
    #                 "savefig": "_VGP" + p + ".png",
    #             }
    #             # list_pert.append(p_VGP)

    Archie_2pert = [ele for ele in scenario["per_name"] if ("Archie" in ele)]
    if len(Archie_2pert) > 0:
        # Initialement, chacun des membres de l’ensemble aurait des paramètres d’Archie
        # légèrement différents, qui seraient utilisés dans la fonction d’observation afin de simuler les
        # données à partir de l’état a priori. Ces paramètres seraient ensuite modifiés à l’étape de
        # mise à jour, de la même façon que les variables de teneur en eau.
        # rFluid=[1.0],a=[1.0],m=[2.0],n=[2.0]
        # simu_DA.Archie_parms
        # literature example:
        # ---------------------------------------------------------------------
        # https://www.sciencedirect.com/science/article/pii/S0169772220302680#tf0005
        # c = (mean = 1.6,sd = 0.5,min = 0.0, max = 2.0),
        # m = (mean = 2.5,sd = 0.8,min = 0.0, max = 3.0)

        # archie_units = ["", "", "", ""]
        archie_units = {
            "porosity": "-",                   # dimensionless (volume fraction)
            "a_Archie": "-",                          # tortuosity factor (dimensionless)
            "m_Archie": "-",                          # cementation exponent (dimensionless)
            "n_Archie": "-",                          # saturation exponent (dimensionless)
            "rFluid_Archie": "ohm.m",           # resistivity of water
        }


        for i, p in enumerate(Archie_2pert):
            
            index = scenario["per_name"].index(p)      
            pert_control_name = None
            if 'pert_control_name' in scenario:
                pert_control_name = scenario['pert_control_name'][index]
                
            
            scenario_nom = scenario["per_nom"][index]
            scenario_mean = scenario["per_mean"][index]
            scenario_sd = scenario["per_sigma"][index]
            scenario_per_type = scenario["per_type"][index]



            if pert_control_name is not None and pert_control_name == 'layers':
                for nstri in range(nlayers):
                    # Adjust bounds and retrieve parameters for the current layer
                    # clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                    layer_nom, layer_mean, layer_sd = scenario_nom[nstri], scenario_mean[nstri], scenario_sd[nstri]
                    # Create porosity entry for this layer and validate distribution
                    p_archie = create_dict_entry(
                        p, scenario, index, nstri, layer_nom, layer_mean, layer_sd, 
                        clip_min, clip_max, NENS,
                        per_type=scenario_per_type[nstri],
                    )
                    check_distribution(p_archie)
                    list_pert.append(p_archie)
            else:
                # pass
                index = scenario["per_name"].index(p)
        
                p_archie = {
                    "type_parm": p,
                    "nominal": scenario["per_nom"][index],  # nominal value
                    "mean": scenario["per_mean"][index],
                    "sd": scenario["per_sigma"][index],
                    "units": archie_units[p],  # units
                    "sampling_type": "normal",
                    "ensemble_size": NENS,  # size of the ensemble
                    "per_type": scenario["per_type"][index],
                    "savefig": "_Archie" + p + ".png",
                    'pert_control_name':pert_control_name

                }
                list_pert.append(p_archie)

    #%% SOIL parameters
    # ------------------------------------------------------------------------


    if "Ks" in scenario["per_name"]:
        # same builder as the van Genuchten parameters: every control (zone, root_map, layers,
        # material, zone_material), the transform of the scenario (log) and the scenario bounds.
        # PERMX, PERMY and PERMZ are all set to the sampled value (param_tables._SOIL_PARM_COLUMNS).
        index = scenario["per_name"].index("Ks")
        list_pert.extend(process_perturbation_block(scenario,
                                                    NENS=NENS,
                                                    index=index,
                                                    nlayers=nlayers,
                                                    zone_raster=zone_raster,
                                                    veg_raster=veg_raster))

    # Check if 'porosity' is a parameter in the scenario and get its index if present
    if "porosity" in scenario["per_name"]:
        index = scenario["per_name"].index("porosity")
    
        # Set perturbation order, defaulting to ['zone', 'root_map', 'layers'] if not specified
        
        
        pert_control_name = None
        if 'pert_control_name' in scenario:
            pert_control_name = scenario['pert_control_name'][index]
        
        
        # Retrieve nominal, mean, and sigma values for the scenario
        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        scenario_per_type = scenario["per_type"][index]
    
        # Define clipping bounds for soil porosity
        clip_min, clip_max = 0.2, 0.7
    
        # Perturb by 'zone' if specified
        if pert_control_name=='zone':
            for zonei in range(_n_classes(zone_raster)):
                # Adjust bounds and retrieve parameters for the current zone
                clip_min, clip_max = check4bounds(
                    scenario, index, None, None, het_size=_n_classes(zone_raster), het_nb=zonei
                )
                zone_nom, zone_mean, zone_sd = scenario_nom[zonei], scenario_mean[zonei], scenario_sd[zonei]
    
                # Create porosity entry for this zone and validate distribution
                porosity = create_dict_entry(
                    'porosity', scenario, index, zonei, zone_nom, zone_mean, zone_sd, 
                    clip_min, clip_max, NENS,
                    per_type=scenario_per_type[zonei],
                    pert_control_name=pert_control_name
                )
                check_distribution(porosity)
                list_pert.append(porosity)
    
        # Perturb by 'layers' if specified
        elif pert_control_name=='layers':
            for nstri in range(nlayers):
                # Adjust bounds and retrieve parameters for the current layer
                clip_min, clip_max = check4bounds(
                    scenario, index, None, None, het_size=nlayers, het_nb=nstri
                )
                layer_nom, layer_mean, layer_sd = scenario_nom[nstri], scenario_mean[nstri], scenario_sd[nstri]
    
                # Create porosity entry for this layer and validate distribution
                porosity = create_dict_entry(
                    'porosity', scenario, index, nstri, layer_nom, layer_mean, layer_sd, 
                    clip_min, clip_max, NENS,
                    per_type=scenario_per_type[nstri],
                    pert_control_name=pert_control_name
                )
                check_distribution(porosity)
                list_pert.append(porosity)
    
        # Default case: apply global perturbation if neither 'zone' nor 'layers' specified
        else:
            clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
            porosity = create_dict_entry(
                'porosity', scenario, index, '', scenario_nom, scenario_mean, scenario_sd, 
                clip_min, clip_max, NENS,
                per_type=scenario["per_type"][index],
            )
            check_distribution(porosity)
            list_pert.append(porosity)
    


    #%% Plant parameters
    # ------------------------------------------------------------------------

    if "PCWLT" in scenario["per_name"]:
        index = scenario["per_name"].index("PCWLT")


        clip_min = -180
        clip_max = None
            
        for nstri in range(len(simu_DA.soil["PCWLT"])):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil["PCWLT"]),
                het_nb=nstri,
            )
            
            PCWLT = {
                "type_parm": "PCWLT" + str(nstri),
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": "$m$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "PCWLT.png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(PCWLT)
            
    if "OMGC" in scenario["per_name"]:
        index = scenario["per_name"].index("OMGC")


        clip_min = 0
        clip_max = 1
            
        for nstri in range(len(simu_DA.soil["OMGC"])):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil["OMGC"]),
                het_nb=nstri,
            )
            
            OMGC = {
                "type_parm": "OMGC" + str(nstri),
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": "$?$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "OMGC.png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(OMGC)
            
            
    if "PZ" in scenario["per_name"]:
        index = scenario["per_name"].index("PZ")


        clip_min = 0
        clip_max = None
            
        for nstri in range(len(simu_DA.soil["PZ"])):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil["PZ"]),
                het_nb=nstri,
            )
            
            PZ = {
                "type_parm": "PZ" + str(nstri),
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": "$?$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "PZ.png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(PZ)
            
    if "PCREF" in scenario["per_name"]:
        index = scenario["per_name"].index("PCREF")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        
        clip_min = -30
        clip_max = -1

        df_SPP, df_FP = simu_DA.read_inputs('soil')

        for nstri in range(len(df_FP.index.unique())):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(df_FP.index.unique()),
                het_nb=nstri,
            )
        
        
            if len(df_FP.index.unique()) > 1:
                scenario_nom = scenario["per_nom"][index][nstri]
                scenario_mean = scenario["per_mean"][index][nstri]
                scenario_sd = scenario["per_sigma"][index][nstri]

            PCREF = {
                "type_parm": "PCREF" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$m$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "PCREF" + str(nstri) + ".png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(PCREF)

    if "ZROOT" in scenario["per_name"]:
        index = scenario["per_name"].index("ZROOT")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        if "sampling_type" in scenario:
            scenario_sampling = scenario["sampling_type"][index]
        else:
            scenario_sampling = "lognormal"

        clip_min = 0
        clip_max = simu_DA.dem_parameters['base']
        
        df_SPP, df_FP = simu_DA.read_inputs('soil')
        
        for nzonesi in range(len(df_FP.index.unique())):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(df_FP.index.unique()),
                het_nb=nzonesi,
            )

            if len(df_FP.index.unique()) > 1:
                scenario_nom = scenario["per_nom"][index][nzonesi]
                scenario_mean = scenario["per_mean"][index][nzonesi]
                scenario_sd = scenario["per_sigma"][index][nzonesi]

            ZROOT = {
                "type_parm": "ZROOT" + str(nzonesi),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$m$",  # units
                "sampling_type": scenario_sampling,
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "ZROOT" + str(nzonesi) + ".png",
                "surf_zones_param": nzonesi,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(ZROOT)
            # nb_surf_nodes = 110

    for entry in list_pert:
        if entry.get("pert_control_name") == "layers":
            entry["nlayers"] = nlayers  # perturbate_parm needs it to sample the layers
        elif entry.get("pert_control_name") == "material":
            groups = scenario.get("material_layers")
            if groups is None:
                raise ValueError("pert_control_name 'material' needs the scenario key 'material_layers' "
                                 "(list of lists of layer numbers, one list per material)")
            entry["layer_groups"] = [list(g) for g in groups]
            entry["nlayers"] = len(groups)  # one sampled value per material
        elif entry.get("pert_control_name") == "zone_material":
            zgroups = scenario.get("zone_material_layers")
            if zgroups is None:
                raise ValueError("pert_control_name 'zone_material' needs the scenario key "
                                 "'zone_material_layers' (list of [zone, [layers]], zone 1-based, layers 0-based)")
            entry["layer_groups"] = [[int(z), [int(l) for l in ls]] for z, ls in zgroups]
            entry["nlayers"] = len(zgroups)  # one sampled value per (zone, material)
    _apply_scenario_bounds(list_pert, scenario)

    _n_lm = sum(e.get("pert_control_name") in LAYER_LIKE for e in list_pert)
    if _n_lm:
        print(f"[DA-PARAM] perturbate: {_n_lm} layer/material entries sampled with the SEEDABLE random "
              "sampler (np.random.seed reproduces the prior; the Sobol sampler is no longer the default)",
              flush=True)

    if check_scenario:
        # authoritative: the names of listUpdateParm must be keys that were really built
        # (a scalar Ks is "Ks", a scalar ZROOT is "ZROOT0": the naming depends on the parameter)
        from pyCATHY.DA.scenario import check_update_parameters

        check_update_parameters(list_pert, scenario.get("listUpdateParm", []))

    return list_pert


#%%


# SAMPLING distribution
# ----------------------

def _std_bounds(clip_min, clip_max, loc, scale):
    """Bounds in units of the standard deviation; an open side (None) is infinite."""
    a = -np.inf if clip_min is None else (clip_min - loc) / scale
    b = np.inf if clip_max is None else (clip_max - loc) / scale
    return a, b


def _has_bounds(entry):
    """True when an entry carries at least one bound (a max-only bound counts)."""
    return entry.get("clip_min") is not None or entry.get("clip_max") is not None


def sampling_dist_trunc(myclip_a, myclip_b, ensemble_size, **kwargs):
    """Truncated normal sample; ``myclip_a`` / ``myclip_b`` may be None (open side)."""
    # https://stackoverflow.com/questions/18441779/how-to-specify-upper-and-lower-limits-when-using-numpy-random-normal
    loc, scale = kwargs["loc"], kwargs["scale"]
    if scale == 0:  # no spread: every member is the mean (kept inside the bounds)
        lo = -np.inf if myclip_a is None else myclip_a
        hi = np.inf if myclip_b is None else myclip_b
        return np.full(ensemble_size, float(np.clip(loc, lo, hi)))
    a, b = _std_bounds(myclip_a, myclip_b, loc, scale)
    return truncnorm(a, b, loc=loc, scale=scale).rvs(ensemble_size)


# Lognormal prior: what does the nominal value ``mean`` mean?
#   "median" (historical default): ln(x) ~ N(ln(mean), sd^2)  -> median = mean, ensemble mean = mean * exp(sd^2/2)
#                                  (+5 % at sd 0.3, +28 % at 0.7, +65 % at 1.0)
#   "mean"                       : ln(x) ~ N(ln(mean) - sd^2/2, sd^2) -> ensemble mean = mean, median = mean * exp(-sd^2/2)
# With clip bounds the law is truncated and the realised mean is no longer exact in either mode
# (see ``param_tables.parameter_summary``: ``prior_mean_over_nominal``).
LOGNORMAL_CENTERS = ("median", "mean")
_lognormal_center_default = "median"


def _check_lognormal_center(center):
    if center not in LOGNORMAL_CENTERS:
        raise ValueError(f"unknown lognormal_center {center!r}: use one of {LOGNORMAL_CENTERS}")
    return center


def set_lognormal_center(center):
    """Set the process-wide default (``"median"`` or ``"mean"``); returns the previous one."""
    global _lognormal_center_default
    previous = _lognormal_center_default
    _lognormal_center_default = _check_lognormal_center(center)
    return previous


@contextmanager
def using_lognormal_center(center):
    """``with using_lognormal_center("mean"): ...`` temporarily centres the lognormal priors on the MEAN."""
    previous = set_lognormal_center(center)
    try:
        yield
    finally:
        set_lognormal_center(previous)


def lognormal_log_location(mean, sd, center=None):
    """Location ``mu`` of ``ln(x) ~ N(mu, sd^2)`` so that ``mean`` is the median (``"median"``) or the mean (``"mean"``).

    ``center=None`` uses the process default (``set_lognormal_center``, ``"median"`` unless changed). The default
    returns exactly ``np.log(mean)``, so existing results are reproduced bit for bit.
    """
    center = _check_lognormal_center(_lognormal_center_default if center is None else center)
    if center == "median":
        return np.log(mean)
    return np.log(mean) - 0.5 * np.asarray(sd, dtype=float) ** 2


def sampling_dist(sampling_type, mean, sd, ensemble_size, **kwargs):
    # No reseeding here: a fixed seed per call gave every parameter the same standard-normal
    # draws (a rank-1 parameter ensemble). Call np.random.seed(n) once beforehand to reproduce.
    # kwargs: ``minmax_uni`` (uniform), ``center`` ("median" | "mean", lognormal; None = process default)
    if sampling_type == "lognormal":
        parm_sampling = np.random.lognormal(lognormal_log_location(mean, sd, kwargs.get("center")),
                                            sigma=sd, size=ensemble_size)
    elif sampling_type == "normal":
        # parm_sampling = np.random.normal(mean, sd, size=ensemble_size)
        parm_sampling = np.random.normal(mean, scale=sd, size=ensemble_size)
    elif sampling_type == "uniform":
        # without explicit limits: the uniform law with this mean and sd (half-width sqrt(3) sd)
        minmax_uni = kwargs.get("minmax_uni")
        if minmax_uni is None:
            half = np.sqrt(3.0) * sd
            minmax_uni = (mean - half, mean + half)
        parm_sampling = np.random.uniform(minmax_uni[0], minmax_uni[1], ensemble_size)
    return parm_sampling


# PERTUBATE distribution
# ----------------------


def perturbate_dist(parm, per_type, parm_sampling, ensemble_size):
    # pertubate
    parm_mat = np.ones(ensemble_size) * parm["nominal"]
    if per_type is None or per_type == "None":
        parm_per_array = parm_sampling
    elif per_type == "multiplicative":
        parm_per_array = parm_mat * parm_sampling
    elif per_type == "additive":
        parm_per_array = parm_mat + parm_sampling
    else:  # used to fall through and raise an UnboundLocalError
        raise ValueError(
            f"{parm.get('type_parm', 'parameter')}: unknown per_type {per_type!r}: use None, "
            "'additive' or 'multiplicative' (a single value, not a list)"
        )
    return parm_per_array


def Carsel_Parrish_VGN_pert():
    cholesky_diag_mat = np.diag(3)
    pass


def Archie_pert_rules(
    parm, type_parm, ensemble_size, mean, sd, per_type, sampling_type
):
    # a : TYPE, optional
    #     Tortuosity factor. The default is [1.0].
    # m : TYPE, optional
    #     Cementation exponent. The default is [2.0]. (usually in the range 1.3 -- 2.5 for sandstones)
    # n : TYPE, optional
    #     Saturation exponent. The default is [2.0].

    # default truncations, replaced side by side by the scenario bounds (per_bounds) of the entry
    defaults = {"rFluid": (0, np.inf), "a": (0, 2.5), "m": (1.3, 3), "n": (1.3, 3)}
    family = "rFluid" if type_parm.startswith("rFluid") else type_parm.split("_")[0]
    if family in defaults:
        lo, hi = defaults[family]
        if parm.get("clip_min") is not None:
            lo = parm["clip_min"]
        if parm.get("clip_max") is not None:
            hi = parm["clip_max"]
        parm_sampling = sampling_dist_trunc(
            myclip_a=lo, myclip_b=hi, ensemble_size=ensemble_size, loc=mean, scale=sd
        )
    else:
        parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size, center=parm.get("lognormal_center"))

    parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    return parm_sampling, parm_per_array


def VG_pert_rules(
    var_per_2add,
    parm,
    type_parm,
    ensemble_size,
    mean,
    sd,
    per_type,
    sampling_type,
    **kwargs,
):

    # if "Carsel_Parrish_VGN_pert" in kwargs:
    #     utils.Carsel_Parrish_1988(soilTexture=None)
    #     Carsel_Parrish_VGN_pert()
    # else:
    #     if "clip_min" in var_per_2add[type_parm].keys():
    #         parm_sampling = sampling_dist_trunc(
    #             myclip_a=var_per_2add[type_parm]["clip_min"],
    #             myclip_b=var_per_2add[type_parm]["clip_max"],
    #             ensemble_size=ensemble_size,
    #             loc=mean,
    #             scale=sd,
    #         )
    #     else:
    #         parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size)
    #     parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    
    """
    Sample a single van Genuchten physical parameter and transform it to Gaussian space via Johnson transformation.

    Parameters:
    - parm: dict
        Must include:
            - 'nominal': float
            - 'sampling_type': 'lognormal' | 'normal' | 'uniform'
            - 'mean': mean of sampling distribution (in physical space)
            - 'sd': standard deviation of sampling
            - 'per_type': 'additive' | 'multiplicative' | None
            - 'type': 'LN' | 'SB' | 'SU' (Johnson transform type)
            - 'A': Johnson transform A parameter (for SB)
            - 'B': Johnson transform B parameter (for SB)
            - Optional: 'minmax_uni': [min, max] (if uniform)

    - ensemble_size: int
    - seed: int

    Returns:
    - v_array: (ensemble_size,) ndarray of physical parameter samples
    - y_array: (ensemble_size,) ndarray of transformed (Gaussian) samples
    """
    
    logger.info(
        "van Genuchten parameters drawn independently here; for their mutual correlation (Botto et al. 2018) call "
        "pyCATHY.DA.vgp_correlation.apply_vgp_correlation on the result"
    )
    
    # (no reseeding: see sampling_dist)

    # --- Sampling ---
    sampling_type = parm["sampling_type"]
    mean = parm["mean"]
    sd = parm["sd"]
    per_type = parm.get("per_type", None)
    transf_type = parm.get("transf_type", None)

    kwargs = {}
    if sampling_type == "uniform" and "minmax_uni" in parm:
        kwargs["minmax_uni"] = parm["minmax_uni"]

    # parm_sampling = sampling_dist(
    #     sampling_type=sampling_type,
    #     mean=mean,
    #     sd=sd,
    #     ensemble_size=ensemble_size,
    #     **kwargs
    # )

    if sampling_type == "lognormal":
        # log-normal prior: ln(x) ~ N(ln(mean), sd^2) truncated to [ln(clip_min), ln(clip_max)]; sd = std of ln(x)
        lo, hi = parm.get("clip_min"), parm.get("clip_max")
        parm_sampling = np.exp(sampling_dist_trunc(
            myclip_a=None if (lo is None or lo <= 0) else np.log(lo),
            myclip_b=None if (hi is None or hi <= 0) else np.log(hi),
            ensemble_size=ensemble_size,
            loc=lognormal_log_location(parm["mean"], sd, parm.get("lognormal_center")),
            scale=sd,
        ))
    else:
        # one call: an open side (None) is infinite, so a missing or max-only bound works too
        parm_sampling = sampling_dist_trunc(
            myclip_a=parm.get("clip_min"),
            myclip_b=parm.get("clip_max"),
            ensemble_size=ensemble_size,
            loc=parm['mean'],
            scale=sd,
        )

    # --- Perturbation (physical space) ---
    v_array = perturbate_dist(
        parm=parm,
        per_type=per_type,
        parm_sampling=parm_sampling,
        ensemble_size=ensemble_size
    )

    # # --- Johnson transform (to Gaussian space) ---
    # y_array = Johnson1970_transform(
    #     v_array,
    #     transform_type=transf_type,
    #     A=parm['clip_min'],
    #     B=parm['clip_max']
    # )

    return parm_sampling, v_array

def atmbc_pert_rules(
        var_per_2add,
    parm,
    type_parm,
    ensemble_size,
    mean,
    sd,
    per_type,
    sampling_type,
    **kwargs,
):

    atmbc_times = parm["data2perturbate"]["time"]
    atmbc_values = parm["data2perturbate"]["VALUE"]
    parm_per_ti = []
    white_noise =  np.random.normal(0, 1, size=(ensemble_size, 
                                                    len(atmbc_times)
                                                    )
                                    )

    Tau = parm["time_decorrelation_len"]

    # Compute perturbed hyetographs
    q = np.zeros_like(white_noise)
    q[:, 0] = white_noise[:, 0]  # Initialize with white noise
    
    deltaT = np.diff(np.r_[atmbc_times,atmbc_times[-1]+86400])
    gamma = np.array([1 - (dt / Tau) for dt in deltaT])
    
    for i in range(1, len(deltaT)):
        q[:, i] = (gamma[i] * q[:, i-1]) + (np.sqrt(1 - gamma[i] ** 2) * white_noise[:, i])

    q_Evenson = q
  
    # Transform using log-normal function
    # m_i, s_i = 1, 0.1  # Example mean and std for transformation
    q_log = np.log(mean) - 0.5 * np.log(((sd / mean) ** 2) + 1) + q_Evenson * np.sqrt(np.log(((sd / mean) ** 2) + 1))
    
    # Back transformation
    q_logback = np.exp(q_log)
    
    # Compute final perturbed forcing
    parm_per_ti = np.multiply(atmbc_values, q_logback)
    
    # Plot perturbed hyetographs
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(ensemble_size):
        ax.step(atmbc_times, parm_per_ti[i, :], where="post", alpha=0.5, label=f"Perturbation {i+1}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Forcing (Perturbed)")
    plt.title("Perturbed Hyetographs")
    # plt.legend()
    # plt.show()

    # var_per_2add[type_parm] = parm_per_ti
    key = "time_variable_perturbation"
    var_per_2add[type_parm][key] = parm_per_ti
    fig.savefig(os.path.join(os.getcwd(), kwargs["savefig"]), dpi=350)


    return white_noise, parm_per_ti


# def Johnson1970():
#     print("not yet implemented - see Botto 2018")


def Johnson1970_transform(V, transform_type, A=None, B=None):
    """
    Apply Johnson transformation to a parameter V. ( see Botto 2018)

    Parameters:
    - V: float or ndarray
        Original parameter(s) to transform.
    - transform_type: str
        One of 'LN', 'SB', or 'SU' indicating which Johnson transformation to use.
    - A, B: float
        Bounds required for the 'SB' transformation.

    Returns:
    - Y: float or ndarray
        Transformed variable (normally distributed).
    """
    transform_type = str(transform_type).replace('Johnson_', '')      # accept 'LN' and 'Johnson_LN' (as the inverse)
    rootName = ''
    if transform_type == rootName+'LN':
        return np.log(V)
    elif transform_type == rootName+'SB':
        if A is None or B is None:
            raise ValueError("A and B must be provided for the SB transformation.")
        U = (V - A) / (B - V)
        return np.log(U)
    elif transform_type == rootName+'SU':
        return np.arcsinh(V)
    else:
        raise ValueError("Invalid transform_type. Choose from 'LN', 'SB', or 'SU'.")

def inverse_Johnson1970_transform(Y, transform_type, A=None, B=None):
    """
    Inverse Johnson transformation from Y to V.

    Parameters:
    - Y: float or ndarray
        Transformed (normal) variable.
    - transform_type: str
        One of 'LN', 'SB', or 'SU'.
    - A, B: float
        Bounds required for 'SB' transformation.

    Returns:
    - V: float or ndarray
        Original variable value.
    """
    transform_type = str(transform_type).replace('Johnson_', '')
    if transform_type == 'LN':
        return np.exp(Y)
    elif transform_type == 'SB':
        if A is None or B is None:
            raise ValueError("A and B must be provided for the SB transformation.")
        U = np.exp(Y)
        return (A + B * U) / (1 + U)
    elif transform_type == 'SU':
        return np.sinh(Y)
    else:
        raise ValueError("Invalid transform_type. Choose from 'LN', 'SB', or 'SU'.")


def Evensen2003(q, wk, deltaT, gamma):
    """
    Ensemble Generation of time-variable atmospheric forcing rates.
    Parameters
    ----------
    wk : np.array([])
        is a sequence of white noise drawn from the standard normal distribution.
    deltaT : float
        assimilation interval in sec
    Tau : np.array([])
        the specified time decorrelation length
    time : int
        assimilation time index

    Returns
    -------
    qki : Ensemble of time-variable atmospheric forcing rates at time i

    """
    for i in range(1, len(deltaT)):
        # if Tau < deltaT[i]:
        #     raise ValueError(
        #         "Time decorrelation length is too small; should be at least>=" + str(deltaT)
        #     )
        # else:
            q[:, i] = (gamma[i] * q[:, i-1]) + (np.sqrt(1 - gamma[i] ** 2) * wk[:, i])

    return q


def build_dict_attributes_pert(
    var_per_2add, type_parm, parm_per_array, parm_sampling, **kwargs
):

    key = "ini_perturbation"
    var_per_2add[type_parm][key] = parm_per_array
    key = "sampling"
    var_per_2add[type_parm][key] = parm_sampling
    # Parameter tranformation
    # --------------------------------------------------------------------
    # The entry built from the scenario (create_dict_entry) already carries its
    # ``transf_type``: an explicit kwarg overrides it, but its absence must NOT reset it
    # to None (that silently disabled the log transform in the analysis step).
    if kwargs.get("transf_type") is not None:
        var_per_2add[type_parm]["transf_type"] = kwargs["transf_type"]
    else:
        var_per_2add[type_parm].setdefault("transf_type", None)
    if "transf_bounds" in kwargs:
        var_per_2add[type_parm]["transf_bounds"] = kwargs["transf_bounds"]
    elif (
        var_per_2add[type_parm].get("transf_type") == "log-ratio"
        and "transf_bounds" not in var_per_2add[type_parm]
    ):
        # bounded log needs A < p < B: default to the scenario bounds of the entry
        lo, hi = var_per_2add[type_parm].get("clip_min"), var_per_2add[type_parm].get("clip_max")
        if lo is None or hi is None:
            raise ValueError(
                f"{type_parm}: transf='log-ratio' needs both per_bounds (min and max) "
                "or an explicit transf_bounds={'A':..,'B':..}"
            )
        var_per_2add[type_parm]["transf_bounds"] = {"A": lo, "B": hi}
    # Parameter spatial extension
    # --------------------------------------------------------------------
    if "surf_zones_param" in kwargs:
        nb_surf_zones = kwargs["surf_zones_param"]
        # parm_per_array = np.tile(parm_per_array,nb_surf_zones)
        var_per_2add[type_parm]["surf_zones_param"] = kwargs["surf_zones_param"]
    return var_per_2add
        
def perturbate_parm_by_layers(type_parm,
                              ensemble_size,
                              nlayers,
                              mean=None,
                              sd=None,
                              per_type=None,
                              nominal=None,
                              clip_min=None,
                              clip_max=None,
                              sampler="random",
                              sampling_type="normal",
                              lognormal_center=None,
                              ):
            
        # sampler "random" (default): uniform draws from the global numpy state, so
        # np.random.seed(...) makes the prior reproducible. "sobol": the previous
        # low-discrepancy sampler, which has its own unseeded generator (and warns when
        # nlayers is not a power of 2).
        # Define the number of layers and the size of the ensemble
        # num_layers = 5

        # Create a Sobol sampler
        # mutable defaults are created per call (a shared default list is a latent bug)
        if mean is None:
            mean = []
        if sd is None:
            sd = []
        if nominal is None:
            nominal = []
        if sampler == "sobol":
            sobol_samples = qmc.Sobol(d=ensemble_size).random(n=nlayers)
        elif sampler == "random":
            # open interval: ppf(0) and ppf(1) are infinite
            sobol_samples = np.random.uniform(np.finfo(float).eps, 1 - np.finfo(float).eps,
                                              size=(nlayers, ensemble_size))
        else:
            raise ValueError(f"unknown sampler '{sampler}' (use 'random' or 'sobol')")
        
        # log-normal prior: the (truncated) normal law of ln(x); mean = the median, sd = std of ln(x)
        if sampling_type == "lognormal":
            lo = None if (clip_min is None or clip_min <= 0) else np.log(clip_min)
            hi = None if (clip_max is None or clip_max <= 0) else np.log(clip_max)
            lmean = lognormal_log_location(mean, sd, lognormal_center)
            if (lo is None and hi is None) or sd == 0:
                return np.exp(norm.ppf(sobol_samples, loc=lmean, scale=sd))
            a_, b_ = _std_bounds(lo, hi, lmean, sd)
            return np.exp(truncnorm.ppf(sobol_samples, a_, b_, loc=lmean, scale=sd))

        # Transform the Sobol sequence to follow a normal distribution
        if (clip_min is None and clip_max is None) or sd == 0:
            pressure_head_withLayers = norm.ppf(sobol_samples, loc=mean, scale=sd)
        else:
            # bounded prior: same low-discrepancy points through the truncated quantile function
            a, b = _std_bounds(clip_min, clip_max, mean, sd)
            pressure_head_withLayers = truncnorm.ppf(sobol_samples, a, b, loc=mean, scale=sd)
        
        parm_mat = np.ones(ensemble_size) * nominal
        if per_type == None:
            pressure_head_withLayers = pressure_head_withLayers
        if per_type == "multiplicative":
            pressure_head_withLayers = parm_mat * pressure_head_withLayers
        elif per_type == "additive":
            pressure_head_withLayers = parm_mat + pressure_head_withLayers

        return pressure_head_withLayers

        

# ---------------------------------------------------------------------------
# Correlated van Genuchten parameters (Botto et al. 2018), applied while the parameters are perturbed
# ---------------------------------------------------------------------------
VGP_CORRELATION = {"spec": None, "sigma_vg": None}      # spec None -> independent draws (default)


def set_vgp_correlation(spec=None, sigma_vg=None):
    """Choose how the VGP (SATCELL_VG, n_VG, thetar_VG) of one group are drawn by ``perturbate_parm``.

    ``spec``: None (default: independent draws), ``'botto:sandy_loam'`` | ``'botto:silt_loam'`` | ``{group: soil_class}``
    (Botto et al. 2018 Table 2; sigma_VG defaults to 0.1), or a 3 x 3 correlation matrix / ``{group: matrix}`` (order
    SATCELL_VG, n_VG, thetar_VG; ``sigma_vg`` defaults to 1). The draw is replaced once all three entries of a group exist."""
    VGP_CORRELATION["spec"], VGP_CORRELATION["sigma_vg"] = spec, sigma_vg


def _apply_vgp_correlation_if_ready(var_per):
    from pyCATHY.DA import vgp_correlation as VC
    spec, sig = VGP_CORRELATION["spec"], VGP_CORRELATION["sigma_vg"]
    groups = sorted({g for key in var_per for n in VC.VGP_NAMES if (g := VC._group_id(key, n)) is not None})
    for g in groups:
        keys = [f"{n}{g}" for n in VC.VGP_NAMES]
        if not all(k in var_per for k in keys) or any("vgp_corr_mode" in var_per[k] for k in keys):
            continue
        sub = {k: var_per[k] for k in keys}                  # same entry objects: modified in place
        if spec is None:
            mode = "independent"
            print(f"[VGP corr] group {g}: DEFAULT (independent draws, no correlation)")
        elif isinstance(spec, str) and spec.startswith("botto:") or (isinstance(spec, dict) and
                                                                      all(isinstance(v, str) for v in spec.values())):
            cls = spec.split(":", 1)[1] if isinstance(spec, str) else spec
            sv = VC.SIGMA_VG_UPDATED if sig in (None, 1.0) else sig
            _, rep = VC.apply_botto_vgp(sub, cls, sigma_vg=sv)
            mode = f"botto:{cls}, sigma_VG={sv}"
            for _g, r in rep.items():
                print(f"[VGP corr] group {g}: NOT default, Botto 2018 {r['soil_class']}, sigma_VG={sv}\n"
                      f"  target corr\n{np.round(r['target'], 2)}\n  sample corr\n{np.round(r['sample'], 2)}")
        else:
            _, rep = VC.apply_vgp_correlation(sub, spec, sigma_vg=1.0 if sig is None else sig)
            mode = "user matrix"
            for _g, r in rep.items():
                print(f"[VGP corr] group {g}: NOT default, user correlation matrix\n  target\n{np.round(r['target'], 2)}"
                      f"\n  sample\n{np.round(r['sample'], 2)}")
        for k in keys:
            var_per[k]["vgp_corr_mode"] = mode
    return var_per

def perturbate_parm(
    var_per,
    parm,
    type_parm,
    mean=None,
    sd=None,
    per_type=None,
    sampling_type="lognormal",
    ensemble_size=128,
    seed=True,
    **kwargs,
):
    """
    Perturbate parameter for ensemble generation

    Possible variable to perturbate:
        - initial conditions
        - hyetograph atmbc
        - van Genuchten retention curves parameters
        - Feddes parameters

    Perturbation:
        - parameters with constrainsts (truncated distribution)
        - parameters time dependant
        - other types of parameters

    Returns
    -------
    var_per : dict
        Perturbated variable dictionnary

    """
    # mutable defaults are created per call (a shared default list is a latent bug)
    if mean is None:
        mean = []
    if sd is None:
        sd = []

    savefig = False
    if 'savefig' in kwargs:
       savefig = kwargs["savefig"]
    
    if per_type=='None':
        per_type = None
    
    # type_parm2check = 
    # match_withLayers = re.search(rf'{type_parm}\d+', type_parm)
    # match_withLayers = re.search(r'\d', type_parm)
    nlayers = None
    if 'nlayers' in kwargs:
        nlayers = kwargs['nlayers']
    elif isinstance(parm, dict):
        nlayers = parm.get('nlayers')  # stored in the entry by perturbate()

    var_per_2add = {}
    
    # pert_control_name

    # copy initiail variable dict and add 'sampling' and 'ini_perturbation' attributes
    # -------------------------------------------------------------------------
    var_per_2add[type_parm] = parm

    key = "sampling_type"
    var_per_2add[type_parm][key] = sampling_type
    key = "sampling_mean"
    var_per_2add[type_parm][key] = mean
    key = "sampling_nominal"
    var_per_2add[type_parm][key] = parm['nominal']
    key = "sampling_sd"
    var_per_2add[type_parm][key] = sd
    key = "per_type"
    var_per_2add[type_parm][key] = per_type

    key_root = re.split(r"(\d+)", type_parm)
    parm_key_nb = 0
    if len(key_root)>1:
        parm_key_nb = int(key_root[1])
       
        
    # Contrainsted perturbation (bounded)
    # --------------------------------------------------------------------
    if "Archie" in type_parm:
        if parm.get('pert_control_name') in LAYER_LIKE:
            _require_nlayers(type_parm, nlayers)
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                        ensemble_size,
                                                        nlayers,
                                                        mean,
                                                        sd,
                                                        per_type=parm['per_type'],
                                                        nominal=parm['nominal'],
                                                        clip_min=parm.get('clip_min'),
                                                        clip_max=parm.get('clip_max')     
                                                        )
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
        else:
            parm_sampling, parm_per_array = Archie_pert_rules(
                parm, type_parm, ensemble_size, mean, sd, per_type, sampling_type
            )
    elif "porosity" in type_parm:
        if parm.get('pert_control_name') in LAYER_LIKE:
            _require_nlayers(type_parm, nlayers)
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                                  ensemble_size,
                                                                  nlayers,
                                                                  mean,
                                                                  sd,
                                                                  per_type=parm['per_type'],
                                                                  nominal=parm['nominal'],
                                                        clip_min=parm.get('clip_min'),
                                                        clip_max=parm.get('clip_max')     
                                                                  )
            key_root = re.split(r"(\d+)", type_parm)
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
    
        else:
            parm_sampling = sampling_dist_trunc(
                myclip_a=0, myclip_b=1, ensemble_size=ensemble_size, loc=mean, scale=sd
            )
            # parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
            
            if len(key_root)>1:
                # Case with different zones
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )
            else:
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )  
    
    
    
    elif "zroot".casefold() in type_parm.casefold():
        parm_sampling = sampling_dist_trunc(
            myclip_a=var_per_2add[type_parm]["clip_min"],
            myclip_b=var_per_2add[type_parm]["clip_max"],
            ensemble_size=ensemble_size,
            loc=mean,
            scale=sd,
        )
        parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)

    # check if parameters in part of van Genuchten retention curves
    # ----------------------------------------------------------------------
    # if type_parm in ['Alpha', 'nVG', 'thethaR']: #van Genuchten retention curves
    elif "VG" in type_parm:  # van Genuchten retention curves

        if parm.get('pert_control_name') in LAYER_LIKE:
            _require_nlayers(type_parm, nlayers)
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                        ensemble_size,
                                                        nlayers,
                                                        mean,
                                                        sd,
                                                        per_type=parm['per_type'],
                                                        nominal=parm['nominal'],
                                                        clip_min=parm.get('clip_min'),
                                                        clip_max=parm.get('clip_max')     
                                                        )
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
        else:
            parm_sampling, parm_per_array = VG_pert_rules(
                                                            var_per_2add,
                                                            parm,
                                                            type_parm,
                                                            ensemble_size,
                                                            mean,
                                                            sd,
                                                            per_type,
                                                            sampling_type,
                                                            **kwargs,
                                                        )



    # Time dependant perturbation
    # --------------------------------------------------------------------
    elif "atmbc" in type_parm:
    # elif any("atmbc" in item for item in type_parm):
        parm_sampling, parm_per_array = atmbc_pert_rules(    
                                                        var_per_2add,
                                                        parm,
                                                        type_parm,
                                                        ensemble_size,
                                                        mean,
                                                        sd,
                                                        per_type,
                                                        sampling_type,
                                                        **kwargs
                                                        )
        savefig = False

    elif "ic" in type_parm:

        if parm.get('pert_control_name') in LAYER_LIKE:
            _require_nlayers(type_parm, nlayers)
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                        ensemble_size,
                                                        nlayers,
                                                        mean,
                                                        sd,
                                                        per_type=parm['per_type'],
                                                        nominal=parm['nominal'],
                                                        clip_min=parm.get('clip_min'),
                                                        clip_max=parm.get('clip_max')
                                                        )
            key_root = re.split(r"(\d+)", type_parm)
    
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
        else:
            if _has_bounds(parm):  # bounded prior (scenario per_bounds)
                parm_sampling = sampling_dist_trunc(
                    myclip_a=parm["clip_min"],
                    myclip_b=parm["clip_max"],
                    ensemble_size=ensemble_size,
                    loc=parm['mean'],
                    scale=sd,
                )
            else:
                parm_sampling = sampling_dist(sampling_type,  parm['mean'], sd, ensemble_size,
                                              center=parm.get('lognormal_center'))
            
            if len(key_root)>1:
                # Case with different zones
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )
            else:
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )  
    
    # For all other types of perturbation
    # --------------------------------------------------------------------
    else:
        if parm.get('pert_control_name') in LAYER_LIKE:
            _require_nlayers(type_parm, nlayers)
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                        ensemble_size,
                                                        nlayers,
                                                        mean,
                                                        sd,
                                                        per_type=parm['per_type'],
                                                        nominal=parm['nominal'],
                                                        clip_min=parm.get('clip_min'),
                                                        clip_max=parm.get('clip_max'),
                                                        sampling_type=(sampling_type if sampling_type == 'lognormal' else 'normal'),
                                                        )
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
        else:
            if _has_bounds(var_per_2add[type_parm]):
                parm_sampling = sampling_dist_trunc(
                    myclip_a=var_per_2add[type_parm]["clip_min"],
                    myclip_b=var_per_2add[type_parm]["clip_max"],
                    ensemble_size=ensemble_size,
                    loc= parm['mean'],
                    scale=sd,
                )
            else:
                parm_sampling = sampling_dist(sampling_type,  parm['mean'], sd, ensemble_size,
                                              center=parm.get('lognormal_center'))
            parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    
                
            

    # build dictionnary of perturbated variable
    # --------------------------------------------------------------------
    var_per_2add = build_dict_attributes_pert(
                                                var_per_2add, 
                                                type_parm, 
                                                parm_per_array, 
                                                parm_sampling, 
                                                **kwargs
                                                )

    # Add to var perturbated stacked dict
    # ----------------------------------
    var_per = var_per | var_per_2add
    if type_parm in ("SATCELL_VG", "n_VG", "thetar_VG") or str(type_parm).startswith(("SATCELL_VG", "n_VG", "thetar_VG")):
        var_per = _apply_vgp_correlation_if_ready(var_per)

    if savefig is not False:
        plt_CT.plot_hist_perturbated_parm(
                                            parm, 
                                            var_per, 
                                            type_parm, 
                                            parm_per_array, 
                                            **kwargs
                                        )

    return var_per



