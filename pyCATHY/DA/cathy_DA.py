"""
Class managing Data Assimilation process i.e.:

    1. Parameters and data perturbation

        This step consists in the generation of the ensemble.

    Steps 2 and 3 below are controlled by the ``DA`` class:

    2.  Open Loop

       Without data run the hydrological model for an ensemble of realisations
       using the `_run_hydro_DA_openLoop` module.

    3. Loop of DA:

        - Running hydro model

        Run iteratively (between two observations) the hydrological model
        using the `_run_ensemble_hydrological_model` module.


        - Map states 2 Observations


        The module `map_states2Observations()` takes care of the mapping.

        For complex mapping such as ERT, it is done via a separate class.
        Example for the ERT we have a `class mappingERT()` which uses the
        Archie's petrophysical relationship.


        - Running the analysis

        The analysis is controlled by `run_analysis()` which takes as argument
        the type of analysis to run i.e. EnkF or Pf.

        - Update ensemble

        The ensemble file update is controlled by `apply_ensemble_updates()`


        - Evaluate performance


"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Callable, Iterator, List, Optional, Sequence, Tuple, Union
import contextlib
import functools
import glob

import multiprocessing

import logging
import os
import sys
import re
import shutil
import subprocess
import warnings
import copy

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd
from rich.progress import track

from pyCATHY.cathy_utils import dictObs_2pd  # noqa: F401  (re-exported, see __all__)
from pyCATHY.cathy_tools import CATHY, subprocess_run_multi
from pyCATHY.DA import enkf, pf, mapper, localisation
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as plt_CT

# Number of CPU cores left free when running the ensemble in parallel.
RESERVED_CPUS = 5

# Public API. `dictObs_2pd`, `mapper` and `localisation` are re-exported on
# purpose: `pyCATHY.DA.observations` and the documentation import them from here.
__all__ = [
    "DA",
    "DAConfig",
    "AnalysisConfig",
    "AnalysisResult",
    "run_analysis",
    "build_obs_node_mask",
    "dictObs_2pd",
    "mapper",
    "localisation",
]

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
) -> None:
    """Set ``value`` in the columns ``update_cols`` of the soil table ``df`` (in place).

    ``df`` is indexed by (zone, layer). Where the value goes depends on how the
    perturbation was built: ``"zone"`` (zone number ``index + 1``), ``"layers"``
    (layer ``index``), ``"root_map"`` (the zones of vegetation ``index + 1``,
    from ``veg_map`` and ``zone``); anything else, or no ``index``: every row.
    """
    if pert_control_name == "zone" and index is not None:
        df.loc[(index + 1, slice(None)), update_cols] = value
    elif pert_control_name == "layers" and index is not None:
        df.loc[(slice(None), index), update_cols] = value
    elif pert_control_name == "root_map" and index is not None:
        zones_in_vegi = zone[index + 1 == veg_map]
        df.loc[zones_in_vegi, update_cols] = value
    else:
        df.loc[:, update_cols] = value


def _log_domain_check(
    data_valid: np.ndarray,
    prediction_valid: np.ndarray,
    data_cov_valid: np.ndarray,
    use_log: bool = False,
) -> None:
    """Log the ranges of the observations, the predictions and the data variances.

    Diagnostic only (it used to be a nested function whose threshold checks
    were commented out): nothing is raised.
    """
    cov_diag = np.diag(data_cov_valid)
    logger.debug(f"\n{'LOG' if use_log else 'LINEAR'} domain check:")
    logger.debug(
        f"  Observation range: {data_valid.min():.2e} – {data_valid.max():.2e}"
    )
    logger.debug(
        f"  Prediction range:  {prediction_valid.min():.2e} – {prediction_valid.max():.2e}"
    )
    logger.debug(f"  Cov diag range:    {cov_diag.min():.2e} – {cov_diag.max():.2e}")


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
    if "rFluid".casefold() in key:
        return values < 0, "negative"
    if "porosity".casefold() in key:  # a volume fraction
        return (values < 0) | (values > 1), "out of range"
    if "a_Archie".casefold() in key:
        return (values < 0) | (values > 3), "out of range"
    if "Ks".casefold() in key:
        return values < 0, "negative Ks"
    if "PCREF".casefold() in key:  # a pressure head: must be negative
        return values > 0, "positive PCREF"
    if "Zroot".casefold() in key:
        return values > root_depth_limit(), "unfeasible root depth"
    return None


@dataclass
class _RunOptions:
    """Settings of one `run_DA_sequential` call, shared by all its time windows."""

    DA_type: str
    list_parm2update: Sequence[str]
    list_assimilated_obs: Union[str, Sequence[str]]
    update_key: str
    parallel: bool
    callexe: str
    threshold_rejected: float


@dataclass
class _EnsembleUpdate:
    """What `apply_ensemble_updates` shares between its keys and members.

    The lists are the ensembles of heterogeneous soil parameters (one entry per
    member), built on first use and then filled member after member.
    """

    cycle_nb: int
    update_key: str
    dict_parm_pert: dict
    zone3d: object
    ks: list = field(default_factory=list)
    porosity: list = field(default_factory=list)
    vg: list = field(default_factory=list)
    feddes: list = field(default_factory=list)


# Errors expected when a CATHY output file is missing, empty or malformed
# (pandas' EmptyDataError / ParserError are subclasses of ValueError).
_READ_ERRORS = (OSError, ValueError, KeyError, IndexError)

_WARNED_ONCE: set = set()

logger = logging.getLogger(__name__)


def _configure_logging(verbose: bool) -> None:
    """Set the verbosity of this module's logger (called when a DA run starts).

    INFO (progress, what is being assimilated) by default, DEBUG (diagnostic
    dumps of ranges, shapes, dictionaries) when ``verbose``; warnings and errors
    are always shown. A plain stdout handler is attached only when the
    application has not configured logging itself, so the messages still appear
    for scripts and notebooks that never call ``logging.basicConfig``.
    """
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    if not logger.handlers and not logging.getLogger().handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
        logger.propagate = False


def _warn_once(key: str, message: str) -> None:
    """Emit a warning only the first time `key` is seen (avoids one warning per analysis step)."""
    if key not in _WARNED_ONCE:
        _WARNED_ONCE.add(key)
        warnings.warn(message)


def _check_positive_pressure_heads(ensemble_psi):
    """Print a notice if any ensemble member has a positive pressure head."""
    if np.any(ensemble_psi > 0):
        logger.warning("!positive pressure heads observed!")


# ---------------------------------------------------------------------------
# ANALYSIS: one adapter per DA_type, selected through the ANALYSES registry
# ---------------------------------------------------------------------------
# An adapter has the signature
#     adapter(data, data_cov, param, state, prediction, config) -> AnalysisResult
# where `state` is the ensemble matrix to update (psi or sw) and `config` is an
# `AnalysisConfig` (Sakov, localisation matrix, inflation, ...).
# The adapters hide the differences between the filters' return values.


@dataclass
class AnalysisResult:
    """Common output of every analysis (EnKF variants and particle filter).

    Attributes
    ----------
    analysis : np.ndarray, shape (n_state, n_members)
        Updated state ensemble.
    analysis_param : np.ndarray, shape (n_members, n_param)
        Updated parameter ensemble (n_param may be 0 for a state-only update).
    diagnostics : dict or None
        EnKF matrices (A, Amean, dA, dD, MeasAvg, S, COV, B, dAS) used for the
        diagnostic plots; None for filters that do not provide them.
    extra : dict
        Filter specific outputs (particle filter: weights, n_eff, resampled, observation).
    """

    analysis: np.ndarray
    analysis_param: np.ndarray
    diagnostics: Optional[dict] = None
    extra: dict = field(default_factory=dict)


_ENKF_DIAGNOSTICS = ("A", "Amean", "dA", "dD", "MeasAvg", "S", "COV", "B", "dAS")


def _from_enkf(result, transpose_param=False):
    """Convert the list returned by an ``enkf.*`` filter into an AnalysisResult.

    The filters return [*diagnostics, analysis, analysis_param]. ``transpose_param``
    is needed for the filters that return the parameters as (n_param, n_members).
    """
    *diagnostics, analysis, analysis_param = result
    if transpose_param:
        analysis_param = analysis_param.T
    return AnalysisResult(
        analysis=analysis,
        analysis_param=analysis_param,
        diagnostics=dict(zip(_ENKF_DIAGNOSTICS, diagnostics)) or None,
    )


@dataclass
class AnalysisConfig:
    """Options of one analysis step (see ``run_analysis``).

    Attributes
    ----------
    sakov : bool
        Use the simplified (Sakov) covariance form where the filter supports it.
    inflation_states : float
        Multiplicative inflation of the ensemble. Legacy names: ``inflation``
        and ``alpha``.
    inflation_params : float
        Multiplicative inflation of the parameter ensemble (localized filter).
    jitter_params : float
        Std of the additive noise on the parameters (localized filter).
    localisation_matrix : np.ndarray or None
        Localisation matrix (n_state x n_obs) of the localized filter; an empty
        list/array means "no localisation" and is stored as None.
    pf_resample_threshold, pf_jitter_std_param, pf_jitter_std_state : float
        Particle filter settings.
    """

    sakov: bool = True
    inflation_states: float = 1.0
    inflation_params: float = 1.0
    jitter_params: float = 0.0
    localisation_matrix: Optional[np.ndarray] = None
    pf_resample_threshold: float = 0.5
    pf_jitter_std_param: float = 0.01
    pf_jitter_std_state: float = 0.05

    def __post_init__(self):
        if self.localisation_matrix is not None:
            self.localisation_matrix = np.asarray(self.localisation_matrix)
            if self.localisation_matrix.size == 0:
                self.localisation_matrix = None

    @classmethod
    def from_legacy_kwargs(cls, kwargs):
        """Build a config from the historical ``run_analysis(**kwargs)`` options.

        The recognised keys are popped from ``kwargs``; any key left over is an
        error for the caller to report.
        """
        renamed = {
            "Sakov": "sakov",
            "localisationMatrix": "localisation_matrix",
            "inflate_params": "inflation_params",
        }
        values = {new: kwargs.pop(old) for old, new in renamed.items() if old in kwargs}
        for name in (
            "jitter_params",
            "pf_resample_threshold",
            "pf_jitter_std_param",
            "pf_jitter_std_state",
        ):
            if name in kwargs:
                values[name] = kwargs.pop(name)
        # `inflation` and `alpha` are two historical names of the same option
        inflation, alpha = kwargs.pop("inflation", None), kwargs.pop("alpha", None)
        if inflation is not None or alpha is not None:
            values["inflation_states"] = inflation if inflation is not None else alpha
        return cls(**values)


# DA-only options historically passed as **kwargs to run_DA_sequential / prepare_DA
# (every other keyword, e.g. IPRT1 or TRAFLAG, belongs to the CATHY processor)
_LEGACY_DA_KEYS = (
    "inflation",
    "localisation",
    "threshold_rejected",
    "use_log_transformed_obs",
    "cleanup_vtk",
    "backup_fort777",
    "seed",
    "save_results",
)


@dataclass
class DAConfig:
    """Options of a sequential DA run (see ``DA.run_DA_sequential``).

    Attributes
    ----------
    inflation : float
        Multiplicative inflation passed to the analysis (``AnalysisConfig.inflation_states``).
    localisation : str or None
        Domain localisation: one of ``"veg_map"``, ``"zones"``, ``"nodes"``. A
        list with a single entry is accepted and converted to that string
        (``localisation.create_mask_localisation`` compares it to a string).
        Currently the observations must be defined on every surface node (maps
        such as ET or EM); point sensors fail in ``_DA_analysis`` (see README).
    threshold_rejected : int
        Ensemble members whose water-balance error exceeds it are rejected.
    use_log_transformed_obs : bool
        Assimilate log-transformed observations (useful for ERT data: the Archie
        operator makes their distribution non Gaussian).
    cleanup_vtk : {False, 'converted', 'all'}
        Delete VTK files produced in ``<workdir>/<project>/DA_Ensemble`` after the run.
    backup_fort777 : bool
        Keep a copy of the ``fort.777`` (ET) file of every DA cycle.
    analysis : AnalysisConfig or None
        Options of the analysis step (Sakov form, inflation of states and parameters,
        jitter, particle filter settings). None: defaults, with ``inflation`` of this
        config. When given it is used as is, except that the localisation matrix of
        the current cycle (``self.localisationMatrix``) is used when there is one.
    seed : int or None
        Seed of the global numpy random generator, set at the start of the run (the
        filters draw their random numbers from it). None: do not touch the generator.
    save_results : bool
        Write ``df_DA`` and ``ET_DA_xr`` to ``<workdir>/<project>/`` at every cycle,
        so that a crash does not lose the results of the finished cycles.
    verbose : bool
        Print the processor output of every member.
    """

    inflation: float = 1.0
    localisation: Optional[str] = None
    threshold_rejected: int = 10
    use_log_transformed_obs: bool = False
    cleanup_vtk: object = False
    backup_fort777: bool = False
    analysis: Optional[AnalysisConfig] = None
    seed: Optional[int] = None
    save_results: bool = False
    verbose: bool = False

    def __post_init__(self):
        if isinstance(self.localisation, (list, tuple)):
            if len(self.localisation) > 1:
                raise ValueError(
                    "Only one localisation type at a time is supported, got "
                    f"{list(self.localisation)}."
                )
            self.localisation = self.localisation[0] if self.localisation else None
        if self.localisation not in (None, "veg_map", "zones", "nodes"):
            raise ValueError(
                "localisation must be None, 'veg_map', 'zones' or 'nodes', got "
                f"{self.localisation!r}."
            )
        if self.cleanup_vtk not in (False, "converted", "all"):
            raise ValueError(
                f"cleanup_vtk must be False, 'converted' or 'all', got {self.cleanup_vtk!r}."
            )

    @classmethod
    def from_legacy_kwargs(cls, kwargs):
        """Build a config from the historical ``**kwargs`` of ``run_DA_sequential``.

        The DA-only keys are popped from ``kwargs`` (the rest belongs to the CATHY
        processor). ``verbose`` is read but left in ``kwargs`` because the processor
        uses it too.
        """
        values = {name: kwargs.pop(name) for name in _LEGACY_DA_KEYS if name in kwargs}
        if values:
            warnings.warn(
                f"Passing {sorted(values)} as keyword arguments is deprecated: "
                "use config=DAConfig(...).",
                DeprecationWarning,
                stacklevel=3,
            )
        if "verbose" in kwargs:
            values["verbose"] = kwargs["verbose"]
        return cls(**values)


def _enkf(filter_function, forward_alpha=False, **fixed):
    """Adapter for the filters called as f(data, cov, param, state, prediction, **kw)."""

    def adapter(data, data_cov, param, state, prediction, config):
        kwargs = dict(fixed)
        if forward_alpha:
            kwargs["alpha"] = config.inflation_states
        return _from_enkf(
            filter_function(data, data_cov, param, state, prediction, **kwargs)
        )

    return adapter


def _enkf_localized_with_inflation(data, data_cov, param, state, prediction, config):
    result = enkf.enkf_analysis_localized_with_inflation(
        data,
        data_cov,
        state,  # NB: ensemble before param in this filter's signature
        param,
        prediction,
        config.localisation_matrix,  # (state x obs), None = no localisation
        Sakov=config.sakov,
        inflate_states=config.inflation_states,
        inflate_params=config.inflation_params,
        jitter_params=config.jitter_params,
    )
    # this filter returns the parameters as (n_param, n_members), all others as (n_members, n_param)
    return _from_enkf(result, transpose_param=True)


def _particle_filter(data, data_cov, param, state, prediction, config):
    analysis = pf.particle_filter_analysis(
        data,
        data_cov,
        param,
        state,
        prediction,
        use_enkf_update=False,
        resample_threshold=config.pf_resample_threshold,
        jitter_std_param=config.pf_jitter_std_param,
        jitter_std_state=config.pf_jitter_std_state,
    )
    logger.info(f"Final N_eff: {analysis['n_eff']:.1f}")
    logger.info(f"Resampled: {analysis['resampled']}")
    return AnalysisResult(
        analysis=analysis["Analysis"],
        analysis_param=analysis["Analysisparam"].T,
        extra={
            k: analysis[k] for k in ("weights", "n_eff", "resampled", "observation")
        },
    )


ANALYSES = {
    "enkf_Evensen2009": _enkf(enkf.enkf_analysis),
    "enkf_Evensen2009_Sakov": _enkf(enkf.enkf_analysis, Sakov=True),
    "enkf_analysis_localized_with_inflation": _enkf_localized_with_inflation,
    "enkf_analysis_inflation": _enkf(enkf.enkf_analysis_inflation, forward_alpha=True),
    "enkf_analysis_inflation_multiparm": _enkf(
        enkf.enkf_analysis_inflation_multiparm, forward_alpha=True
    ),
    "enkf_dual_analysis": _enkf(enkf.enkf_dual_analysis, Sakov=False),
    "enkf_dual_analysis_Sakov": _enkf(enkf.enkf_dual_analysis, Sakov=True),
    "pf": _particle_filter,
}


def run_analysis(
    DA_type: str,
    data: np.ndarray,
    data_cov: np.ndarray,
    param: np.ndarray,
    list_update_parm: Sequence[str],
    ensembleX: Sequence[np.ndarray],
    prediction: np.ndarray,
    default_state: str = "psi",
    config: Optional[AnalysisConfig] = None,
    **legacy_kwargs: Any,
) -> AnalysisResult:
    """
    Perform the DA analysis step

    Parameters
    ----------
    DA_type : str
        type of analysis: one of the keys of ``ANALYSES``.
    data : np.array([])
        measured data.
    data_cov : np.array([]) or int
        measured data covariance matrice.
    param : np.array([])
        model parameters to update, shape (n_param, n_members); use an empty
        (0, n_members) array for a state-only update.
    list_update_parm : list
        names of the updated parameters (kept for API compatibility, unused here).
    ensembleX : list of np.array
        [pressure head, saturation water] ensemble matrices at each node.
    prediction : np.array([])
        predicted observation (after mapping).
    default_state : {"psi", "sw"}
        which element of ``ensembleX`` is updated.
    config : AnalysisConfig, optional
        analysis options (defaults if omitted).
    **legacy_kwargs
        Historical way of passing the options (Sakov, localisationMatrix,
        inflation / alpha, inflate_params, jitter_params, pf_*); converted with
        ``AnalysisConfig.from_legacy_kwargs``. Cannot be combined with ``config``.

    Returns
    -------
    AnalysisResult
    """
    if DA_type not in ANALYSES:
        raise ValueError(f"Unknown DA_type {DA_type!r}. Available: {sorted(ANALYSES)}")

    if config is None:
        config = AnalysisConfig.from_legacy_kwargs(legacy_kwargs)
    if legacy_kwargs:
        raise TypeError(
            "run_analysis got unexpected keyword argument(s) "
            f"{sorted(legacy_kwargs)} (or legacy options combined with `config`)."
        )

    state = ensembleX[1 if default_state == "sw" else 0]
    return ANALYSES[DA_type](data, data_cov, param, state, prediction, config)


# ---------------------------------------------------------------------------
# SMALL HELPERS (pure functions, no `self`)
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def _working_directory(path: str) -> Iterator[None]:
    """Change the working directory inside the ``with`` block, restore it afterwards (even on error).

    CATHY files are addressed relative to the directory of an ensemble member,
    and some pyCATHY routines read ``os.getcwd()``; this is the only way
    ``cathy_DA`` changes the directory of the process.
    """
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _restores_cwd(method: Callable) -> Callable:
    """Decorator: the working directory is restored when the method returns."""

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        with _working_directory(os.getcwd()):
            return method(*args, **kwargs)

    return wrapper


def build_obs_node_mask(
    surface_nodes_mask: np.ndarray,
    mesh_nodes_observed: np.ndarray,
    zone_nodes: np.ndarray,
    n_repeat: Optional[int] = None,
) -> np.ndarray:
    """Boolean mask over the observation vector selecting the observations of one zone.

    Parameters
    ----------
    surface_nodes_mask : np.ndarray of bool, shape (n_surface_nodes,)
        True for the surface nodes of the zone.
    mesh_nodes_observed : np.ndarray of int, shape (n_obs,)
        Mesh node of every observation.
    zone_nodes : np.ndarray of int
        Mesh nodes (all depths) of the zone.
    n_repeat : int, optional
        Number of consecutive sets of surface nodes in the observation vector (EM:
        one set per coil); None for a single set.

    Returns
    -------
    np.ndarray of bool, shape (n_obs,)
        When the observations are one value per surface node (maps) the surface mask
        is used as is; otherwise (point sensors) an observation belongs to the zone if
        its mesh node is one of ``zone_nodes``.
    """
    mask = surface_nodes_mask
    if n_repeat is not None:
        mask = np.tile(surface_nodes_mask, n_repeat)
    if len(mesh_nodes_observed) != len(mask):
        mask = np.isin(mesh_nodes_observed, zone_nodes)
    return mask


# observation types whose predictions may be negative (pressure heads)
_NEGATIVE_VALUES_ALLOWED = ("tensio",)


def _update_key(cycle_nb: int) -> str:
    """Key of `dict_parm_pert[parameter]` holding the ensemble of a DA cycle.

    Cycle 0 (and below) is the initial perturbation, cycle n > 0 the n-th update.
    """
    return "ini_perturbation" if cycle_nb <= 0 else f"update_nb{cycle_nb}"


# ---------------------------------------------------------------------------
# PERFORMANCE: RMSE / NRMSE bookkeeping (pure functions, no `self`)
# ---------------------------------------------------------------------------


def _compute_diff_matrix(data, prediction):
    """Absolute difference between data (n_obs,) and prediction (n_obs, n_ens)."""
    data = np.atleast_1d(np.asarray(data, dtype=float))  # () -> (1,)
    prediction = np.asarray(prediction, dtype=float)
    if prediction.ndim == 1:  # (n_ens,) -> (1, n_ens)
        prediction = prediction.reshape(len(data), -1)
    if data.shape[0] != prediction.shape[0]:  # never broadcast one datum over all rows
        raise ValueError(
            f"{data.shape[0]} observations but {prediction.shape[0]} predictions"
        )
    return np.abs(data[:, None] - prediction)  # broadcast, no loops


def _compute_rmse(diff_matrix):
    """Mean absolute difference: over ensemble members (axis 1), then over observations."""
    return float(np.mean(np.nanmean(diff_matrix, axis=1)))


def _get_rmse_for_sensor(obs2eval, prediction):
    """Computes RMSE for a given sensor."""
    diff_mat = _compute_diff_matrix(obs2eval, prediction)
    return _compute_rmse(diff_mat)


def _append_to_performance_df(
    df_performance,
    t_obs,
    name_sensor,
    rmse_sensor,
    rmse_avg,
    nrmse_sensor,
    nrmse_avg,
    ol_bool,
):
    """Appends performance metrics to the main DataFrame.

    The wide layout and the ``NMRMSE`` column names are legacy: they are read by
    ``pyCATHY.plotters.cathy_plots.DA_RMS`` and by existing parquet files.
    Change them only together with that plotter.
    """
    data = {
        "time": [t_obs],
        "ObsType": [name_sensor],
        "RMSE" + name_sensor: [rmse_sensor],
        "RMSE_avg": [rmse_avg],
        "NMRMSE" + name_sensor: [nrmse_sensor],
        "NMRMSE_avg": [nrmse_avg],
        "OL": [ol_bool],
    }
    return pd.concat([df_performance, pd.DataFrame(data)], ignore_index=True)


def _save_performance_to_parquet(df_performance, file_path="performance_data.parquet"):
    """Saves the performance DataFrame to Parquet format, appending if the file exists."""
    if os.path.exists(file_path):
        df_existing = pd.read_parquet(file_path, engine="pyarrow")
        df_performance = pd.concat([df_existing, df_performance], ignore_index=True)

    df_performance.to_parquet(
        file_path, engine="pyarrow", compression="snappy", index=False
    )


def _load_performance_from_parquet(file_path="performance_data.parquet"):
    """Loads performance data from a Parquet file if it exists."""
    if os.path.exists(file_path):
        return pd.read_parquet(file_path, engine="pyarrow")
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no file exists


# %%
# DA class
# ---------


class DA(CATHY):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._init_da_state()

    def _init_da_state(self) -> None:
        """Defaults of the attributes that belong to the DA (the rest comes from `CATHY`).

        They used to appear whenever some method first assigned them, and the code
        guarded the reads with `getattr(self, name, None)`. A run resets the
        per-run ones in `_DA_init`; the others (`Archie_parms`, `VGN_parms`,
        `localisationMatrix`, `ens_valid`) may be set by the user before the run.
        """
        # set by the user (or by `set_Archie_parm`) before the run, else defaulted in the run
        self.Archie_parms: Optional[dict] = None
        self.VGN_parms: Optional[dict] = None
        self.localisationMatrix: Optional[Any] = (
            None  # legacy: one localisation matrix per DA cycle
        )
        self.ens_valid: Optional[Sequence[int]] = (
            None  # indices of the valid members (all of them by default)
        )
        # per-run state (reset by `_DA_init`). `dict_obs` and `count_DA_cycle` are NOT
        # defaulted here: `CATHY.__init__` defines them, and CATHY tests
        # `count_DA_cycle is not None` to know whether a DA is running.
        self.dict_parm_pert: dict = {}
        self.df_DA = pd.DataFrame()
        self.df_Archie = pd.DataFrame()
        self.ET_DA_xr = None
        self.isET = False
        self.localisationBool = False
        self.Archie_parms_mat_ens: Any = []  # list, then an array (n_param, n_members)
        self.count_atmbc_cycle = 0

    def _DA_init(self, dict_obs, dict_parm_pert, list_parm2update="all"):
        """
        Initial preparation for DA


        .. note::

            0. check dictionnaries validity
            1. update cathyH file given NENS, ENS_times + flag self.parm['DAFLAG']==True
            2. create the processor cathy_origin.exe
            3. create directory tree for the ensemble
            4. create panda dataframe for each perturbated variable
            5. update input files

        Parameters
        ----------
        NENS : int
            # Number of realizations in the ensemble
        ENS_times : int
            # Observation times

        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        # Initiate
        # -------------------------------------------------------------------
        update_key = "ini_perturbation"

        # check if dict_obs is ordered in time
        # ------------------------------------
        if (
            all(i < j for i, j in zip(list(dict_obs.keys()), list(dict_obs.keys())[1:]))
        ) is False:
            raise ValueError("Observation List is not sorted.")

        # if hasattr(self,'dict_obs') is False:
        self.dict_obs = dict_obs  # self.dict_obs is already assigned (in read observatio! change it to self.obs
        self.dict_parm_pert = dict_parm_pert
        self.df_DA = pd.DataFrame()
        self.df_Archie = pd.DataFrame()
        # per-run state, reset together with df_DA (it used to be created lazily and
        # an ET_DA_xr of a previous run was silently extended by the next run)
        self.ET_DA_xr = None
        self.isET = False
        self.localisationBool = False

        # Infer ensemble size NENS from perturbated parameter dictionnary
        # -------------------------------------------------------------------
        for name in self.dict_parm_pert:
            self.NENS = len(self.dict_parm_pert[name]["ini_perturbation"])

        # Infer ensemble update times ENS_times from observation dictionnary
        # -------------------------------------------------------------------
        ENS_times = []
        for ti in self.dict_obs:
            ENS_times.append(float(ti))

        # start DA cycle counter
        # -------------------------------------------------------------------
        self.count_DA_cycle = 0
        self.count_atmbc_cycle = 0
        # (the counter is incremented during the update analysis)

        if self.ens_valid is None:
            self.ens_valid = list(np.arange(0, self.NENS))

        # create sub directories for each ensemble
        # ---------------------------------------------------------------------
        self._create_subfolders_ensemble()

        # update list of updated parameters based on problem heterogeneity
        # ---------------------------------------------------------------------
        newlist_parm2update = ["St. var."]
        for d in self.dict_parm_pert.keys():
            for l in list_parm2update:
                if l in d:
                    newlist_parm2update.append(d)

        return ENS_times, newlist_parm2update, update_key

    def prepare_DA(
        self,
        dict_obs,
        dict_parm_pert,
        list_parm2update,
        list_assimilated_obs,
        parallel,
        config=None,
        **kwargs,
    ):

        if config is None:
            config = DAConfig.from_legacy_kwargs(kwargs)
        elif any(key in kwargs for key in _LEGACY_DA_KEYS):
            raise TypeError(
                "DA options given both as `config` and as keyword arguments."
            )
        self.config = config
        if config.seed is not None:
            np.random.seed(config.seed)

        self.run_processor(
            DAFLAG=1, runProcess=False, **kwargs
        )  # to create the processor exe
        callexe = "./" + self.processor_name

        self.inflation = config.inflation
        self.localisation = config.localisation
        threshold_rejected = config.threshold_rejected
        self.verbose = config.verbose
        _configure_logging(config.verbose)
        self.use_log_transformed_obs = config.use_log_transformed_obs
        # backup all spatial ET iterations (the ET files are large!)
        self.backupfort777 = config.backup_fort777

        # list_parm2update
        # initiate DA
        # -------------------------------------------------------------------
        ENS_times, list_parm2update, update_key = self._DA_init(
            dict_obs,
            dict_parm_pert,
            list_parm2update,
        )
        # initiate mapping petro
        # -------------------------------------------------------------------
        self._mapping_petro_init()

        # update the ensemble files with pertubarted parameters
        # -------------------------------------------------------------------
        self.apply_ensemble_updates(
            dict_parm_pert,
            list_parm2update="all",
            cycle_nb=self.count_DA_cycle,
        )

        all_atmbc_times = np.copy(self.atmbc["time"])

        return (callexe, all_atmbc_times, threshold_rejected, ENS_times, update_key)

    def run_DA_smooth(
        self,
        DA_type,
        dict_obs,
        list_parm2update,
        dict_parm_pert,
        list_assimilated_obs,
        parallel=True,
        **kwargs,
    ):
        """

        Run Data Assimilation

        .. note::

            Steps are:
            1. DA init (create subfolders)
            2. run CATHY hydrological model stacked using ALL DA times
            3. check before analysis
            4. analysis
            5. check after analysis

        Parameters
        ----------
        parallel : bool
            True for multiple ensemble folder simulations.
        DA_type : str
            type of Data Assimilation.
        list_assimilated_obs : TYPE
            list of names of observations to assimilate.
        dict_obs : dict
            observations dictionnary.
        list_update_parm : list
            list of names of parameters to update.
        dict_parm_pert : dict
            parameters perturbated dictionnary.
        """

        (
            callexe,
            self.all_atmbc_times,
            threshold_rejected,
            self.ENS_times,
            update_key,
        ) = self.prepare_DA(
            dict_obs,
            dict_parm_pert,
            list_parm2update,
            list_assimilated_obs,
            parallel,
            **kwargs,
        )

        # map states to observation = apply H operator to state variable
        # ----------------------------------------------------------------
        prediction = self.map_states2Observations(
            list_assimilated_obs,
            default_state="psi",
            parallel=parallel,
        )
        # self.atmbc

        # DA analysis
        # ----------------------------------------------------------------
        self.console.print(r""":chart_with_upwards_trend: [b]Analysis[/b]:
                               - DA type: {}
                               - Inflation: {}
                           """.format(DA_type, self.inflation))

        (
            ensemble_psi_valid_bef_update,
            ensemble_sw_valid_bef_update,
            ens_size,
            sim_size,
        ) = self._read_state_ensemble()

        (
            ensemble_psi,
            ensemble_sw,
            data,
            analysis,
            analysis_param,
        ) = self._DA_analysis(
            prediction,
            DA_type,
            list_parm2update,
            list_assimilated_obs,
            ens_valid=self.ens_valid,
        )

        pass

    def get_mesh_mask_nodes_observed(self, list_assimilated_obs):
        mesh_mask_nodes_observed = np.zeros(int(self.grid3d["nnod3"]), dtype=bool)
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        mesh_nodes_observed = []
        if obskey2map[0] == "ERT":
            raise ValueError("Cannot localise with apparent resistivity data.")
        elif obskey2map[0] == "EM":
            nbofcoils = len(obs2map[0]["coils"])
            for node_nb in obs2map[0]["data"]["mesh_node"]:
                mesh_mask_nodes_observed[int(node_nb)] = True
                mesh_nodes_observed.append(int(node_nb))
            mesh_mask_nodes_observed = [mesh_mask_nodes_observed] * nbofcoils
            mesh_nodes_observed = [mesh_nodes_observed] * nbofcoils
        else:
            for obsi in obs2map:
                mesh_mask_nodes_observed[obsi["mesh_nodes"]] = True
                mesh_nodes_observed.append(obsi["mesh_nodes"])
        return mesh_nodes_observed, mesh_mask_nodes_observed

    @_restores_cwd
    def run_DA_sequential(
        self,
        DA_type,
        dict_obs,
        list_parm2update,
        dict_parm_pert,
        list_assimilated_obs,
        parallel=True,
        config=None,
        **kwargs,
    ):
        """

        Run Data Assimilation

        .. note::

            Steps are:
            1. DA init (create subfolders)
            2. run CATHY hydrological model recursively using DA times
                <- Loop-->
                    3. check before analysis
                    4. analysis
                    5. check after analysis
                    6. update input files
                <- Loop-->

        Parameters
        ----------
        parallel : bool
            True for multiple ensemble folder simulations.
        DA_type : str
            type of Data Assimilation.
        list_assimilated_obs : TYPE
            list of names of observations to assimilate.
        dict_obs : dict
            observations dictionnary.
        list_update_parm : list
            list of names of parameters to update.
        dict_parm_pert : dict
            parameters perturbated dictionnary.
        config : DAConfig, optional
            DA options (inflation, localisation, threshold_rejected,
            use_log_transformed_obs, cleanup_vtk, verbose). The same options are
            still accepted as keyword arguments (deprecated).
        **kwargs
            Passed to the CATHY processor (IPRT1, TRAFLAG, VTKF, ...).

        """
        if config is None:
            config = DAConfig.from_legacy_kwargs(kwargs)

        (
            callexe,
            self.all_atmbc_times,
            threshold_rejected,
            self.ENS_times,
            update_key,
        ) = self.prepare_DA(
            dict_obs,
            dict_parm_pert,
            list_parm2update,
            list_assimilated_obs,
            parallel,
            config=config,
            **kwargs,
        )

        options = _RunOptions(
            DA_type,
            list_parm2update,
            list_assimilated_obs,
            update_key,
            parallel,
            callexe,
            threshold_rejected,
        )

        # -------------------------------------------------------------------
        # (time-windowed) update input files ensemble again
        # -------------------------------------------------------------------
        self._update_input_ensemble(list_parm2update)
        # -----------------------------------
        # Run hydrological model sequentially
        # = Loop over atmbc times (including assimilation observation times)
        # -----------------------------------
        # the callees below restore the working directory they change: the loop
        # always runs from the project working directory
        os.chdir(self.workdir)
        for t_atmbc in self.all_atmbc_times[:-1]:  # N times -> N-1 windows
            self._run_window(t_atmbc, options)

        self._cleanup_vtk_files(config.cleanup_vtk)

    def _run_window(self, t_atmbc: float, options: _RunOptions) -> None:
        """One time window of the sequential DA: run the ensemble, analyse, update the inputs.

        1. run CATHY for the window ending at the next atmbc time; 2. drop the members
        whose run failed; 3. analysis if there is an observation at ``t_atmbc``;
        4. record the results; 5. update the input files for the next window.
        """
        self._run_ensemble_hydrological_model(options.parallel, options.callexe)

        # check scenario (result of the hydro simulation)
        # ----------------------------------------------------------------
        rejected_ens = self._check_before_analysis(
            self.ens_valid,
            options.threshold_rejected,
        )
        id_valid = ~np.array(rejected_ens)
        self.ens_valid = list(np.arange(0, self.NENS)[id_valid])
        # check if there is an observation at the given atmbc time
        # --------------------------------------------------------
        if t_atmbc in self.ENS_times:
            state_bef_update, analysis_valid = self._analysis_cycle(
                options.DA_type,
                options.list_parm2update,
                options.list_assimilated_obs,
                options.update_key,
                options.parallel,
            )
        else:
            state_bef_update, analysis_valid = self._no_analysis_cycle()

        self.count_atmbc_cycle = self.count_atmbc_cycle + 1

        self._log_end_of_time_step(t_atmbc)

        # create dataframe _DA_var_pert_df holding the results of the DA update
        # ---------------------------------------------------------------------
        self._save_step_outputs(
            state=state_bef_update,
            state_analysis=analysis_valid,
            rejected_ens=rejected_ens,
            list_assimilated_obs=options.list_assimilated_obs,
            list_parm2update=options.list_parm2update,
        )

        logger.info(f"Parameters to update are: {options.list_parm2update}")
        # overwrite input files ensemble (perturbated variables)
        # ---------------------------------------------------------------------
        if self.count_atmbc_cycle < (
            len(self.all_atmbc_times) - 1
        ):  # -1 cause all_atmbc_times include TMAX
            self._update_input_ensemble(
                options.list_parm2update, analysis=analysis_valid
            )
        else:
            self.console.rule(":red_circle: end of DA" ":red_circle:", style="yellow")

        self._log_end_of_da_update()
        plt.close("all")

    def _analysis_cycle(
        self,
        DA_type: str,
        list_parm2update: Sequence[str],
        list_assimilated_obs: Union[str, Sequence[str]],
        update_key: str,
        parallel: bool,
    ) -> Tuple[list, np.ndarray]:
        """One DA cycle at an observation time: map, analyse, check, record.

        Returns
        -------
        (list, array)
            ``[psi, sw]`` of the ensemble before the update, and the valid
            analysed ensemble (states, plus parameters when updated).
        """
        self.console.print(":world_map: [b] map states to observations [/b]")
        # H operator applied to the state variable
        prediction = self.map_states2Observations(
            list_assimilated_obs,
            default_state="psi",
            parallel=parallel,
        )
        self._validate_prediction(prediction, list_assimilated_obs)

        self.console.print(r""":chart_with_upwards_trend: [b]Analysis[/b]:
                               - DA type: {}
                               - Inflation: {}
                           """.format(DA_type, self.inflation))
        psi_bef_update, sw_bef_update, _, _ = self._read_state_ensemble()

        # global or domain-localised analysis
        ensemble_psi, ensemble_sw, data, analysis, analysis_param = self._analyse(
            prediction, DA_type, list_parm2update, list_assimilated_obs
        )

        (
            prediction_valid,
            ensemble_psi_valid,
            _,
            analysis_valid,
            analysis_param_valid,
        ) = self._mark_invalid_ensemble(
            self.ens_valid,
            prediction,
            ensemble_psi,
            ensemble_sw,
            analysis,
            analysis_param,
        )

        self.console.print(":face_with_monocle: [b]check analysis performance[/b]")
        self._performance_assessment_pq(
            list_assimilated_obs,
            data,
            prediction_valid,
            t_obs=self.count_DA_cycle,
        )
        self.count_DA_cycle = self.count_DA_cycle + 1

        _check_positive_pressure_heads(ensemble_psi_valid)
        self.update_pert_parm_dict(update_key, list_parm2update, analysis_param_valid)
        return [psi_bef_update, sw_bef_update], analysis_valid

    def _no_analysis_cycle(self) -> Tuple[list, np.ndarray]:
        """Time step without observation: only the hydrological model was run.

        Returns the same pair as `_analysis_cycle`; the ensemble is not
        validated and the "analysis" is the unchanged pressure-head ensemble.
        """
        self.console.print(
            ":confused: No observation for this time - run hydrological model only"
        )
        psi, sw, _, _ = self._read_state_ensemble()
        return [psi, sw], psi

    def _validate_prediction(self, prediction, list_assimilated_obs=None):
        """Raise ValueError if the predicted observations contain invalid values.

        NaN is always invalid. Observation types in ``_NEGATIVE_VALUES_ALLOWED``
        (pressure heads of tensiometers) accept any finite value; zero is accepted
        for ET (``self.isET``); every other observation type must be strictly positive.
        The rule is applied row by row when there is one row per observation.
        """
        prediction = np.asarray(prediction, dtype=float)
        nan = np.isnan(prediction)

        allow_negative = np.zeros(prediction.shape[0], dtype=bool)
        if list_assimilated_obs is not None:
            _, obs2map = self._obs_key_select(list_assimilated_obs)
            types = [obs.get("data_type") for obs in obs2map]
            if len(types) == prediction.shape[0]:
                allow_negative = np.array(
                    [t in _NEGATIVE_VALUES_ALLOWED for t in types]
                )
            elif types and all(t in _NEGATIVE_VALUES_ALLOWED for t in types):
                allow_negative[:] = True
        allow_negative = allow_negative.reshape((-1,) + (1,) * (prediction.ndim - 1))

        non_positive = prediction < 0 if self.isET else prediction <= 0
        non_positive = non_positive & ~allow_negative
        invalid = non_positive | nan

        if np.any(invalid):
            logger.warning(
                f"⚠️ Invalid prediction values detected: {np.sum(invalid)} entries total"
            )
            logger.warning(f"   - NaN values: {np.sum(nan)}")
            logger.warning(
                f"   - Negative or zero values: {np.sum(non_positive & ~nan)}"
            )
            raise ValueError(
                "prediction contains invalid values. Please clean or replace them before proceeding."
            )

    def _log_end_of_time_step(self, t_atmbc):
        """Console progress after a time step of the loop."""
        self.console.rule(
            ":round_pushpin: end of time step (s) "
            + str(int(t_atmbc))
            + "/"
            + str(int(self.all_atmbc_times[-1]))
            + " :round_pushpin:",
            style="yellow",
        )
        self.console.rule(
            ":round_pushpin: end of atmbc update # "
            + str(self.count_atmbc_cycle)
            + "/"
            + str(len(self.all_atmbc_times))
            + " :round_pushpin:",
            style="yellow",
        )

    def _log_end_of_da_update(self):
        """Console progress after a DA update."""
        self.console.rule(
            ":red_circle: end of DA update "
            + str(self.count_DA_cycle)
            + "/"
            + str(len(self.ENS_times))
            + " :red_circle:",
            style="yellow",
        )
        self.console.rule(
            ":dart: % of valid ensemble: "
            + str((len(self.ens_valid) * 100) / (self.NENS))
            + ":dart:",
            style="yellow",
        )

    def _save_step_outputs(
        self,
        state,
        state_analysis,
        rejected_ens,
        list_assimilated_obs,
        list_parm2update,
    ):
        """Store the results of the current step (dataframes, ET, backups)."""
        self.console.rule("Back up DA step")
        self._DA_df(
            state=state,
            state_analysis=state_analysis,
            rejected_ens=rejected_ens,
        )
        self._DA_ET_xr()

        # export summary results of DA
        meta_DA = {
            "listAssimilatedObs": list_assimilated_obs,
            "listUpdatedparm": list_parm2update,
        }
        # backup ET map for each iteration of ensemble numbers
        if self.backupfort777:
            self.backup_fort777()
        self.backup_results_DA(meta_DA)
        self.backup_simu()
        if self.config.save_results:
            self._save_results_to_disk()

    def _save_results_to_disk(self):
        """Write the in-memory results of the DA (``df_DA``, ``ET_DA_xr``) next to the project."""
        out_dir = os.path.join(self.workdir, self.project_name)
        self.df_DA.to_pickle(os.path.join(out_dir, self.project_name + "_df_DA.pkl"))
        if self.ET_DA_xr is not None:
            self.ET_DA_xr.to_netcdf(
                os.path.join(out_dir, self.project_name + "_ET_DA.nc")
            )

    def _cleanup_vtk_files(self, mode):
        """Delete the VTK files of the ensemble members (``mode``: False, 'converted', 'all')."""
        if not mode:
            return
        ensemble_root = os.path.join(self.workdir, self.project_name, "DA_Ensemble")
        patterns = {
            "converted": os.path.join(ensemble_root, "**", "*converted*.vtk"),
            "all": os.path.join(ensemble_root, "**", "*.vtk"),
        }
        for f in glob.glob(patterns[mode], recursive=True):
            os.remove(f)

    # backup ET map file at each DA iteration
    def backup_fort777(self):

        for ensi in range(self.NENS):
            dst_dir = os.path.join(
                self.workdir,
                self.project_name,
                "DA_Ensemble/cathy_" + str(ensi + 1),
                "DAnb" + str(self.count_DA_cycle) + "fort.777",
            )

            shutil.copy(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ensi + 1),
                    "fort.777",
                ),
                dst_dir,
            )

    def _analyse(self, prediction, DA_type, list_parm2update, list_assimilated_obs):
        """Run the analysis of one DA cycle, globally or zone by zone (domain localisation).

        Returns
        -------
        ensemble_psi, ensemble_sw, data, analysis, analysis_param
            as returned by ``_DA_analysis``; with a localisation the zone analyses are
            assembled on the whole mesh (nodes without observation keep their prior).
        """
        if self.localisation is not None:
            self.localisationBool = self.localisation in (
                "veg_map",
                "zones",
                "nodes",
            )

        if self.localisationBool:
            return self._analyse_by_zones(
                prediction, DA_type, list_parm2update, list_assimilated_obs
            )
        return self._analyse_global(
            prediction, DA_type, list_parm2update, list_assimilated_obs
        )

    def _analyse_global(
        self, prediction, DA_type, list_parm2update, list_assimilated_obs
    ):
        """Analysis without localisation: every node and every observation is valid."""
        obs_valid = np.ones(prediction.shape[0], dtype=bool)
        mesh3d_nodes_valid = np.ones(int(self.grid3d["nnod3"]), dtype=bool)
        return self._DA_analysis(
            prediction,
            DA_type,
            list_parm2update,
            list_assimilated_obs,
            ens_valid=self.ens_valid,
            obs_valid=obs_valid,
            obs_mesh_nodes_valid=mesh3d_nodes_valid,
        )

    def _analyse_by_zones(
        self, prediction, DA_type, list_parm2update, list_assimilated_obs
    ):
        """Domain localisation: one analysis per zone, assembled on the whole mesh."""
        self.console.print(
            f":carpentry_saw: [b] Domain localisation: {self.localisation} [/b]"
        )

        if getattr(self, "zone", None) is None:
            zone, _ = self.read_inputs("zone")
            self.zone = zone
        mesh3d_nodes_valid, mask_surface_nodes = localisation.create_mask_localisation(
            self.localisation,
            self.veg_map,
            self.zone,
            self.hapin,
            self.grid3d,
            self.DEM,
        )

        # position of the sensors, to make sure each one is located in the right zone
        mesh_nodes_observed, mesh_mask_nodes_observed = (
            self.get_mesh_mask_nodes_observed(list_assimilated_obs)
        )
        mesh_nodes_observed = np.hstack(mesh_nodes_observed)
        # EM: one set of surface nodes per coil
        n_repeat = (
            np.shape(mesh_mask_nodes_observed)[0]
            if "EM" in list_assimilated_obs
            else None
        )

        analysis_with_localisation = []
        analysis_param = []
        last_zone_ensemble = None
        # Loop over unique map nodes and perform the DA analysis
        for loci, node_value in enumerate(np.unique(mask_surface_nodes)):
            if node_value == -9999:
                continue

            # surface nodes of this zone
            obs_surface_nodes_bool = np.ones(int(self.grid3d["nnod"]), dtype=bool)
            valid_surface_nodes = mask_surface_nodes.ravel()
            obs_surface_nodes_bool[valid_surface_nodes != node_value] = False

            # observations of this zone
            obs_surface_nodes_valid = build_obs_node_mask(
                obs_surface_nodes_bool,
                mesh_nodes_observed,
                mesh3d_nodes_valid[loci],
                n_repeat=n_repeat,
            )
            # if None of the observation are within the localisation area:
            #     continue without analysis
            if np.sum(obs_surface_nodes_valid) == 0:
                continue

            selected_parm2update = self._zone_parameters(list_parm2update, loci)
            self.console.print(f"""[b]
                               SHORTCUT: ALL PARAMS ARE UPDATED but
                               are not necessaraly associated
                               with observation localisation [/b]
                               """)
            self.console.print(
                f"[b] Params to update (domain loc): {selected_parm2update} [/b]"
            )

            # Perform DA analysis for each mask of localisation
            (
                ensemble_psi_loci,
                ensemble_sw_loci,
                _,
                analysis_loci,
                analysis_param_loci,
            ) = self._DA_analysis(
                prediction,
                DA_type,
                selected_parm2update,
                list_assimilated_obs,
                ens_valid=self.ens_valid,
                obs_valid=obs_surface_nodes_valid,
                obs_mesh_nodes_valid=mesh3d_nodes_valid[loci],
            )
            analysis_with_localisation.append(
                np.c_[mesh3d_nodes_valid[loci], analysis_loci]
            )
            analysis_param.append(analysis_param_loci)
            last_zone_ensemble = (ensemble_psi_loci, ensemble_sw_loci)

        if last_zone_ensemble is None:
            raise ValueError(
                f"no observation lies in any '{self.localisation}' localisation zone: "
                "nothing to analyse"
            )

        ensemble_psi, ensemble_sw = last_zone_ensemble
        analysis = self._assemble_localised_analysis(
            np.vstack(analysis_with_localisation), ensemble_psi
        )
        # the performance is assessed on ALL the observations (as without
        # localisation), not on those of one zone: the predictions cover all of them
        data, _ = self._get_data2assimilate(list_assimilated_obs)
        return (
            ensemble_psi,
            ensemble_sw,
            np.asarray(data, dtype=float),
            analysis,
            np.hstack(analysis_param),
        )

    def _zone_parameters(self, list_parm2update, loci):
        """Parameters updated by the analysis of zone number ``loci``.

        With a vegetation-map or zone localisation each zone updates the state
        and its own parameter (the one at position ``loci + 1``).
        """
        if self.localisation in ("veg_map", "zones"):
            if len(list_parm2update) > 1:
                return [list_parm2update[0], list_parm2update[loci + 1]]
            return [list_parm2update[0]]
        return list_parm2update

    def _assemble_localised_analysis(self, analysis_with_localisation, ensemble_psi):
        """Analysis on the whole mesh from the stacked zone analyses.

        ``analysis_with_localisation``: rows ``[node index, analysed values...]``.
        The nodes that no zone analysed keep the state read from the ensemble.
        """
        analysis_with_localisation = analysis_with_localisation[
            analysis_with_localisation[:, 0].argsort()
        ]
        if len(analysis_with_localisation) == len(ensemble_psi):
            return analysis_with_localisation[:, 1:]

        self.console.print(f"""
                               [b] !! All the mesh nodes havent been analysed with DA !! [/b]
                               Taking state values where missing analysis nodes
                           """)
        analysis, _, _, _ = self._read_state_ensemble()
        indice_nodes_with_analysis = [
            int(ii) for ii in analysis_with_localisation[:, 0]
        ]
        analysis[indice_nodes_with_analysis] = analysis_with_localisation[:, 1:]
        return analysis

    def _DA_analysis(
        self,
        prediction,
        DA_type="enkf_Evensen2009_Sakov",
        list_update_parm=None,
        list_assimilated_obs="all",
        ens_valid=None,
        obs_valid=None,
        obs_mesh_nodes_valid=None,
    ):
        """
        Analysis ensemble using DA

        1. map state variable 2 Observations

        2. apply filter

        3. checkings

        Parameters
        ----------
        typ : str
            type of Data Assimilation
        NUDN : int
            #  NUDN  = number of observation points for nudging or EnKF (NUDN=0 for no nudging)
        ENS_times : np.array([])
            # Observation time for ensemble kalman filter (EnKF).
        rejected_ens : list
            list of ensemble to reject
        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """
        if list_update_parm is None:
            list_update_parm = ["St. var."]
        if ens_valid is None:
            ens_valid = []
        if obs_valid is None:
            obs_valid = []
        if obs_mesh_nodes_valid is None:
            obs_mesh_nodes_valid = []

        update_key = _update_key(self.count_DA_cycle)

        obskey2map, _ = self._obs_key_select(list_assimilated_obs)
        ensemble_psi, ensemble_sw, _, _ = self._read_state_ensemble()

        # select ONLY the data to assimilate
        data, _ = self._get_data2assimilate(list_assimilated_obs)
        self.console.print(r""":
                               - Data size: {}
                                 --> Observations --> {}
                           """.format(str(np.shape(data)), list_assimilated_obs))
        data_cov = self._data_cov_of_cycle(len(data))
        analysis_config = self._analysis_config_of_cycle()

        # filter non-valid ensemble members before the analysis
        # ---------------------------------------------------------------------
        ensemble_psi_valid = ensemble_psi[:, ens_valid]
        ensemble_sw_valid = ensemble_sw[:, ens_valid]

        # prepare parameters
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updating
        # the parameters
        param = self.transform_parameters(list_update_parm, update_key)
        logger.debug("parm size: " + str(len(param)))

        # empty (0, n_members) array when only the states are updated: the inflation
        # variants of the filters fail on a plain empty list
        if len(param) > 0:
            param_valid = param[:, ens_valid]
        else:
            param_valid = np.zeros((0, len(ens_valid)))

        # filter non-valid observations before the analysis
        # (particularly for the DA with localisation)
        # ---------------------------------------------------------------------
        data = np.asarray(data, dtype=float)
        data_valid = data[obs_valid]
        data_cov_valid = data_cov[obs_valid, :][:, obs_valid]
        prediction_valid = prediction[obs_valid, :]
        ensemble_psi_valid = ensemble_psi_valid[obs_mesh_nodes_valid, :]
        ensemble_sw_valid = ensemble_sw_valid[obs_mesh_nodes_valid, :]

        # observations and predictions as seen by the filter
        if self.use_log_transformed_obs:
            data_in, prediction_in = self._log_transform_obs(
                data_valid, prediction_valid, data_cov_valid
            )
        else:
            self._save_obs_pred_boxplot(data_valid, prediction_valid)
            data_in, prediction_in = data_valid, prediction_valid

        result_analysis = run_analysis(
            DA_type,
            data_in,
            data_cov_valid,
            param_valid,
            list_update_parm,
            [ensemble_psi_valid, ensemble_sw_valid],
            prediction_in,
            config=analysis_config,
        )
        analysis = result_analysis.analysis
        analysis_param = result_analysis.analysis_param

        # plot ensemble covariance matrices and changes (only for ENKF)
        if result_analysis.diagnostics is not None:
            self._plot_analysis_diagnostics(
                result_analysis.diagnostics,
                ensemble_psi_valid,
                data_valid,
                analysis,
                obskey2map,
            )

        # Back-transformation of the parameters
        # ---------------------------------------------------------------------
        self.transform_parameters(
            list_update_parm, param=np.hstack(analysis_param), back=True
        )
        return (ensemble_psi, ensemble_sw, data_valid, analysis, analysis_param)

    def _data_cov_of_cycle(self, n_data: int) -> np.ndarray:
        """Stacked data covariance of the current DA cycle, checked against the data size."""
        data_cov = self.stacked_data_cov[self.count_DA_cycle]
        actual_shape = np.shape(data_cov)
        if actual_shape[0] != n_data:
            raise ValueError(
                f"Wrong stacked_data_cov shape: expected ({n_data}, {n_data}), "
                f"but got {actual_shape}. "
                "Consider re-computing data covariance."
            )
        return data_cov

    def _analysis_config_of_cycle(self) -> AnalysisConfig:
        """`AnalysisConfig` of the current cycle: the user's one (or the default built from
        the inflation) completed with the localisation matrix of the cycle."""
        if self.localisationMatrix is not None:
            localisation_matrix = self.localisationMatrix[self.count_DA_cycle]
        else:
            localisation_matrix = None
        if self.config.analysis is None:
            return AnalysisConfig(
                inflation_states=self.inflation,
                localisation_matrix=localisation_matrix,
            )
        return replace(
            self.config.analysis,
            localisation_matrix=(
                localisation_matrix
                if localisation_matrix is not None
                else self.config.analysis.localisation_matrix
            ),
        )

    def _log_transform_obs(
        self,
        data_valid: np.ndarray,
        prediction_valid: np.ndarray,
        data_cov_valid: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """log10 of the observations and of the predictions (e.g. apparent resistivities).

        Logs the value ranges before and after, warns about non-positive
        predictions and saves a histogram observed vs predicted. Returns
        ``(data_log, prediction_log)``.
        """
        eps = 1e-3
        logger.debug("Before log-transform:")
        logger.debug(
            f"Observed ERa: min = {data_valid.min():.2f}, max = {data_valid.max():.2f}"
        )
        logger.debug(
            f"Predicted ERa: min = {prediction_valid.min():.2f}, max = {prediction_valid.max():.2f}"
        )

        invalid = (prediction_valid <= 0) | np.isnan(prediction_valid)
        if np.any(invalid):
            logger.warning(f"""
                      ⚠️ Invalid values found in prediction valid (before log-transform): 
                          {np.sum(prediction_valid <= 0)} entries ≤ 0 
                          {np.sum(np.isnan(prediction_valid))}  nan"
                        """)

        data_valid_log = np.log10(data_valid + eps)
        prediction_valid_log = np.log10(prediction_valid + eps)

        logger.debug("\nAfter log-transform:")
        logger.debug(
            f"Observed log10(ERa): min = {data_valid_log.min():.2f}, max = {data_valid_log.max():.2f}"
        )
        logger.debug(
            f"Predicted log10(ERa): min = {prediction_valid_log.min():.2f}, max = {prediction_valid_log.max():.2f}"
        )
        _log_domain_check(data_valid, prediction_valid, data_cov_valid, use_log=True)

        fig, ax = plt.subplots()
        ax.hist(
            np.log10(data_valid + 1e-3),
            bins=30,
            alpha=0.5,
            label="Observed",
            density=True,
        )
        ax.hist(
            np.log10(prediction_valid + 1e-3).flatten(),
            bins=30,
            alpha=0.5,
            label="Predicted",
            density=True,
        )
        ax.legend()
        plt.title("log10(ERa): Observed vs Predicted")
        fig.savefig(self._da_ensemble_file("log10_Obs_Pred_hist") + ".png", dpi=300)
        plt.close(fig)
        return data_valid_log, prediction_valid_log

    def _save_obs_pred_boxplot(
        self, data_valid: np.ndarray, prediction_valid: np.ndarray
    ) -> None:
        """Save a boxplot observed vs predicted in the `DA_Ensemble` folder."""
        fig, ax = plt.subplots()
        ax.boxplot(
            [data_valid.flatten(), prediction_valid.flatten()],
            labels=["Observed", "Predicted"],
            showfliers=True,  # Show outliers
            notch=True,  # Optional: adds a notch for median confidence
            patch_artist=True,  # Optional: for filled boxes
        )
        ax.set_title("Observed vs Predicted")
        fig.savefig(self._da_ensemble_file("Obs_Pred_hist") + ".png", dpi=300)
        plt.close(fig)

    def _da_ensemble_file(self, prefix: str) -> str:
        """Path (without extension) of a file of the current cycle in `DA_Ensemble`."""
        return os.path.join(
            self.workdir,
            self.project_name,
            "DA_Ensemble",
            prefix + str(self.count_DA_cycle),
        )

    def _plot_analysis_diagnostics(
        self,
        diag: dict,
        ensemble_psi_valid: np.ndarray,
        data_valid: np.ndarray,
        analysis: np.ndarray,
        obskey2map: Sequence[str],
    ) -> None:
        """Plot the ensemble covariance matrices and their change (diagnostic only)."""
        self.console.print("[b]Plotting COV matrices[/b]")
        try:
            if np.shape(diag["COV"])[0] < 1e9:
                plt_CT.show_DA_process_ens(
                    ensemble_psi_valid,
                    data_valid,
                    diag["COV"],
                    diag["dD"],
                    diag["dAS"],
                    diag["B"],
                    analysis,
                    savefig=True,
                    savename=self._da_ensemble_file("DA_Matrices_t"),
                    label_sensor=str([*obskey2map]),
                )
            else:
                self.console.print("[b]Impossible to plot matrices[/b]")
        except Exception as exc:  # plotting is diagnostic only: never abort the DA run
            self.console.print(
                f"[b]Impossible to plot matrices[/b] ({type(exc).__name__}: {exc})"
            )

    def _mark_invalid_ensemble(
        self, ens_valid, prediction, ensemble_psi, ensemble_sw, analysis, analysis_param
    ):
        """mark invalid ensemble - invalid ensemble are filled with NaN values"""

        ensemble_psi_valid = np.empty(ensemble_psi.shape)
        ensemble_psi_valid[:] = np.nan
        ensemble_psi_valid[:, ens_valid] = ensemble_psi[:, ens_valid]

        ensemble_sw_valid = np.empty(ensemble_psi.shape)
        ensemble_sw_valid[:] = np.nan
        ensemble_sw_valid[:, ens_valid] = ensemble_sw[:, ens_valid]

        analysis_valid = np.empty(ensemble_psi.shape)
        analysis_valid[:] = np.nan
        analysis_valid[:, ens_valid] = analysis

        prediction_valid = np.empty([prediction.shape[0], ensemble_psi.shape[1]])
        prediction_valid[:] = np.nan
        prediction_valid[:, ens_valid] = prediction

        analysis_param_valid = []
        if len(analysis_param[0]) > 0:
            analysis_param_valid = np.empty(
                [analysis_param.shape[1], ensemble_psi.shape[1]]
            )
            analysis_param_valid[:] = np.nan
            analysis_param_valid[:, ens_valid] = analysis_param.T

        return (
            prediction_valid,
            ensemble_psi_valid,
            ensemble_sw_valid,
            analysis_valid,
            analysis_param_valid,
        )

    def _run_ensemble_hydrological_model(self, parallel, callexe):
        """multi run CATHY hydrological model from the independant folders composing the ensemble"""

        self.console.print(":athletic_shoe: [b]Run hydrological model[/b]")
        # ----------------------------------------------------------------
        if parallel == True:
            pathexe_list = []
            # for ens_i in range(self.NENS):

            for ens_i in self.ens_valid:
                path_exe = os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
                pathexe_list.append(path_exe)
            with multiprocessing.Pool(
                processes=max(1, multiprocessing.cpu_count() - RESERVED_CPUS)
            ) as pool:
                result = pool.map(subprocess_run_multi, pathexe_list)

                if self.verbose:
                    self.console.print(result)
        # process each ensemble folder one by one
        # ----------------------------------------------------------------
        else:
            # Loop over ensemble realisations
            # ------------------------------------------------------------
            for ens_i in track(
                range(self.NENS), description="Run hydrological fwd model..."
            ):

                self.console.print(
                    ":keycap_number_sign: [b]ensemble nb:[/b]"
                    + str(ens_i + 1)
                    + "/"
                    + str(self.NENS)
                )
                member_dir = os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
                p = subprocess.run(
                    [callexe], text=True, capture_output=True, cwd=member_dir
                )

                if self.verbose:
                    self.console.print(p.stdout)
                    self.console.print(p.stderr)

        pass

    def _check_before_analysis(
        self, ens_valid: Optional[Sequence[int]] = None, threshold_rejected: float = 10
    ) -> List[bool]:
        """
        Filter is applied only on selected ensemble

        Parameters
        ----------
        update_key : TYPE
            DESCRIPTION.

        Returns
        -------
        rejected_ens

        """
        if ens_valid is None:
            ens_valid = []

        self.console.print(":white_check_mark: [b]check scenarii before analysis[/b]")

        ###CHECK which scenarios are OK and which are to be rejected and then create and applied the filter
        ###only on the selected scenarios which are to be considered.

        rejected_ens_new = []
        for n in range(self.NENS):
            if n not in ens_valid:
                rejected_ens_new.append(True)
                continue
            try:
                reason = self._member_rejection_reason(n)
            except _READ_ERRORS as exc:
                rejected_ens_new.append(True)
                logger.warning(
                    f"cannot read mbeconv for ensemble {n + 1} ({type(exc).__name__}: {exc})"
                )
                continue
            if reason is not None:
                self.console.print(f":x: [b]{reason}[/b], ens_nb:{n + 1}")
            rejected_ens_new.append(reason is not None)

        self._raise_if_too_many_rejected(rejected_ens_new, threshold_rejected)
        return rejected_ens_new

    def _member_rejection_reason(self, n: int) -> Optional[str]:
        """Why the CATHY run of member ``n`` (0-based) cannot be used, else ``None``.

        A run is rejected when its mass-balance file (`mbeconv`) is empty, has
        missing cumulative values, or stops before the last printing time of
        the current window. Raises one of `_READ_ERRORS` if it cannot be read.
        """
        df_mbeconv = out_CT.read_mbeconv(
            os.path.join(
                self.workdir,
                self.project_name,
                "DA_Ensemble/cathy_" + str(n + 1),
                "output/mbeconv",
            )
        )
        if len(df_mbeconv) == 0:
            return "len(df_mbeconv)['TIME']" + str(len(df_mbeconv))
        if df_mbeconv["CUM."].isnull().values.any():
            return "df_mbeconv['CUM.']" + str(df_mbeconv["CUM."].isnull().values.any())
        last_time = df_mbeconv["TIME"].iloc[-1]
        if np.round(last_time) < (
            np.round(self.parm["(TIMPRT(I),I=1,NPRT)"][-1]) - self.parm["DELTAT"]
        ):
            return "df_mbeconv['time']" + str(last_time)
        return None

    def _raise_if_too_many_rejected(
        self, rejected: Sequence[bool], threshold_rejected: float
    ) -> None:
        """Stop the run when the share of rejected members reaches the threshold (%)."""
        if sum(rejected) >= (threshold_rejected / 100) * self.NENS:
            logger.error("new parameters (if updated)")
            logger.error("%s", self.dict_parm_pert)
            raise ValueError(
                "% number of rejected ensemble is too high:"
                + str((sum(rejected) * 100) / (self.NENS))
            )

    def _check_after_analysis(
        self, update_key: str, list_update_parm: Sequence[str]
    ) -> set:
        """
        Check which scenarios parameters are OK and which one to discard after analysis.
        Test the range of each physical properties.

        Parameters
        ----------
        update_key : str
            Update cycle key identifier.
        list_update_parm : list
            list of updated model parameters (after analysis).
        """

        self.console.print(":white_check_mark: [b]check scenarii post update[/b]")

        id_nonvalid = set()
        for name in list_update_parm:
            if "St. var." in name:
                continue
            values = np.asarray(self.dict_parm_pert[name][update_key])
            rule = _nonphysical_values(name, values, lambda: self._max_root_depth())
            if rule is None:  # no physical bound known for this parameter
                continue
            mask, label = rule
            bad = np.where(mask)[0]
            if len(bad):
                self.console.print(
                    f":x: [b]{label} {name}:[/b]{values[bad]}, ens_nb:{list(bad)}"
                )
            id_nonvalid.update(bad.tolist())
        return set(np.arange(self.NENS)).difference(id_nonvalid)

    def _max_root_depth(self) -> float:
        """Deepest feasible root depth: the thickness of the 3d mesh."""
        depth = abs(min(self.grid3d["mesh3d_nodes"][:, -1]))
        logger.debug("deepest feasible root depth: %s", depth)
        return depth

    def transform_parameters(
        self, list_update_parm, update_key=None, back=False, param=None
    ):
        """
        Parameter (back) transformation before analysis (log, bounded log, ...)

        ..note:: parameters are log transform so that they become Gaussian distributed
                 to be updated by the Ensemble Kalman Filter
        """
        param_new = []
        for pp in enumerate(list_update_parm[:]):
            param_trans = []

            if "St. var." in pp[1]:
                pass
            else:
                if update_key:
                    param = self.dict_parm_pert[pp[1]][update_key]
                if self.dict_parm_pert[pp[1]]["transf_type"] == "log":
                    logger.debug("log transformation")
                    if back:
                        param_trans = np.exp(param)
                    else:
                        param_trans = np.log(param)
                elif self.dict_parm_pert[pp[1]]["transf_type"] == "log-ratio":
                    logger.debug("bounded log transformation")
                    range_tr = self.dict_parm_pert[pp[1]]["transf_bounds"]
                    A = range_tr["A"]  # min range
                    B = range_tr["B"]  # max range
                    if back:
                        param_trans = np.exp((param - A) / (B - param))
                    else:
                        param_trans = np.log((param - A) / (B - param))
                elif self.dict_parm_pert[pp[1]]["transf_type"] == "hyperbolic":
                    logger.debug("hyperbolic transformation")
                    # Y = sinh−1(U )
                    param_trans = np.log(self.dict_parm_pert[pp[1]][update_key])

                else:
                    param_trans = param
                param_new.append(param_trans)
            if len(param_new) > 0:
                np.vstack(param_new)

        return np.array(param_new)

    def _mapping_petro_init(self):
        """
        Initiate Archie and VGP petro/pedophysical parameters

        If Archie and VGP dictionnaries are not existing fill with default values
        """

        warnings_petro = []
        # check that porosity is a list
        # -------------------------------------

        if getattr(self, "soil_SPP", None) is None:
            self.read_inputs("soil")

        porosity = self.soil_SPP["SPP"]["POROS"].mean()
        if not isinstance(porosity, list):
            porosity = [porosity]

        # check if Archie parameters exists
        # ---------------------------------
        if self.Archie_parms is None:
            warnings_petro.append(
                "Archie parameters not defined - Falling back to defaults"
            )
            self.Archie_parms = {
                "porosity": porosity,
                "rFluid_Archie": [1.0],
                "a_Archie": [1.0],
                "m_Archie": [2.0],
                "n_Archie": [2.0],
                "pert_sigma_Archie": [0],
            }

        # check if VGN parameters exists
        # ---------------------------------
        if self.VGN_parms is None:
            warnings_petro.append(
                "VGN parameters not defined - Falling back to defaults"
            )
            self.VGN_parms = {}

        self.console.rule(":warning: warning messages above :warning:", style="yellow")
        for ww in warnings_petro:
            self.console.print(ww, style="yellow")
        self.console.rule("", style="yellow")

        pass

    def _update_soil_params_from_perturbations(
        self,
        dict_parm_pert: dict,
        df_SPP_2fill: List[pd.DataFrame],
        key_root: List[str],
        key: str,
        update_key: str,
        ens_nb: int,
        shellprint_update: bool,
    ) -> List[pd.DataFrame]:
        """Set one soil parameter of one ensemble member and, once complete, write its soil file.

        The perturbed value of member ``ens_nb`` is written into the member's soil
        table (one table per member, in ``df_SPP_2fill``) where ``pert_control_name``
        says (``zone``, ``layers``, ``root_map``, else everywhere). The soil file is
        only rewritten when no value is missing any more: with several zones or layers
        the file is written by the call that fills the last one.

        Parameters
        ----------
        dict_parm_pert : dict
            Perturbed parameters; ``dict_parm_pert[key]`` holds ``pert_control_name``
            and the values under ``update_key``.
        df_SPP_2fill : list of pandas.DataFrame
            Soil tables of the members, filled in place. Empty at the first call:
            it is then built from the soil files with the updated columns set to NaN.
        key_root : list
            Parameter name and zone / layer number, e.g. ``["Ks", "0"]``.
        key : str
            Key of the parameter in ``dict_parm_pert``, e.g. ``"Ks0"``.
        update_key : str
            Which ensemble of values to use (e.g. ``"update_nb1"``).
        ens_nb : int
            Member index (zero-based).
        shellprint_update : bool
            Echo the update of the soil file on the console.

        Returns
        -------
        list of pandas.DataFrame
            ``df_SPP_2fill``, for the next call.

        Raises
        ------
        ValueError
            If the parameter name is unknown.

        Must run inside the directory of the member (the soil file and the vtk
        are addressed from the current directory).
        """
        pert_control_name = dict_parm_pert[key]["pert_control_name"]
        update_value = dict_parm_pert[key][update_key][ens_nb]
        param_name = key_root[0]
        index = int(key_root[1]) if len(key_root) > 1 else None

        update_cols = _SOIL_PARM_COLUMNS.get(param_name)
        if update_cols is None:
            raise ValueError(f"Unknown parameter '{param_name}'.")

        if len(df_SPP_2fill) == 0:
            df_SPP_2fill.extend(self._soil_tables_to_fill(update_cols))
            logger.info(
                "\nAssuming that the dictionary of perturbed parameters is built "
                "using 'pert_control_name', i.e. per layer or per zone like ks0, ks1, etc...\n"
            )

        df = df_SPP_2fill[ens_nb]
        _assign_soil_value(
            df,
            pert_control_name,
            index,
            update_cols,
            update_value,
            self.veg_map if pert_control_name == "root_map" else None,
            self.zone if pert_control_name == "root_map" else None,
        )

        # the soil file is only written when nothing is missing any more
        if not df.isna().any().any():
            self._write_member_soil(df, pert_control_name, shellprint_update)

        return df_SPP_2fill

    def _soil_tables_to_fill(self, update_cols: List[str]) -> List[pd.DataFrame]:
        """Soil table of every member, with the columns to update emptied (NaN)."""
        tables = []
        for i in range(self.NENS):
            df_SPP, _ = in_CT.read_soil(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    f"DA_Ensemble/cathy_{i+1}",
                    "input/soil",
                ),
                self.dem_parameters,
                self.cathyH["MAXVEG"],
            )
            df = df_SPP.copy()
            df[update_cols] = np.nan
            tables.append(df)
        return tables

    def _write_member_soil(
        self, df: pd.DataFrame, pert_control_name: str, shellprint_update: bool
    ) -> None:
        """Write the soil file of the member (current directory); for layers, also map the porosity."""
        self.update_soil(
            SPP_map=df,
            FP_map=self.soil_FP["FP_map"],
            verbose=self.verbose,
            filename=os.path.join(os.getcwd(), "input/soil"),
            shellprint_update=shellprint_update,
            backup=True,
        )
        if pert_control_name == "layers":
            saveMeshPath = os.path.join(os.getcwd(), "vtk/", self.project_name + ".vtk")
            _, POROS_mesh_nodes = self.map_prop_2mesh_markers(
                "POROS",
                df["POROS"].groupby("layer").mean().to_list(),
                to_nodes=False,
                saveMeshPath=saveMeshPath,
            )
            self.map_prop2mesh({"POROS": POROS_mesh_nodes})

    def update_archie_parameter(
        self,
        Archie_p_names,
        key_root,
        update_key,
        ens_nb,
        dict_parm_pert,
        key,
    ):
        """
        Update a single Archie parameter for a given ensemble member.
        """
        self.console.print(":arrows_counterclockwise: [b]Update Archie parameters[/b]")

        idArchie = Archie_p_names.index(key_root[0])

        if len(self.Archie_parms_mat_ens) == 0:
            for p in Archie_p_names:
                if len(self.Archie_parms[p]) != self.NENS:
                    self.Archie_parms[p] = list(
                        self.Archie_parms[p] * np.ones(self.NENS)
                    )
            self.Archie_parms_mat_ens = np.zeros([len(Archie_p_names), self.NENS])
            self.Archie_parms_mat_ens[:] = np.nan

        self.Archie_parms_mat_ens[idArchie, ens_nb] = dict_parm_pert[key][update_key][
            ens_nb
        ]

        if np.isnan(self.Archie_parms_mat_ens[idArchie, :]).any() == False:
            self.Archie_parms[key_root[0]] = list(
                self.Archie_parms_mat_ens[idArchie, :]
            )
            logger.debug(f"New Archie parm: {self.Archie_parms[key_root[0]]}ˇ")

    def _update_atmbc_member(self, dict_parm_pert, ens_nb, key):
        """Update the atmbc of one ensemble member (only the first perturbed node updates the stacked parameters)."""
        if key.casefold() in "atmbc0".casefold():
            times = dict_parm_pert[key]["data2perturbate"]["time"]
            VALUE = dict_parm_pert[key]["time_variable_perturbation"][ens_nb]
            self.update_atmbc(
                HSPATM=self.atmbc["HSPATM"],
                IETO=self.atmbc["IETO"],
                time=list(times),
                netValue=VALUE,
                filename=os.path.join(os.getcwd(), "input/atmbc"),
                verbose=self.verbose,
            )

        else:
            pass

    def _update_wtposition_member(
        self, cycle_nb, dict_parm_pert, ens_nb, key, shellprint_update
    ):
        """Update the initial water table position of one ensemble member (initial perturbation only)."""
        if cycle_nb == 0:
            self.update_ic(
                INDP=4,
                IPOND=0,
                WTPOSITION=dict_parm_pert[key]["ini_perturbation"][ens_nb],
                filename=os.path.join(os.getcwd(), "input/ic"),
                backup=False,
                shellprint_update=shellprint_update,
            )

    def _update_ic_member(
        self, cycle_nb, dict_parm_pert, ens_nb, key, shellprint_update
    ):
        """Update the initial pressure head of one ensemble member (cycle 0 only, or by layers)."""
        if dict_parm_pert[key]["pert_control_name"] == "layers":

            # if match_ic_withLayers:
            DApath = os.getcwd()
            saveMeshPath = os.path.join(DApath, "vtk", self.project_name + ".vtk")
            ic_values_by_layers = dict_parm_pert[key][f"ini_{key}_withLayers"][
                :, ens_nb
            ]
            ic_nodes = self.map_prop_2mesh_markers(
                "ic",
                ic_values_by_layers,
                to_nodes=True,
                saveMeshPath=saveMeshPath,
            )
            self.update_ic(
                INDP=1,
                IPOND=0,
                pressure_head_ini=ic_nodes,
                filename=os.path.join(DApath, "input/ic"),
                backup=False,
                shellprint_update=shellprint_update,
            )
        else:
            if cycle_nb == 0:
                self.update_ic(
                    INDP=0,
                    IPOND=0,
                    pressure_head_ini=dict_parm_pert[key]["ini_perturbation"][ens_nb],
                    filename=os.path.join(os.getcwd(), "input/ic"),
                    backup=False,
                    shellprint_update=shellprint_update,
                )

    def apply_ensemble_updates(self, dict_parm_pert, list_parm2update, **kwargs):
        """
        Applies updates to ensemble files by overwriting input files.

        Updates may include:
        - State variables (e.g., pressure head from data assimilation)
        - Model parameters: Ks, porosity, Feddes, Archie, initial conditions, atmbc
        - Boundary and initial conditions

        Parameters
        ----------
        dict_parm_pert : dict
            Dictionary of perturbation data per parameter.
        list_parm2update : list or "all"
            List of parameters to update. Use "all" to update everything.
        **kwargs : dict
            cycle_nb: int, optional
                Data assimilation cycle number (0 = prior).
            analysis: np.ndarray, optional
                State analysis results to update the pressure head.

        Returns
        -------
        None
        """
        self.console.print(":arrows_counterclockwise: [b]Update ensemble[/b]")

        # replace dictionnary key with cycle nb
        # -------------------------------------
        cycle_nb = kwargs.get("cycle_nb", 0)
        analysis = kwargs.get("analysis", None)

        update_key = _update_key(cycle_nb)

        if cycle_nb > 0 and "St. var." in list_parm2update and analysis is None:
            raise ValueError(
                "`analysis` must be provided to update state variables when cycle_nb > 0."
            )

        if list_parm2update == "all":
            list_parm2update = ["St. var."]
            for pp in dict_parm_pert:
                list_parm2update.append(pp)
        self._update_states_from_analysis(list_parm2update, cycle_nb, analysis)

        # accumulators shared by all keys and members of this call
        self.Archie_parms_mat_ens = []  # matrice of ensemble for Archie parameters
        ctx = _EnsembleUpdate(
            cycle_nb=cycle_nb,
            update_key=update_key,
            dict_parm_pert=dict_parm_pert,
            zone3d=getattr(self, "zone3d", []),  # Soil 3d if existing
        )

        # loop over dict of perturbated variables, then over the ensemble members
        # ----------------------------------------------------------------------
        for key in list_parm2update:
            key_root = re.split(r"(\d+)", key)
            if len(key_root) == 1:
                key_root.append("0")

            for ens_nb in range(self.NENS):
                with _working_directory(self._member_dir(ens_nb)):
                    self._update_parameter_member(ctx, key, key_root, ens_nb)

    def _member_dir(self, ens_nb: int) -> str:
        """Directory of the CATHY run of ensemble member ``ens_nb`` (0-based)."""
        return os.path.join(
            self.workdir,
            self.project_name,
            "./DA_Ensemble/cathy_" + str(ens_nb + 1),
        )

    def _update_states_from_analysis(
        self, list_parm2update: Sequence[str], cycle_nb: int, analysis: np.ndarray
    ) -> None:
        """Write the analysed pressure heads into the `input/ic` of the valid members.

        Always done: the pressure head is the state variable. Nothing is written
        in the prior cycle (``cycle_nb == 0``).
        """
        for i, ens_nb in enumerate(self.ens_valid):
            with _working_directory(self._member_dir(ens_nb)):
                if "St. var." in list_parm2update:
                    if cycle_nb > 0:
                        self.update_ic(
                            INDP=1,
                            IPOND=0,
                            pressure_head_ini=analysis[:, ens_nb],
                            filename=os.path.join(os.getcwd(), "input/ic"),
                            backup=True,
                            shellprint_update=(
                                i == 0
                            ),  # echo for the first member only
                        )
                elif cycle_nb > 0:
                    raise ValueError(
                        "no state variable update - use last iteration as initial conditions ?"
                    )

    def _update_parameter_member(
        self, ctx: _EnsembleUpdate, key: str, key_root: List[str], ens_nb: int
    ) -> None:
        """Update the input files of one member for one perturbed variable.

        Must run inside the directory of the member. The parameter family is
        chosen from ``key_root[0]`` (name without the zone / layer number).
        """
        name = key_root[0]
        dict_parm_pert, update_key = ctx.dict_parm_pert, ctx.update_key

        # atmbc update
        # if non-homogeneous atmbc, take only the first node perturbated parameters
        # to avoid looping over all the perturbated atmbc nodes as the update is done
        # only once over the stacked perturbated atmbc parameters
        if name.casefold() in "atmbc".casefold():
            self._update_atmbc_member(dict_parm_pert, ens_nb, key)

        # ic update (i.e. water table position update)
        elif name.casefold() in "WTPOSITION".casefold():
            self._update_wtposition_member(
                ctx.cycle_nb, dict_parm_pert, ens_nb, key, True
            )

        # ic update (i.e. initial pressure head update)
        elif name.casefold() in "ic".casefold():
            self._update_ic_member(ctx.cycle_nb, dict_parm_pert, ens_nb, key, True)

        # saturated hydraulic conductivity update
        elif name.casefold() in "ks":
            ctx.ks = self._update_soil_params_from_perturbations(
                dict_parm_pert, ctx.ks, key_root, key, update_key, ens_nb, True
            )

        elif name.casefold() in "porosity":
            ctx.porosity = self._update_soil_params_from_perturbations(
                dict_parm_pert, ctx.porosity, key_root, key, update_key, ens_nb, True
            )
            self.update_archie_parameter(
                _ARCHIE_PARM_NAMES, key_root, update_key, ens_nb, dict_parm_pert, key
            )

        # VG parameters update
        elif name in _VG_PARM_NAMES:
            ctx.vg = self._update_soil_params_from_perturbations(
                dict_parm_pert, ctx.vg, key_root, key, update_key, ens_nb, True
            )

        # Feddes parameters update
        elif name in _FEDDES_PARM_NAMES:
            self._update_feddes_member(ctx, key, key_root, ens_nb)

        # Archie parameters update
        elif name in _ARCHIE_PARM_NAMES:
            self.update_archie_parameter(
                _ARCHIE_PARM_NAMES, key_root, update_key, ens_nb, dict_parm_pert, key
            )

        # the state is updated apart (`_update_states_from_analysis`)
        elif key != "St. var.":
            raise ValueError("key:" + str(key) + " to update is not existing!")

    def _update_feddes_member(
        self, ctx: _EnsembleUpdate, key: str, key_root: List[str], ens_nb: int
    ) -> None:
        """Set one Feddes parameter of one member and rewrite its `input/soil`."""
        df_SPP, _ = in_CT.read_soil(
            os.path.join(
                self.workdir,
                self.project_name,
                f"DA_Ensemble/cathy_{ens_nb+1}",
                "input/soil",
            ),
            self.dem_parameters,
            self.cathyH["MAXVEG"],
        )

        if len(ctx.feddes) == 0:  # one Feddes table per member, built on first use
            ctx.feddes = [
                copy.deepcopy(feddes) for feddes in [self.soil_FP["FP_map"]] * self.NENS
            ]

        new_Feddes_parm = ctx.dict_parm_pert[key][ctx.update_key]
        ctx.feddes[ens_nb][key_root[0]].iloc[int(key_root[1])] = new_Feddes_parm[ens_nb]
        # solution to check to avoid warning: .loc[int(key_root[1])] = new_Feddes_parm[ens_nb]

        self.update_soil(
            FP_map=ctx.feddes[ens_nb],
            SPP_map=df_SPP,
            zone3d=ctx.zone3d,
            verbose=self.verbose,
            filename=os.path.join(os.getcwd(), "input/soil"),
            shellprint_update=True,
            backup=True,
        )

    def _get_data2assimilate(self, list_assimilated_obs, time_ass=None, match=False):
        """
        Loop over observation dictionnary and select corresponding data for a given assimilation time

        Parameters
        ----------
        list_assimilated_obs : list
            list of observation to assimilate.

        Returns
        -------
        data : list
            list of data values.
        """

        # find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        logger.info(f"Assimilated Observations are: {list_assimilated_obs}")

        def extract_data(sensor, time_ass, data):
            data2add = []
            if "ERT" in sensor:
                if "pygimli" in items_dict[time_ass][1][sensor]["data_format"]:
                    data2add.append(
                        np.array(items_dict[time_ass][1][sensor]["data"]["rhoa"])
                    )
                else:
                    data2add.append(
                        items_dict[time_ass][1][sensor]["data"]["resist"].to_numpy()
                    )
                data2add = np.hstack(data2add)
            elif "EM" in sensor:
                data2add.append(
                    np.hstack(
                        items_dict[time_ass][1][sensor]["data"][
                            obs2map[0]["coils"]
                        ].to_numpy()
                    )
                )
            else:
                data2add.append(items_dict[time_ass][1][sensor]["data"])

            return data2add

        if time_ass is None:
            time_ass = (
                self.count_DA_cycle
            )  # + 1 since we compare the predicted observation at ti +1
            logger.warning("!!!!Take care here not sure I am taking the right time!!!")

        data = []
        items_dict = list(self.dict_obs.items())

        # rebuild obskey2map from zero only if match is True (and not in the "all" case)
        rebuild_obskey = match and "all" not in list_assimilated_obs
        if rebuild_obskey:
            obskey2map = []

        for sensor in items_dict[time_ass][1].keys():
            if "all" in list_assimilated_obs:
                data2add = extract_data(sensor, time_ass, data)
                data.append(data2add)
            else:
                if match:
                    if sensor in list_assimilated_obs:
                        data2add = extract_data(sensor, time_ass, data)
                        data.append(data2add)
                        obskey2map.append(sensor)
                else:
                    str_sensor_rootname = re.sub(r"\d", "", sensor)
                    if str_sensor_rootname in list_assimilated_obs:
                        data2add = extract_data(sensor, time_ass, data)
                        data.append(data2add)

        # atleast_1d: a single observation must stay a vector of length 1, not a 0-d array
        return np.atleast_1d(np.squeeze(data)), obskey2map

    def _obs_key_select(self, list_assimilated_obs):
        """find data to map from dictionnary of observations"""
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obskey2map = []  # data type 2 map
        obs2map = []  # data dict 2 map
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[self.count_DA_cycle][1].keys():
            if "all" in list_assimilated_obs:
                obskey2map.append(sensor)
                obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
            else:
                str_sensor_rootname = re.sub(r"\d", "", sensor)
                if str_sensor_rootname in list_assimilated_obs:
                    obskey2map.append(sensor)
                    obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
        return obskey2map, obs2map

    def _check_prediction_shape(self, Hx_ens, obskey2map):
        if isinstance(Hx_ens, list):
            try:
                Hx_ens = np.array(Hx_ens)
            except (
                ValueError
            ):  # ragged list: keep the list, shape check below reports it
                pass
        if not isinstance(Hx_ens, float):
            Hx_ens_shape = np.shape(Hx_ens)
            expected_shape = (len(obskey2map), len(self.ens_valid))
            if Hx_ens_shape != expected_shape:
                if Hx_ens_shape == (len(self.ens_valid), 1, len(obskey2map)):
                    Hx_ens = np.vstack(Hx_ens).T
                elif Hx_ens_shape[0] == expected_shape[1]:
                    Hx_ens = np.array(Hx_ens).T

        return Hx_ens

    def map_states2Observations(
        self, list_assimilated_obs="all", parallel=False, default_state="psi", **kwargs
    ):
        """
        Translate (map) the state values (pressure head or saturation) to observations (or predicted) values


        In other words: apply H mapping operator to convert model to predicted value Hx

        In practice, for each assimilation time, loop over the observation
        dictionnary and infer observation type. Stack Hx.
        Hx = [Hx0,Hx1,Hx_nb_of_observation]

        .. note::

            - For ERT data, H is Archie law, SWC --> ER0 --> ERapp
            - For SWC data, H is a multiplicator (porosity)
            - For tensiometer data, no mapping needed or VGP model if model state used is saturation
            - For discharge, no mapping needed
            - For scale data: no mapping needed


        .. note::

            - For ERT data each ERT mesh node is an observation

        Parameters
        ----------
        state : np.array([])
            [pressure heads, saturation] for each mesh nodes. The default is [None, None].
        list_assimilated_obs : TYPE, optional
            DESCRIPTION. The default is 'all'.
        parallel : Bool, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        Hx_ens : np.array (size: nb of observations x nb of ensembles)
            Ensemble of the simulated (predicted) observations.
        """

        write2shell_map = kwargs.get("write2shell_map", True)

        # observations to map at this DA cycle (the same for every member)
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        self.isET = any("ETact" in item for item in obskey2map)

        Hx_ens = []  # matrice of predicted observation for each ensemble realisation
        for ens_nb in self.ens_valid:  # loop over ensemble files
            Hx_stacked = self._predict_member(ens_nb, obskey2map, obs2map)
            if len(Hx_stacked) > 0:
                Hx_ens.append(np.hstack(Hx_stacked))

        Hx_ens = self._check_prediction_shape(Hx_ens, obskey2map)

        # ERT is mapped for the whole ensemble at once (sequential assimilation). The
        # mapper takes the ERT survey of the cycle (`dict_obs[...]["ERT"]`) whatever the
        # key that triggered it: more than one key would repeat the same forward model
        # and stack the same predictions twice.
        ert_keys = [obs_key for obs_key in obskey2map if "ERT" in obs_key]
        if len(ert_keys) > 1:
            raise ValueError(
                "only one ERT survey per assimilation time is supported, got "
                + str(ert_keys)
            )
        if ert_keys:
            ert_inputs = self._ert_ensemble_inputs() if parallel else None
            Hx_ens = self._map_ERT(Hx_ens, parallel, ert_inputs, write2shell_map)

        return Hx_ens  # meas_size * ens_size

    def _predict_member(
        self, ens_nb: int, obskey2map: Sequence[str], obs2map: Sequence[dict]
    ) -> list:
        """Predicted observations of one ensemble member (list of arrays, one per observation).

        Applies the mapping operator H of each observation type to the last
        state (pressure head, saturation) of the member.
        """
        path_fwd_CATHY = os.path.join(
            self.workdir, self.project_name, f"DA_Ensemble/cathy_{ens_nb + 1}"
        )

        # soil properties of this member
        SPP_ensi, _ = in_CT.read_soil(
            os.path.join(path_fwd_CATHY, self.input_dirname, "soil"),
            self.dem_parameters,
            self.cathyH["MAXVEG"],
        )
        _, POROS_mesh_nodes_ensi = self.map_prop_2mesh_markers(
            "POROS",
            SPP_ensi["POROS"].to_list(),
            to_nodes=False,
        )
        logger.debug("ens nb:" + str(ens_nb))
        logger.debug("mean porosity:" + str(np.mean(POROS_mesh_nodes_ensi)))

        # last state of the member
        df_psi = self.read_outputs(
            filename="psi", path=os.path.join(path_fwd_CATHY, self.output_dirname)
        )
        df_sw, _ = self.read_outputs(
            filename="sw", path=os.path.join(path_fwd_CATHY, self.output_dirname)
        )
        state = [df_psi.iloc[-1], df_sw.iloc[-1].values]

        Hx_stacked = []  # stacked predicted observation
        if self.isET:
            Hx_stacked.append(self._predict_ET(path_fwd_CATHY, obs2map))

        # Loop over observations to map
        # ---------------------------------------------------------------------
        for i, obs_key in enumerate(obskey2map):

            if "tensio" in obs_key:
                obs2map_node = obs2map[i]["mesh_nodes"]
                Hx_stacked.append(mapper.tensio_mapper(state, obs2map_node))

            if "swc" in obs_key:
                # sw assimilation (Hx_SW)
                obs2map_node = obs2map[i]["mesh_nodes"]
                Hx_stacked.append(
                    mapper.swc_mapper(
                        state,
                        obs2map_node,
                        self.grid3d["mesh3d_nodes"],
                        self.dem_parameters,
                        SPP_ensi,
                        self.console,
                    )
                )

            # Atmact-vf(13) / -v (14) / -r (15) / -d (16): infiltration or
            # exfiltration at atmospheric BC nodes; discharge: read `hgsfdet`
            if "stemflow" in obs_key:
                logger.warning("Not yet implemented")
            if "scale" in obs_key:
                logger.warning("Not yet implemented")
            if "discharge" in obs_key:
                logger.warning("Not yet implemented")

            if "EM" in obs_key:
                # Electromagnetic surveys
                Hx_EM, _ = mapper._map_EM(
                    self.dict_obs,
                    self.project_name,
                    self.Archie_parms,
                    self.count_DA_cycle,
                    path_fwd_CATHY,
                    ens_nb,
                    POROS_mesh_nodes_ensi,
                    self.grid3d,
                    verbose=True,
                )
                Hx_stacked.append(Hx_EM)

        return Hx_stacked

    def _ert_ensemble_inputs(self) -> Tuple[List[str], np.ndarray]:
        """Paths and mesh-node porosities of the valid members, for the parallel ERT forward model."""
        ensemble_fwd_paths = []  # list of ensemble paths
        ensemble_porosity_arrays = (
            np.ones(self.NENS, dtype=object) * self.Archie_parms["porosity"]
        )  # preallocate object array
        logger.debug("valid members: %s", self.ens_valid)

        for ens_nb in self.ens_valid:  # loop over ensemble files
            forward_path_i = os.path.join(
                self.workdir, self.project_name, f"DA_Ensemble/cathy_{ens_nb + 1}"
            )
            ensemble_fwd_paths.append(forward_path_i)

            SPP_ensi, _ = in_CT.read_soil(
                os.path.join(forward_path_i, self.input_dirname, "soil"),
                self.dem_parameters,
                self.cathyH["MAXVEG"],
            )
            # porosity of this member, mapped to the mesh nodes
            _, POROS_mesh_nodes_ensi = self.map_prop_2mesh_markers(
                "POROS", SPP_ensi["POROS"].to_list(), to_nodes=False
            )
            ensemble_porosity_arrays[ens_nb] = POROS_mesh_nodes_ensi

        first_porosity = ensemble_porosity_arrays[0]
        unique_count = (
            0 if first_porosity is None else len(np.unique(np.round(first_porosity)))
        )
        logger.debug(
            "Dimension of porosity used in Archie (first ensemble): %s", unique_count
        )
        return ensemble_fwd_paths, ensemble_porosity_arrays

    def _map_ERT(
        self,
        Hx_ens: np.ndarray,
        parallel: bool,
        ert_inputs: Optional[tuple],
        write2shell_map: bool,
    ) -> np.ndarray:
        """Append the ERT predictions of the ensemble to ``Hx_ens`` (rows = observations)."""
        if write2shell_map:
            self.console.print(":zap: Mapping to ERT prediction")
            self.console.print(
                r"""
                                    Fwd modelling ERT
                                    Noise level: {??}
                                    """,
                style="green",
            )
            self.console.rule(
                ":octopus: Parameter perturbation :octopus:", style="green"
            )
            self.console.print(
                r"""
                                    Archie perturbation for DA analysis
                                    Archie rFluid: {}
                                    Pert. Sigma Archie: {}
                                    Nb of zones (?): {}
                                    """.format(
                    self.Archie_parms["rFluid_Archie"],
                    self.Archie_parms["pert_sigma_Archie"],
                    len(self.Archie_parms["rFluid_Archie"]),
                ),
                style="green",
            )
            self.console.rule("", style="green")

        if not parallel:
            raise NotImplementedError(
                "ERT observations are only mapped with parallel=True"
            )

        ensemble_fwd_paths, ensemble_porosity_arrays = ert_inputs
        Hx_ens_ERT, df_Archie, _ = mapper._map_ERT_parallel(
            self.dict_obs,
            self.count_DA_cycle,
            self.project_name,
            self.Archie_parms,
            ensemble_porosity_arrays,
            ensemble_fwd_paths,
            list_assimilated_obs="all",
            default_state="psi",
            DA_cnb=self.count_DA_cycle,
            savefig=True,
            verbose=True,
        )
        self.df_Archie = mapper._add_2_ensemble_Archie(self.df_Archie, df_Archie)
        if len(Hx_ens) > 0:
            return np.vstack([Hx_ens, Hx_ens_ERT])
        return Hx_ens_ERT

    def _predict_ET(self, path_fwd_CATHY, obs2map):
        """Predicted actual evapotranspiration of one ensemble member at the observed nodes.

        Reads ``fort.777`` of the member and returns the ``ACT. ETRA`` column of the
        last time step at the surface nodes of the ``ETact`` observations.
        """
        df_fort777 = out_CT.read_fort777(os.path.join(path_fwd_CATHY, "fort.777"))
        df_fort777 = df_fort777.set_index("time_sec")
        last_time = df_fort777.index.unique()[-1]

        nodes = [obs["mesh_nodes"] for obs in obs2map if obs["data_type"] == "ETact"]
        # flat integer positions (a list of arrays is an error with recent pandas)
        selected = df_fort777.loc[last_time]["ACT. ETRA"].iloc[np.hstack(nodes)]
        return selected.values

    def _read_state_ensemble(self):
        # read grid to infer dimension of the ensemble X
        # --------------------------------------------------------------------
        if len(self.grid3d) == 0:
            self.grid3d = out_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name, "output", "grid3d")
            )

        M_rows = np.shape(self.grid3d["mesh3d_nodes"])[0]
        N_col = self.NENS

        # Ensemble matrix X of M rows and N  + fill with psi values
        # --------------------------------------------------------------------
        psi = np.full([M_rows, N_col], 1e99)
        sw = np.full([M_rows, N_col], 1e99)
        for j in range(self.NENS):
            # for j in range(len(ens_valid)):
            try:
                df_psi = out_CT.read_psi(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(j + 1),
                        "output/psi",
                    )
                )
                psi[:, j] = df_psi.iloc[-1, :]
            except _READ_ERRORS:
                pass  # member left at its initial value; see the consistency check below

            try:
                df_sw, _ = out_CT.read_sw(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(j + 1),
                        "output/sw",
                    )
                )
                sw[:, j] = df_sw.iloc[-1, :]
            except _READ_ERRORS:
                pass  # member left at its initial value; see the consistency check below
        # members still at the fill value could not be read (rejected members are ignored)
        if np.count_nonzero(psi[:, self.ens_valid] == 1e99) != 0:
            logger.warning("!!!unconsistent filled X, missing values!!!")

        return psi, sw, N_col, M_rows

    def update_pert_parm_dict(self, update_key, list_update_parm, analysis_param_valid):
        """update dict of perturbated parameters i.e. add a new entry with new params"""

        index_update = 0
        if self.count_DA_cycle > 0:
            for pp in list_update_parm:
                if "St. var." in pp:
                    index_update = index_update + 1
                    pass
                else:
                    update_key = _update_key(self.count_DA_cycle)
                    self.dict_parm_pert[pp][update_key] = analysis_param_valid[
                        index_update - 1
                    ]
                    index_update = index_update + 1

            # check after analysis
            # ----------------------------------------------------------------
            id_valid_after = self._check_after_analysis(update_key, list_update_parm)
            intersection_set = set.intersection(
                set(id_valid_after), set(self.ens_valid)
            )
            self.ens_valid = list(intersection_set)

        pass

    def _update_input_ensemble(self, list_parm2update=None, analysis=None, **kwargs):
        """
        Select the time window of the hietograph.
        Write new variable perturbated/updated value into corresponding input file

        Parameters
        ----------
        parm_pert : dict
            dict of all perturbated variables holding values and metadata.

        Returns
        -------
        New file written/overwritten into the input dir.

        """
        if list_parm2update is None:
            list_parm2update = ["St. var."]
        if analysis is None:
            analysis = []

        self.console.print(":sponge: [b]update input ensemble[/b]")

        # resample atmbc file for the given DA window
        # ---------------------------------------------------------------------
        self.selec_atmbc_window()
        # ---------------------------------------------------------------------
        self.apply_ensemble_updates(
            self.dict_parm_pert,
            list_parm2update,
            cycle_nb=self.count_DA_cycle,
            analysis=analysis,
        )

        pass

    def selec_atmbc_window(self):
        if len(self.grid3d) == 0:
            self.run_processor(IPRT1=3, DAFLAG=0, verbose=True)
            self.grid3d = out_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name, "output", "grid3d")
            )

        df_atmbc_full, HSPATM, IETO = in_CT.read_atmbc(
            os.path.join(self.workdir, self.project_name, "input", "atmbc"),
            grid=self.grid3d,
        )

        if self.count_atmbc_cycle is None:
            raise ValueError("count_atmbc_cycle must be set to select an atmbc window")

        for ens_nb in range(self.NENS):
            cwd = os.path.join(
                self.workdir, self.project_name, "./DA_Ensemble/cathy_" + str(ens_nb + 1)
            )

            df_atmbc = df_atmbc_full  # reset for each member

            # perturbed atmbc (homogeneous case only)
            if "atmbc0" in self.dict_parm_pert.keys():
                if HSPATM != 0:
                    VALUE = np.array(
                        self.dict_parm_pert["atmbc0"]["time_variable_perturbation"]
                    )[ens_nb]
                    times = self.dict_parm_pert["atmbc0"]["data2perturbate"]["time"]
                    df_atmbc = pd.DataFrame(np.c_[times, VALUE], columns=["time", "value"])
                else:
                    print("atmbc perturbation not supported for spatially variable atmbc: ignored")

            times_all = df_atmbc.time.unique()
            i = self.count_atmbc_cycle
            if i + 1 >= len(times_all):
                raise IndexError(
                    f"count_atmbc_cycle={i} is beyond the last atmbc interval "
                    f"({len(times_all) - 1} intervals available)"
                )
            time_window_atmbc = [int(times_all[i]), int(times_all[i + 1])]

            df_atmbc_window = df_atmbc[df_atmbc["time"].isin(time_window_atmbc)]
            diff_time = time_window_atmbc[1] - time_window_atmbc[0]

            self.update_parm(
                TIMPRTi=[0, diff_time],
                TMAX=diff_time,
                filename=os.path.join(cwd, "input/parm"),
                IPRT1=2,
                backup=True,
            )

            VALUE = []
            for t in time_window_atmbc:
                values = df_atmbc_window.loc[df_atmbc_window["time"] == t, "value"].values
                VALUE.append(values.item() if values.size == 1 else values)

            self.update_atmbc(
                HSPATM=HSPATM,
                IETO=IETO,
                time=[0, diff_time],
                netValue=VALUE,
                filename=os.path.join(cwd, "input/atmbc"),
            )

        pass

    def _create_empty_ET_xr(self):
        x = np.unique(self.grid3d["mesh3d_nodes"][:, 0])
        y = np.unique(self.grid3d["mesh3d_nodes"][:, 1])
        nan_dataset = xr.Dataset(
            {
                "time_sec": self.all_atmbc_times[self.count_atmbc_cycle],
                "SURFACE NODE": (("x", "y"), np.full((len(x), len(y)), np.nan)),
                "ACT. ETRA": (("x", "y"), np.full((len(x), len(y)), np.nan)),
            },
            coords={
                "x": x,
                "y": y,
            },
        )
        return nan_dataset

    def _DA_ET_xr(self):

        ET_ftime_ensi = []
        for nensi in range(self.NENS):
            if nensi in self.ens_valid:
                try:
                    df_fort777 = out_CT.read_fort777(
                        os.path.join(
                            self.workdir,
                            self.project_name,
                            f"DA_Ensemble/cathy_{nensi+1}",
                            "fort.777",
                        )
                    )
                    df_fort777 = df_fort777.drop_duplicates()
                    df_fort777.columns = [
                        col.lower() if col.lower() in ["x", "y"] else col
                        for col in df_fort777.columns
                    ]
                    df_fort777 = df_fort777.set_index(["time", "x", "y"]).to_xarray()
                    df_fort777_selec_ti = df_fort777.isel(time=-1)
                    df_fort777_selec_ti = df_fort777_selec_ti.rio.set_spatial_dims(
                        "x", "y"
                    )
                    df_fort777_selec_ti = df_fort777_selec_ti.reset_coords(
                        "time", drop=True
                    )
                    df_fort777_selec_ti["time_sec"] = self.all_atmbc_times[
                        self.count_atmbc_cycle
                    ]
                    ET_ftime_ensi.append(df_fort777_selec_ti)

                except _READ_ERRORS as exc:
                    logger.warning(
                        f"Impossible to read ET - exclude scenario ensemble {nensi} ({type(exc).__name__}: {exc})"
                    )
                    nan_dataset = self._create_empty_ET_xr()
                    ET_ftime_ensi.append(nan_dataset)
            else:
                nan_dataset = self._create_empty_ET_xr()
                ET_ftime_ensi.append(nan_dataset)
        ET_DA_xr_time_i = xr.concat(ET_ftime_ensi, dim="ensemble")

        # Concatenate along the "assimilation" dimension if `self.ET_DA_xr` exists, or create it otherwise
        if self.ET_DA_xr is not None:
            self.ET_DA_xr = xr.concat(
                [
                    self.ET_DA_xr,
                    ET_DA_xr_time_i.expand_dims(assimilation=[self.count_atmbc_cycle]),
                ],
                dim="assimilation",
            )
        else:
            self.ET_DA_xr = ET_DA_xr_time_i.expand_dims(
                assimilation=[self.count_atmbc_cycle]
            )

        pass

    def _DA_df(self, state=None, state_analysis=None, rejected_ens=None, **kwargs):
        """
        Build a dataframe container to keep track of the changes before and after update for each ensemble.
        This dataframe is the support for DA results visualisation.


        Parameters
        ----------
        sw : np.array([])
            State variable mesh nodes values (before update).
        sw_analysis : np.array([])
            State variable mesh nodes values after DA analysis.

        Returns
        -------
        None.

        """
        if state is None:
            state = [None, None]
        if rejected_ens is None:
            rejected_ens = []

        df_DA_ti_ni = pd.DataFrame()  # nested df for a given ensemble/given time
        df_DA_ti = pd.DataFrame()  # nested df for a given time

        cols_root = [
            "time",
            "Ensemble_nb",
            "psi_bef_update",
            "sw_bef_update_",
            "analysis",
            "OL",  # Open loop boolean
            "rejected",
        ]

        # Open loop kwargs
        # --------------------------
        t_ass = []
        if "t_ass" in kwargs:
            t_ass = kwargs["t_ass"]

        if "ens_nb" in kwargs:
            ens_nb = kwargs["ens_nb"]
        # --------------------------

        if "openLoop" in kwargs:

            data_df_root = [
                t_ass * np.ones(len(state[0])),
                ens_nb * np.ones(len(state[0])),
                state[0],
                state[1],
                state[0],
                True * np.ones(len(state[0])),
                False * np.ones(len(state[0])),
            ]
            df_DA_ti = pd.DataFrame(np.transpose(data_df_root), columns=cols_root)

        else:

            for n in range(self.NENS):

                data_df_root = [
                    self.count_atmbc_cycle * np.ones(len(state[0])),
                    n * np.ones(len(state[0])),
                    state[0][:, n],
                    state[1][:, n],
                    state_analysis[:, n],
                    False * np.ones(len(state[0])),
                    int(rejected_ens[n]) ** np.ones(len(state[0])),
                ]
                df_DA_ti_ni = pd.DataFrame(
                    np.transpose(data_df_root), columns=cols_root
                )

                df_DA_ti = pd.concat([df_DA_ti, df_DA_ti_ni], axis=0)

        self.df_DA = pd.concat([self.df_DA, df_DA_ti], axis=0, ignore_index=True)

        pass

    def _add_to_perturbated_dict(self, var_per_2add):
        """
        Update dict of perturbated variable.
        """
        self.var_per_dict = self.var_per_dict | var_per_2add
        return self.var_per_dict

    def _create_subfolders_ensemble(self):
        """
        Create ensemble subfolders
        """

        if not os.path.exists(
            os.path.join(self.workdir, self.project_name, "DA_Ensemble/cathy_origin")
        ):
            os.makedirs(
                os.path.join(
                    self.workdir, self.project_name, "DA_Ensemble/cathy_origin"
                )
            )

            # copy input, output and vtk dir
            for dir2copy in enumerate(["input", "output", "vtk"]):
                shutil.copytree(
                    os.path.join(self.workdir, self.project_name, dir2copy[1]),
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_origin",
                        dir2copy[1],
                    ),
                )

        try:
            # for ff in [self.processor_name,'cathy.fnames','prepro']
            # copy exe into cathy_origin folder
            shutil.copy(
                os.path.join(self.workdir, self.project_name, self.processor_name),
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_origin",
                    self.processor_name,
                ),
            )
            # copy cathy.fnames into cathy_origin folder
            shutil.copy(
                os.path.join(self.workdir, self.project_name, "cathy.fnames"),
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_origin/cathy.fnames",
                ),
            )
            # copy prepro into cathy_origin folder
            shutil.copytree(
                os.path.join(self.workdir, self.project_name, "prepro"),
                os.path.join(
                    self.workdir, self.project_name, "DA_Ensemble/cathy_origin/prepro"
                ),
            )

        except OSError as exc:  # FileNotFoundError, FileExistsError, shutil.Error
            self.console.print(
                f":worried_face: [b]could not copy processor files into cathy_origin[/b] ({type(exc).__name__}: {exc})"
            )

        path_origin = os.path.join(
            self.workdir, self.project_name, "DA_Ensemble/cathy_origin"
        )

        # copy origin folder to each ensemble subfolders
        for i in range(self.NENS):
            path_nudn_i = os.path.join(
                self.workdir, self.project_name, "DA_Ensemble/cathy_" + str(i + 1)
            )

            if not os.path.exists(path_nudn_i):
                shutil.copytree(path_origin, path_nudn_i)

    def set_Archie_parm(
        self,
        porosity=None,
        rFluid_Archie=None,
        a_Archie=None,
        m_Archie=None,
        n_Archie=None,
        pert_sigma_Archie=None,
    ):
        """
        Dict of Archie parameters. Each type of soil is describe within a list
        Note that if pert_sigma is not None a normal noise is
        added during the translation of Saturation Water to ER

        Parameters
        ----------
        rFluid : TYPE, list
            Resistivity of the pore fluid. The default is [1.0].
        a : TYPE, list
            Tortuosity factor. The default is [1.0].
        m : TYPE, list
            Cementation exponent. The default is [2.0].
            (usually in the range 1.3 -- 2.5 for sandstones)
        n : TYPE, list
            Saturation exponent. The default is [2.0].
        pert_sigma_Archie : TYPE, list
            Gaussian noise to add. The default is None.

        ..note:
            Field procedure to obtain tce he covariance structure of the model
            estimates is described in Tso et al () - 10.1029/2019WR024964
            "Fit a straight line for log 10 (S) and log 10 (ρ S ) using the least-squares criterion.
            The fitting routine returns the covariance structure of the model estimates, which can be used to de-
            termine the 68% confidence interval (1 standard deviation) of the model estimates.""

        """
        if porosity is None:
            porosity = []
        if rFluid_Archie is None:
            rFluid_Archie = [1.0]
        if a_Archie is None:
            a_Archie = [1.0]
        if m_Archie is None:
            m_Archie = [2.0]
        if n_Archie is None:
            n_Archie = [2.0]
        if pert_sigma_Archie is None:
            pert_sigma_Archie = [0]
        if len(porosity) == 0:
            porosity = self.soil_SPP["SPP"]["POROS"].mean()

        if not isinstance(porosity, list):
            porosity = [porosity]
        self.Archie_parms = {
            "porosity": porosity,
            "rFluid_Archie": rFluid_Archie,
            "a_Archie": a_Archie,
            "m_Archie": m_Archie,
            "n_Archie": n_Archie,
            "pert_sigma_Archie": pert_sigma_Archie,
        }
        pass

    def _performance_assessment_pq(
        self, list_assimilated_obs, data, prediction, t_obs, **kwargs
    ):
        """Calculates (N)RMSE for each observation in a data assimilation cycle."""

        ol_bool = kwargs.get("openLoop", False)
        file_path = os.path.join(
            self.workdir,
            self.project_name,
            self.project_name + "_df_performance.parquet",
        )
        # Rows already on disk (empty frame if the file does not exist yet).
        # Only the rows computed in this call are written back below, otherwise
        # the existing rows would be duplicated at every cycle.
        df_existing = _load_performance_from_parquet(file_path)
        df_new = pd.DataFrame()
        # History used for the cumulative (N)RMSE: previous cycles of the same
        # mode only (open loop and DA runs share the parquet file).
        df_past = df_existing
        if "OL" in df_past:
            df_past = df_past[df_past["OL"] == ol_bool]
        if "time" in df_past:
            df_past = df_past[df_past["time"] < t_obs]
        # RMSE_avg is repeated on every sensor row of a cycle: keep one value per cycle
        if "RMSE_avg" in df_past:
            past_rmse_avg_sum = df_past.drop_duplicates("time")["RMSE_avg"].sum()
        else:
            past_rmse_avg_sum = 0.0

        all_obs_diff_mat = _compute_diff_matrix(data, prediction)
        all_obs_rmse_avg = _compute_rmse(all_obs_diff_mat)
        nrmse_avg = (past_rmse_avg_sum + all_obs_rmse_avg) / (t_obs + 1)

        # Example function to retrieve list of keys for assimilated observations
        obs2eval_key = self._get_data2assimilate(list_assimilated_obs)[1]
        start_line_obs = 0

        for name_sensor in obs2eval_key:
            obs2eval, _ = self._get_data2assimilate([name_sensor], match=True)
            # rows of this observation in the prediction: its number of DATA (one per
            # sensor for point sensors, but e.g. 91 for one ERT survey under one key)
            n_obs = len(obs2eval)
            prediction2eval = prediction[start_line_obs : start_line_obs + n_obs]

            rmse_sensor = _get_rmse_for_sensor(obs2eval, prediction2eval)

            col = "RMSE" + name_sensor
            if "ObsType" in df_past and col in df_past:
                past_rmse_sum = df_past.loc[
                    df_past["ObsType"] == name_sensor, col
                ].sum()
            else:
                past_rmse_sum = 0.0
            nrmse_sensor = (past_rmse_sum + rmse_sensor) / (t_obs + 1)

            df_new = _append_to_performance_df(
                df_new,
                t_obs,
                name_sensor,
                rmse_sensor,
                all_obs_rmse_avg,
                nrmse_sensor,
                nrmse_avg,
                ol_bool,
            )

            start_line_obs += n_obs

        # Save to Parquet and clear memory if too large
        # if len(self.df_performance) > 10000:  # Adjust batch size as needed
        _save_performance_to_parquet(df_new, file_path)  # appends to the file
        self.df_performance = (
            pd.DataFrame()
        )  # kept for backward compatibility (always empty)
