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
from typing import Any, List, Optional, Sequence, Tuple, Union
import glob

import multiprocessing

import logging
import json
import os
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
from pyCATHY.DA import mapper, localisation
from pyCATHY.ERT import ERT_checks
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as plt_CT
from pyCATHY.DA.analysis import (  # noqa: F401  (re-exported: the code moved to its own module)
    ANALYSES,
    _ENKF_DIAGNOSTICS,
    _LEGACY_DA_KEYS,
    AnalysisConfig,
    AnalysisResult,
    DAConfig,
    _enkf,
    _enkf_localized_with_inflation,
    _from_enkf,
    _particle_filter,
    run_analysis,
)
from pyCATHY.DA.param_tables import (  # noqa: F401
    _ARCHIE_PARM_NAMES,
    _FEDDES_PARM_NAMES,
    _SOIL_PARM_COLUMNS,
    _VG_PARM_NAMES,
    _assign_soil_value,
    _check_positive_pressure_heads,
    _has_name_token,
    _nonphysical_values,
    _set_feddes_value,
    _update_key,
    _valid_zone_values,
)
from pyCATHY.DA.performance_tables import (  # noqa: F401
    _append_to_performance_df,
    _compute_diff_matrix,
    _compute_rmse,
    _get_rmse_for_sensor,
    _load_performance_from_parquet,
    _save_performance_to_parquet,
)
from pyCATHY.DA.runtime import (  # noqa: F401
    _LOG_DATEFMT,
    _LOG_FORMAT,
    _Tee,
    _configure_logging,
    _restores_cwd,
    _tee_run_output,
    _working_directory,
)

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


_LOG_EPS = 1e-3  # offset added before log10 (same as in _log_transform_obs)


def _log10_cov(data_valid: np.ndarray, data_cov_valid: np.ndarray, eps: float = _LOG_EPS):
    """Data covariance in the log10 domain (first-order propagation).

    R_log = J R J with J = diag(1 / (ln(10) (d + eps))). Non-positive d + eps
    (which cannot be log-transformed) get J = 0 and a warning.
    """
    # non-positive observations are clipped to 0 first (as in _log_transform_obs),
    # so that the Jacobian is defined: log10(max(d, 0) + eps)
    d = np.clip(np.asarray(data_valid, dtype=float), 0.0, None) + eps
    jac = 1.0 / (np.log(10.0) * d)
    return np.asarray(data_cov_valid, dtype=float) * np.outer(jac, jac)






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





def _warn_once(key: str, message: str) -> None:
    """Emit a warning only the first time `key` is seen (avoids one warning per analysis step)."""
    if key not in _WARNED_ONCE:
        _WARNED_ONCE.add(key)
        warnings.warn(message)






def _log_once(key: str, message: str, level: int = logging.WARNING) -> None:
    """Log a message only the first time `key` is seen (once per process, not once per cycle)."""
    if key not in _WARNED_ONCE:
        _WARNED_ONCE.add(key)
        logger.log(level, message)




# ---------------------------------------------------------------------------
# ANALYSIS: one adapter per DA_type, selected through the ANALYSES registry
# ---------------------------------------------------------------------------
# An adapter has the signature
#     adapter(data, data_cov, param, state, prediction, config) -> AnalysisResult
# where `state` is the ensemble matrix to update (psi or sw) and `config` is an
# `AnalysisConfig` (Sakov, localisation matrix, inflation, ...).
# The adapters hide the differences between the filters' return values.
























# ---------------------------------------------------------------------------
# SMALL HELPERS (pure functions, no `self`)
# ---------------------------------------------------------------------------










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
        self.ert_check_reports: List[dict] = []  # ERT checks (see ERT_checks), saved every cycle
        self.rejection_log: List[dict] = []  # who was rejected, when and why (saved as json)
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
        self.ert_check_reports = []
        self.rejection_log = []
        mapper.pop_ert_reports()  # drop reports of a previous run

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
        self.log_transform_cov = config.log_transform_cov
        self.mask_bad_ert = config.mask_bad_ert
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

        self._log_perturbation_layout(list_parm2update)
        self._log_time_alignment()

        # update the ensemble files with pertubarted parameters
        # -------------------------------------------------------------------
        self.apply_ensemble_updates(
            dict_parm_pert,
            list_parm2update="all",
            cycle_nb=self.count_DA_cycle,
        )

        all_atmbc_times = np.copy(self.atmbc["time"])

        return (callexe, all_atmbc_times, threshold_rejected, ENS_times, update_key)

    def describe_perturbation(self, list_parm2update=()):
        """What every perturbed parameter touches (soil rows, Feddes row), as a DataFrame.

        Uses ``self.dict_parm_pert`` (set by ``run_DA_sequential`` or by the user) and the
        soil table, zone map and vegetation map of the project. ``df.attrs["warnings"]``
        lists writes that override each other. See ``pyCATHY.DA.scenario.describe_ensemble``.
        """
        from pyCATHY.DA.scenario import describe_ensemble

        entries = [
            {"type_parm": key, **val}
            for key, val in (self.dict_parm_pert or {}).items()
            if isinstance(val, dict)
        ]
        soil_index = zone = veg = None
        try:
            soil_index = self.read_inputs("soil")[0].index
        except Exception as exc:  # no soil table yet: the rows are then not described
            logger.debug("describe_perturbation: no soil table (%s)", exc)
        for name, attr in (("zone", "zone"), ("root_map", "veg")):
            try:
                value = self.read_inputs(name)[0]
            except Exception as exc:
                logger.debug("describe_perturbation: no %s map (%s)", name, exc)
                continue
            zone, veg = (value, veg) if attr == "zone" else (zone, value)
        return describe_ensemble(entries, list(list_parm2update), soil_index, zone, veg)

    def _log_perturbation_layout(self, list_parm2update):
        """Log once, at the start of a run, which rows each perturbed parameter writes.

        Diagnostics only: whatever goes wrong here must not stop the run.
        """
        try:
            table = self.describe_perturbation(list_parm2update)
        except Exception as exc:
            logger.debug("perturbation layout not described (%s)", exc)
            return
        if table.empty:
            return
        cols = ["key", "control", "target", "mean", "sd", "min", "max", "updated"]
        logger.info("Perturbed parameters:\n%s", table[cols].to_string(index=False))
        for message in table.attrs.get("warnings", []):
            logger.warning("perturbation layout: %s", message)

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
        logger.info("Analysis: DA type=%s, inflation=%s", DA_type, self.inflation)

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
    def _run_DA_sequential_impl(
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

    def run_DA_sequential(self, *args, **kwargs):
        """Run the sequential Data Assimilation (see ``_run_DA_sequential_impl``).

        Everything printed during the run is also written, live, to
        ``<workdir>/<project_name>/console_log.txt``.
        """
        log_path = os.path.join(self.workdir, self.project_name, "console_log.txt")
        with _tee_run_output(log_path, getattr(self, "console", None)):
            return self._run_DA_sequential_impl(*args, **kwargs)

    run_DA_sequential.__doc__ = (run_DA_sequential.__doc__ or "") + (
        "\n" + (_run_DA_sequential_impl.__doc__ or "")
    )

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
        t_end = self.all_atmbc_times[self.count_atmbc_cycle + 1]
        if t_end in self.ENS_times:          # was: if t_atmbc in self.ENS_times   
            state_bef_update, analysis_valid = self._analysis_cycle(
                options.DA_type,
                options.list_parm2update,
                options.list_assimilated_obs,
                options.update_key,
                options.parallel,
                t_window_start=t_atmbc,
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
            logger.info("End of DA")

        self._log_end_of_da_update()
        plt.close("all")

    def _analysis_cycle(
        self,
        DA_type: str,
        list_parm2update: Sequence[str],
        list_assimilated_obs: Union[str, Sequence[str]],
        update_key: str,
        parallel: bool,
        t_window_start: Optional[float] = None,
    ) -> Tuple[list, np.ndarray]:
        """One DA cycle at an observation time: map, analyse, check, record.

        Returns
        -------
        (list, array)
            ``[psi, sw]`` of the ensemble before the update, and the valid
            analysed ensemble (states, plus parameters when updated).
        """
        self._log_assimilation_times(t_window_start)
        logger.info("Map states to observations")
        # H operator applied to the state variable
        prediction = self.map_states2Observations(
            list_assimilated_obs,
            default_state="psi",
            parallel=parallel,
        )
        self._validate_prediction(prediction, list_assimilated_obs)

        logger.info("Analysis: DA type=%s, inflation=%s", DA_type, self.inflation)
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

        logger.info("Check analysis performance")
        self._performance_assessment_pq(
            list_assimilated_obs,
            data,
            prediction_valid,
            t_obs=self.count_DA_cycle,
        )
        self.count_DA_cycle = self.count_DA_cycle + 1

        _check_positive_pressure_heads(ensemble_psi_valid, "forecast")
        _check_positive_pressure_heads(analysis_valid, "analysed")  # effect of the update (B10)
        self.update_pert_parm_dict(update_key, list_parm2update, analysis_param_valid)
        return [psi_bef_update, sw_bef_update], analysis_valid

    def _log_time_alignment(self):
        """Log, before the first window, how the observation times line up with the model windows.

        Diagnostics only (see ``time_alignment``): an observation that cannot be assimilated, or
        the one-window offset between an observation and the state it is compared with, is
        reported here instead of showing up as a surprising result.
        """
        from pyCATHY.DA.time_alignment import check_time_alignment

        try:
            findings = check_time_alignment(self.ENS_times, self.all_atmbc_times)
        except Exception as exc:  # never stop a run for a diagnostic
            logger.debug("time alignment not checked (%s)", exc)
            return
        for finding in findings:
            log = logger.error if finding.level == "error" else logger.warning
            log("observation times: %s", finding.message)

    def _log_assimilation_times(self, t_window_start):
        """Log the times involved in one analysis, to check the time alignment (B4).

        The model has just been run over ``[window start, window end]`` (consecutive
        atmbc times). The observation of this cycle is stamped ``t_obs`` (``ENS_times``).
        The three numbers are logged so that the observation time can be compared with
        the time of the states it is confronted with.
        """
        t_obs = self.ENS_times[self.count_DA_cycle] if self.count_DA_cycle < len(self.ENS_times) else None
        i = self.count_atmbc_cycle
        t_end = self.all_atmbc_times[i + 1] if i + 1 < len(self.all_atmbc_times) else None
        logger.info(
            "Analysis %d: observation stamped t=%s s; model window [%s, %s] s",
            self.count_DA_cycle + 1, t_obs, t_window_start, t_end,
        )
        if t_obs is not None and t_end is not None and t_obs != t_end:
            logger.info(
                "observation time (%s) differs from the end of the model window (%s): "
                "check that this is the intended alignment",
                t_obs, t_end,
            )

    def _no_analysis_cycle(self) -> Tuple[list, np.ndarray]:
        """Time step without observation: only the hydrological model was run.

        Returns the same pair as `_analysis_cycle`; the ensemble is not
        validated and the "analysis" is the unchanged pressure-head ensemble.
        """
        logger.info("No observation at this time: running the hydrological model only")
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
        logger.info("End of time step %d/%d s", int(t_atmbc), int(self.all_atmbc_times[-1]))
        logger.info("End of atmbc update %d/%d", self.count_atmbc_cycle, len(self.all_atmbc_times))

    def _log_end_of_da_update(self):
        """Console progress after a DA update."""
        logger.info("End of DA update %d/%d", self.count_DA_cycle, len(self.ENS_times))
        logger.info("Valid ensemble members: %d/%d (%.0f%%)", len(self.ens_valid), self.NENS, 100 * len(self.ens_valid) / self.NENS)

    def _save_step_outputs(
        self,
        state,
        state_analysis,
        rejected_ens,
        list_assimilated_obs,
        list_parm2update,
    ):
        """Store the results of the current step (dataframes, ET, backups)."""
        logger.debug("Back up DA step")
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
        self._save_ert_check_reports()
        if self.config.save_results:
            self._save_results_to_disk()

    def _save_ert_check_reports(self):
        """Write the ERT check reports next to the project and log how many failed."""
        if not self.ert_check_reports:
            return
        path = os.path.join(
            self.workdir, self.project_name, self.project_name + "_ERT_checks.json"
        )
        ERT_checks.save_reports(self.ert_check_reports, path)
        n_failed = ERT_checks.count_failed(self.ert_check_reports)
        if n_failed:
            logger.warning(f"{n_failed} ERT check report(s) failed - see {path}")

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
        logger.info("Domain localisation: %s", self.localisation)

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
        # loci indexes the zones WITHOUT the nodata value: it is the position in
        # mesh3d_nodes_valid (built without -9999) and the zone's parameter in
        # list_parm2update. Enumerating np.unique() including -9999 shifted it by one.
        for loci, node_value in enumerate(_valid_zone_values(mask_surface_nodes)):

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
            _log_once("loc_shortcut", "Domain localisation: all the parameters are updated, whether or not they are associated with the observation localisation")
            logger.debug("Params to update (domain loc): %s", selected_parm2update)

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
        and its own parameter, paired **by name** (``ZROOT0`` belongs to zone 0; see
        ``scenario.select_zone_parameters``, which keeps the old positional rule for names
        without a numeric suffix).
        """
        if self.localisation in ("veg_map", "zones"):
            from pyCATHY.DA.scenario import select_zone_parameters

            return select_zone_parameters(list_parm2update, loci)
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

        logger.warning("Not all the mesh nodes were analysed with DA: taking the state values where the analysis is missing")
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
        logger.info("Assimilating %s: %d data", list_assimilated_obs, np.size(data))
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

        # ERT: drop the unusable quadrupoles from the filter inputs only (the full-length
        # data are returned, the performance assessment compares them to the full prediction)
        data_full = data_valid
        keep = self._ert_quality_mask(obskey2map, len(data_valid), obs_valid)
        if keep is not None and not keep.all():
            data_valid = data_valid[keep]
            data_cov_valid = data_cov_valid[keep, :][:, keep]
            prediction_valid = prediction_valid[keep, :]
            analysis_config = self._mask_localisation_matrix(analysis_config, keep)

        # ERT: data covariance, predictions and innovation as they enter the filter
        if any("ERT" in key for key in obskey2map):
            # the check transforms data/predictions itself: give it the R the filter will use
            cov_checked = (
                _log10_cov(data_valid, data_cov_valid)
                if self.use_log_transformed_obs and getattr(self, "log_transform_cov", True)
                else data_cov_valid
            )
            ert_report = ERT_checks.check_assimilation_inputs(
                data_valid,
                prediction_valid,
                cov_checked,
                use_log=self.use_log_transformed_obs,
                ens_ids=[int(i) + 1 for i in ens_valid],
                name=f"ERT analysis cycle {self.count_DA_cycle}",
            ).warn()
            self.ert_check_reports.append(
                {
                    "stage": "analysis",
                    "cycle": self.count_DA_cycle,
                    "report": ert_report.to_dict(),
                }
            )

        # observations and predictions as seen by the filter
        if self.use_log_transformed_obs:
            data_in, prediction_in = self._log_transform_obs(
                data_valid, prediction_valid, data_cov_valid
            )
            if getattr(self, "log_transform_cov", True):
                data_cov_valid = _log10_cov(data_valid, data_cov_valid)
                logger.debug(
                    "R propagated to log10: diag %.2e-%.2e",
                    np.diag(data_cov_valid).min(),
                    np.diag(data_cov_valid).max(),
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
        # (the result used to be discarded: parameters with a transf_type were
        # stored, and then transformed again at the next cycle, in the
        # transformed space)
        analysis_param = self._back_transform_analysis_param(
            list_update_parm, analysis_param
        )
        return (ensemble_psi, ensemble_sw, data_full, analysis, analysis_param)

    def _mask_localisation_matrix(self, analysis_config, keep):
        """Drop from the covariance-localisation matrix the columns of the excluded quadrupoles.

        The matrix has one column per observation. It is given for the observations before the
        quality mask (``keep.size``); the filter sees only the kept ones, and a matrix with the
        old width fails there with a shape error.
        """
        matrix = getattr(analysis_config, "localisation_matrix", None)
        if matrix is None:
            return analysis_config
        matrix = np.asarray(matrix)
        if matrix.ndim == 2 and matrix.shape[1] == keep.size:
            return replace(analysis_config, localisation_matrix=matrix[:, keep])
        if matrix.ndim == 2 and matrix.shape[1] == int(keep.sum()):
            return analysis_config  # already built for the kept observations
        logger.warning(
            "localisation matrix has %s columns but %d observations were available (%d kept after "
            "the ERT quality mask): it is left unchanged and will probably be rejected",
            matrix.shape[1] if matrix.ndim == 2 else matrix.shape, keep.size, int(keep.sum()),
        )
        return analysis_config

    def _ert_quality_mask(self, obskey2map, n_data, obs_valid=None):
        """Mask (True = keep) of the usable ERT data among the ``n_data`` filter inputs.

        Non-ERT data are always kept. Returns None when masking is off, when there
        is no ERT data, or when the layout of the stacked data cannot be matched
        (then nothing is dropped and a warning is logged).
        """
        if not getattr(self, "mask_bad_ert", True):
            return None
        cycle_obs = list(self.dict_obs.items())[self.count_DA_cycle][1]
        segments, any_ert = [], False
        for key in obskey2map:
            entry = cycle_obs[key]
            df = entry["data"]
            if "ERT" in key and "rhoa" in getattr(df, "columns", []):
                keep_k, reasons = ERT_checks.quadrupole_quality_mask(df)
                any_ert = True
                if reasons:
                    logger.warning(
                        "ERT %s: %d/%d quadrupoles excluded from the analysis %s",
                        key, int((~keep_k).sum()), keep_k.size, reasons,
                    )
                # flag, never delete: the data stay complete in dict_obs and in the outputs,
                # this record says which quadrupoles the filter did not use, and why
                self.ert_check_reports.append(
                    {
                        "stage": "mask",
                        "cycle": self.count_DA_cycle,
                        "sensor": key,
                        "n_total": int(keep_k.size),
                        "n_excluded": int((~keep_k).sum()),
                        "reasons": reasons,
                        "excluded_idx": np.flatnonzero(~keep_k).tolist(),
                    }
                )
                segments.append(keep_k)
            else:
                segments.append(np.ones(int(np.size(entry["data"])), dtype=bool))
        if not any_ert:
            return None
        keep = np.concatenate(segments)
        n_before = keep.size
        if obs_valid is not None and np.size(obs_valid) == n_before:
            keep = keep[np.asarray(obs_valid, dtype=bool)]
        if keep.size != n_data:
            logger.warning(
                "ERT quality mask (%d) does not match the data (%d): no quadrupole excluded",
                keep.size, n_data,
            )
            return None
        return keep

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
        if np.ndim(data_cov) == 1:
            # stored as the diagonal (make_data_cov(..., diagonal=True)): the dense matrix of
            # THIS cycle only is rebuilt, so every filter and diagnostic sees what it always saw
            data_cov = np.diag(np.asarray(data_cov, dtype=float))
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

        # non-positive values (already warned about above) are clipped to 0: they get
        # log10(eps) = -3 instead of NaN / -inf
        data_valid_log = np.log10(np.clip(data_valid, 0.0, None) + eps)
        prediction_valid_log = np.log10(np.clip(prediction_valid, 0.0, None) + eps)

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
            data_valid_log,
            bins=30,
            alpha=0.5,
            label="Observed",
            density=True,
        )
        ax.hist(
            prediction_valid_log.flatten(),
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
        logger.debug("Plotting COV matrices")
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
                logger.warning("Impossible to plot the COV matrices (too large)")
        except Exception as exc:  # plotting is diagnostic only: never abort the DA run
            logger.warning("Impossible to plot matrices (%s: %s)", type(exc).__name__, exc)

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

        logger.info("Run hydrological model (%s)", "parallel" if parallel else "sequential")
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
                    logger.debug("%s", result)
        # process each ensemble folder one by one
        # ----------------------------------------------------------------
        else:
            # Loop over ensemble realisations
            # ------------------------------------------------------------
            for ens_i in track(
                range(self.NENS), description="Run hydrological fwd model..."
            ):

                logger.debug("Run ensemble member %d/%d", ens_i + 1, self.NENS)
                member_dir = os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
                p = subprocess.run(
                    [callexe], text=True, capture_output=True, cwd=member_dir
                )

                if self.verbose:
                    logger.debug("%s", p.stdout)
                    logger.debug("%s", p.stderr)

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

        logger.debug("Check scenarii before analysis")

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
                self._record_rejection(
                    "after_run", [n], f"cannot read mbeconv ({type(exc).__name__}: {exc})"
                )
                continue
            if reason is not None:
                logger.warning("Ensemble member %d rejected: %s", n + 1, reason)
                self._record_rejection("after_run", [n], reason)
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
        t_end = np.round(self.parm["(TIMPRT(I),I=1,NPRT)"][-1])
        if np.round(last_time) < (t_end - self.parm["DELTAT"]):
            msg = (
                f"CATHY stopped at t={last_time:g} s before the end of the window "
                f"(last printing time {t_end:g} s)"
            )
            if last_time <= 10 * float(self.parm.get("DTMIN", 0) or 0) + 1e-12:
                msg += (
                    f"; the time step collapsed to DTMIN={self.parm.get('DTMIN')} right at the "
                    "start: the analysed initial state/parameters are probably non-physical "
                    "(check the ERT innovation / geometry reports)"
                )
            return "df_mbeconv['time']" + str(last_time) + " - " + msg
        return None

    def _record_rejection(self, stage: str, members, reason: str) -> None:
        """Add rejected members (0-based indices) to the ledger and save it.

        ``stage`` is ``"after_run"`` (the CATHY run of the window failed) or
        ``"after_analysis"`` (an analysed parameter is not physical). The ledger answers
        who was rejected, in which cycle and why, which a bare percentage does not.
        """
        for member in members:
            self.rejection_log.append(
                {
                    "cycle": self.count_DA_cycle,
                    "stage": stage,
                    "member": int(member) + 1,
                    "reason": str(reason),
                }
            )
        path = os.path.join(
            self.workdir, self.project_name, self.project_name + "_rejections.json"
        )
        try:
            with open(path, "w") as f:
                json.dump(self.rejection_log, f, indent=1)
        except OSError as exc:  # diagnostics only: never stop a run for them
            logger.debug("cannot save %s (%s)", path, exc)

    def rejection_summary(self) -> dict:
        """Members (1-based, sorted) rejected so far, per stage."""
        summary: dict = {}
        for rec in self.rejection_log:
            summary.setdefault(rec["stage"], set()).add(rec["member"])
        return {stage: sorted(members) for stage, members in summary.items()}

    def _raise_if_too_many_rejected(
        self, rejected: Sequence[bool], threshold_rejected: float
    ) -> None:
        """Stop the run when the share of rejected members reaches the threshold (%)."""
        if sum(rejected) >= (threshold_rejected / 100) * self.NENS:
            logger.error("new parameters (if updated)")
            logger.error("%s", self.dict_parm_pert)
            summary = self.rejection_summary()
            for rec in self.rejection_log:
                logger.error(
                    "rejected: cycle %s, stage %s, member %s: %s",
                    rec["cycle"], rec["stage"], rec["member"], rec["reason"],
                )
            raise ValueError(
                "% number of rejected ensemble is too high:"
                + str((sum(rejected) * 100) / (self.NENS))
                + f". Rejected members by stage: {summary}"
                + f" ({self.project_name}_rejections.json has the reasons)"
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

        logger.debug("Check scenarii post update")

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
                logger.warning("%s %s: %s (members %s)", label, name, values[bad], [int(i) + 1 for i in bad])
                for i in bad:
                    self._record_rejection("after_analysis", [i], f"{label} {name} = {values[i]:.6g}")
            id_nonvalid.update(bad.tolist())
        return set(np.arange(self.NENS)).difference(id_nonvalid)

    def _max_root_depth(self) -> float:
        """Deepest feasible root depth: the thickness of the 3d mesh."""
        depth = abs(min(self.grid3d["mesh3d_nodes"][:, -1]))
        logger.debug("deepest feasible root depth: %s", depth)
        return depth

    def _back_transform_analysis_param(self, list_update_parm, analysis_param):
        """Back-transform the analysed parameters, one column per parameter.

        ``analysis_param`` is (n_members, n_param) with the "St. var." entries
        of ``list_update_parm`` excluded. Returned with the same shape.
        """
        names = [pp for pp in list_update_parm if "St. var." not in pp]
        ap = np.asarray(analysis_param, dtype=float)
        if ap.size == 0 or not names:
            return analysis_param
        if ap.ndim != 2 or ap.shape[1] != len(names):
            logger.warning(
                "analysed parameters have shape %s for %d parameter(s): "
                "back-transformation skipped",
                ap.shape,
                len(names),
            )
            return analysis_param
        out = ap.copy()
        for j, name in enumerate(names):
            out[:, j] = self.transform_parameters([name], param=ap[:, j], back=True)[0]
        return out

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
                        # inverse of y = log((p - A) / (B - p))
                        param_trans = (A + B * np.exp(param)) / (1.0 + np.exp(param))
                    else:
                        param_trans = np.log((param - A) / (B - param))
                elif self.dict_parm_pert[pp[1]]["transf_type"] == "hyperbolic":
                    logger.debug("hyperbolic transformation")
                    # Y = arcsinh(U); back-transformation U = sinh(Y)
                    param_trans = np.sinh(param) if back else np.arcsinh(param)

                else:
                    param_trans = param
                param_new.append(param_trans)

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

        for ww in warnings_petro:
            logger.warning(ww)

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
            dict_parm_pert[key].get("layer_groups"),
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
        if pert_control_name == "layers":  # porosity is never a material parameter
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
        logger.debug("Update Archie parameters")

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
                shellprint_update=False,
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
        logger.debug("Update ensemble")

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
                ctx.cycle_nb, dict_parm_pert, ens_nb, key, False
            )

        # ic update (i.e. initial pressure head update)
        elif name.casefold() in "ic".casefold():
            self._update_ic_member(ctx.cycle_nb, dict_parm_pert, ens_nb, key, False)

        # saturated hydraulic conductivity update
        elif name.casefold() in "ks":
            ctx.ks = self._update_soil_params_from_perturbations(
                dict_parm_pert, ctx.ks, key_root, key, update_key, ens_nb, False
            )

        elif name.casefold() in "porosity":
            ctx.porosity = self._update_soil_params_from_perturbations(
                dict_parm_pert, ctx.porosity, key_root, key, update_key, ens_nb, False
            )
            self.update_archie_parameter(
                _ARCHIE_PARM_NAMES, key_root, update_key, ens_nb, dict_parm_pert, key
            )

        # VG parameters update
        elif name in _VG_PARM_NAMES:
            ctx.vg = self._update_soil_params_from_perturbations(
                dict_parm_pert, ctx.vg, key_root, key, update_key, ens_nb, False
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
        _set_feddes_value(
            ctx.feddes[ens_nb], key_root[0], int(key_root[1]), new_Feddes_parm[ens_nb]
        )

        self.update_soil(
            FP_map=ctx.feddes[ens_nb],
            SPP_map=df_SPP,
            zone3d=ctx.zone3d,
            verbose=self.verbose,
            filename=os.path.join(os.getcwd(), "input/soil"),
            shellprint_update=False,
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
        logger.debug("Assimilated Observations are: %s", list_assimilated_obs)

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
            _log_once(
                "da_time_default",
                "No assimilation time given: selecting the observations of the current DA cycle",
            )

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
            ert_obs = list(self.dict_obs.values())[self.count_DA_cycle].get("ERT", {})
            logger.info(
                "Mapping states to ERT apparent resistivities (fwd noise level: %s%%)",
                ert_obs.get("fwdNoiseLevel", "n/a"),
            )
            logger.debug(
                "Archie perturbation for DA analysis: rFluid=%s, pert_sigma=%s, nb of zones=%d",
                self.Archie_parms["rFluid_Archie"],
                self.Archie_parms["pert_sigma_Archie"],
                len(self.Archie_parms["rFluid_Archie"]),
            )

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
        for rep in mapper.pop_ert_reports():
            self.ert_check_reports.append({"stage": "forward", **rep})
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

        if self.count_DA_cycle > 0:
            # one row of analysis_param_valid per parameter, "St. var." excluded
            # (the previous counter also counted "St. var." and so read row -1
            # for the first parameter when no "St. var." entry came first)
            index_update = 0
            for pp in list_update_parm:
                if "St. var." in pp:
                    continue
                update_key = _update_key(self.count_DA_cycle)
                self.dict_parm_pert[pp][update_key] = self._clip_updated_parameter(
                    pp, analysis_param_valid[index_update]
                )
                index_update = index_update + 1

            # check after analysis
            # ----------------------------------------------------------------
            id_valid_after = self._check_after_analysis(update_key, list_update_parm)
            intersection_set = set.intersection(
                set(id_valid_after), set(self.ens_valid)
            )
            # sorted: a set has no guaranteed order, and ens_valid is used to index arrays
            self.ens_valid = sorted(intersection_set)

    def _clip_updated_parameter(self, name: str, values) -> np.ndarray:
        """Project the analysed values of a parameter onto its bounds.

        The bounds are the ``clip_min`` / ``clip_max`` of ``dict_parm_pert[name]``
        (scenario ``per_bounds``), already used to sample the prior. The analysis
        is linear and unbounded: without this projection it can return values such
        as a negative root depth. Clipping keeps the member alive (rejecting it
        shrinks the ensemble); NaN (invalid members) are left untouched.
        """
        values = np.array(values, dtype=float)
        entry = self.dict_parm_pert.get(name, {})
        lo, hi = entry.get("clip_min"), entry.get("clip_max")
        if lo is None and hi is None:
            return values
        clipped = np.clip(values, lo, hi)
        with np.errstate(invalid="ignore"):
            changed = np.where(clipped != values)[0]
        if len(changed):
            logger.warning(
                "%s: %d analysed value(s) clipped to [%s, %s] (members %s): %s -> %s",
                name, len(changed), lo, hi,
                [int(i) + 1 for i in changed], values[changed], clipped[changed],
            )
        return clipped

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

        logger.debug("Update input ensemble")

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
                    logger.warning("atmbc perturbation not supported for spatially variable atmbc: ignored")

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
                shellprint_update=False,
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
                shellprint_update=False,
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
            logger.warning("Could not copy processor files into cathy_origin (%s: %s)", type(exc).__name__, exc)

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
