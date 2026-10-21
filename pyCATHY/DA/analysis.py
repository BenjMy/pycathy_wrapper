"""Analysis step of the data assimilation: result container, configuration, and the registry
of filters (``ANALYSES``) with their adapters. Moved out of ``cathy_DA`` without change; every
name is still importable from ``pyCATHY.DA.cathy_DA``."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Sequence
import logging
import warnings

import numpy as np

from pyCATHY.DA import enkf, pf

logger = logging.getLogger(__name__)


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
        Solve the analysis in ensemble space (exact, cheap for many data) where the
        filter supports it.
    perturb_obs : bool
        Stochastic EnKF: every member assimilates the data plus a draw of N(0, R).
        False (default, historical) gives an analysis spread smaller than the
        posterior by K R K^T and the ensemble tends to collapse over the cycles;
        True is recommended for ``enkf_Evensen2009*``. Seeded by ``DAConfig.seed``.
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
    pf_temper_neff_frac : float or None
        Temper the likelihood (weights ~ exp(beta*loglik)) so that N_eff >= frac * N. None: exact
        likelihood (degenerates to N_eff = 1 with thousands of ERT data).
    pf_liu_west_a : float or None
        Liu-West kernel shrinkage factor (e.g. 0.98) for the parameters after resampling; replaces
        the additive ``pf_jitter_std_param``. None: legacy jitter.
    """

    sakov: bool = True
    perturb_obs: bool = False
    inflation_states: float = 1.0
    inflation_params: float = 1.0
    jitter_params: float = 0.0
    localisation_matrix: Optional[np.ndarray] = None
    pf_resample_threshold: float = 0.5
    pf_jitter_std_param: float = 0.01
    pf_jitter_std_state: float = 0.05
    pf_temper_neff_frac: Optional[float] = None
    pf_liu_west_a: Optional[float] = None

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
            "perturb_obs",
            "jitter_params",
            "pf_resample_threshold",
            "pf_jitter_std_param",
            "pf_jitter_std_state",
            "pf_temper_neff_frac",
            "pf_liu_west_a",
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
    "log_transform_cov",
    "mask_bad_ert",
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
    log_transform_cov : bool
        With ``use_log_transformed_obs``: propagate the data covariance to the
        log10 domain (first-order, R_log = J R J, J = diag(1/(ln10 (d+eps)))).
        False reproduces the old behaviour (linear-domain R applied to log data,
        which makes the observations effectively ignored). Default True.
    mask_bad_ert : bool
        Drop from the analysis the ERT quadrupoles flagged unusable (k outliers or
        too large, non-positive rhoa, invalid flag; see
        ``ERT_checks.quadrupole_quality_mask``). They are removed from the data, R
        and predictions together; the full-length data are still used for the
        performance assessment. Default True.
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
    psi_bounds : (float or None, float or None) or None
        ``(psi_min, psi_max)`` [m]: the ANALYSED pressure heads are clipped to this range
        before they are written as the initial condition of the next window (None on one
        side = no bound on that side; None = no clipping, default). Use e.g.
        ``(-1e3, 0.0)`` for a purely unsaturated column: a large innovation can otherwise
        push the whole ensemble to saturation (resistivity floor), which collapses the
        spread and can make the next CATHY window fail (time step collapse). The number of
        clipped values is logged at every cycle.
    """

    inflation: float = 1.0
    localisation: Optional[str] = None
    threshold_rejected: int = 10
    use_log_transformed_obs: bool = False
    log_transform_cov: bool = True
    mask_bad_ert: bool = True
    cleanup_vtk: object = False
    backup_fort777: bool = False
    analysis: Optional[AnalysisConfig] = None
    seed: Optional[int] = None
    save_results: bool = False
    verbose: bool = False
    psi_bounds: Optional[tuple] = None

    def __post_init__(self):
        if self.psi_bounds is not None:
            if len(self.psi_bounds) != 2:
                raise ValueError(
                    f"psi_bounds must be (psi_min, psi_max), got {self.psi_bounds!r}."
                )
            lo, hi = self.psi_bounds
            if lo is not None and hi is not None and lo > hi:
                raise ValueError(f"psi_bounds: psi_min {lo} > psi_max {hi}.")
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


def _enkf(filter_function, forward_alpha=False, forward_perturb=False, **fixed):
    """Adapter for the filters called as f(data, cov, param, state, prediction, **kw)."""

    def adapter(data, data_cov, param, state, prediction, config):
        kwargs = dict(fixed)
        if forward_alpha:
            kwargs["alpha"] = config.inflation_states
        if forward_perturb:
            kwargs["perturb_obs"] = config.perturb_obs
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
        perturb_obs=config.perturb_obs,
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
        temper_neff_frac=config.pf_temper_neff_frac,
        liu_west_a=config.pf_liu_west_a,
    )
    logger.info(f"Final N_eff: {analysis['n_eff']:.1f} (exact likelihood: {analysis['n_eff_raw']:.2f}, beta = {analysis['beta']:.2e})")
    logger.info(f"Resampled: {analysis['resampled']}")
    return AnalysisResult(
        analysis=analysis["Analysis"],
        analysis_param=analysis["Analysisparam"].T,
        extra={
            k: analysis[k] for k in ("weights", "n_eff", "n_eff_raw", "beta", "resampled", "observation")
        },
    )


ANALYSES = {
    "enkf_Evensen2009": _enkf(enkf.enkf_analysis, forward_alpha=True, forward_perturb=True),
    "enkf_Evensen2009_Sakov": _enkf(
        enkf.enkf_analysis, forward_alpha=True, forward_perturb=True, Sakov=True
    ),
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
