"""Ensemble Smoother with Multiple Data Assimilation (ES-MDA).

Emerick & Reynolds (2013), Computers & Geosciences 55, 3-15.

The same data are assimilated ``n_iter`` times with an inflated observation error
``alpha_i * R``; the forward model is re-run from t = 0 between iterations. The
inflation coefficients must satisfy ``sum_i 1/alpha_i = 1`` so that, for a linear
Gaussian problem, the result equals one Kalman update with the original ``R``.

This module only holds the numerics (no CATHY run): the driver that loops the
forward model lives in ``cathy_DA.DA.run_DA_esmda``.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from pyCATHY.DA import enkf

logger = logging.getLogger(__name__)

_SCHEDULES = ("constant", "geometric", "user")
OBS_TIME_CONVENTIONS = ("state_time", "window_start")


def compute_alphas(
    n_iter: int,
    schedule: str = "constant",
    user_alphas: Optional[Sequence[float]] = None,
    geometric_ratio: float = 0.7,
    tol: float = 1e-8,
) -> np.ndarray:
    """Inflation coefficients ``alpha_i`` with ``sum(1/alpha_i) == 1``.

    Parameters
    ----------
    n_iter : int
        Number of assimilations (>= 1).
    schedule : {"constant", "geometric", "user"}
        ``constant``: alpha_i = n_iter. ``geometric``: alpha_i proportional to
        ``geometric_ratio**i`` (decreasing), rescaled to satisfy the constraint.
        ``user``: ``user_alphas`` as given (must satisfy the constraint within ``tol``).
    geometric_ratio : float
        Ratio between successive alphas for the geometric schedule, in (0, 1].
    """
    if not isinstance(n_iter, (int, np.integer)) or n_iter < 1:
        raise ValueError(f"n_iter must be an integer >= 1, got {n_iter!r}.")
    if schedule not in _SCHEDULES:
        raise ValueError(f"schedule must be one of {_SCHEDULES}, got {schedule!r}.")

    if schedule == "constant":
        alphas = np.full(n_iter, float(n_iter))
    elif schedule == "geometric":
        if not 0.0 < geometric_ratio <= 1.0:
            raise ValueError(f"geometric_ratio must be in (0, 1], got {geometric_ratio}.")
        raw = geometric_ratio ** np.arange(n_iter, dtype=float)
        # sum(1/alpha) = 1  with alpha_i = c * raw_i  ->  c = sum(1/raw)
        alphas = raw * np.sum(1.0 / raw)
    else:
        if user_alphas is None:
            raise ValueError("schedule='user' requires user_alphas.")
        alphas = np.asarray(user_alphas, dtype=float)
        if alphas.shape != (n_iter,):
            raise ValueError(
                f"user_alphas must have n_iter={n_iter} entries, got shape {alphas.shape}."
            )
        if np.any(alphas <= 0):
            raise ValueError("user_alphas must be strictly positive.")

    if abs(np.sum(1.0 / alphas) - 1.0) > tol:
        raise ValueError(
            f"sum(1/alpha) must be 1, got {np.sum(1.0 / alphas):.10f} for alphas={alphas}."
        )
    return alphas


@dataclass
class ESMDAConfig:
    """Options of an ES-MDA run.

    Attributes
    ----------
    n_iter : int
        Number of data assimilations.
    alpha_schedule : {"constant", "geometric", "user"}
        See ``compute_alphas``.
    user_alphas : sequence of float or None
        Coefficients when ``alpha_schedule='user'``.
    geometric_ratio : float
        Ratio of the geometric schedule.
    seed : int or None
        Seed of the generator used for the observation perturbations.
    update_initial_conditions : bool
        Also update the initial state (default False: parameters only).
    inflation_params : float
        Multiplicative inflation of the parameter anomalies applied after each update.
    localisation_matrix : np.ndarray or None
        (n_param x n_obs) taper applied to the cross covariance; None = no localisation.
    cleanup_between_iterations : bool
        Delete heavy model output of an iteration once its predictions are read.
    obs_time_convention : {"state_time", "window_start"}
        Which model state an observation stamped ``t`` is compared with (CATHY adapter only).
        ``"state_time"`` (default): the state at ``t`` itself, which is what
        ``run_DA_sequential`` does since commit d37036f9 (the analysis is made when the window
        ENDS at ``t``). ``"window_start"``: the state at the atmbc time that follows ``t``, the
        convention of ``run_DA_sequential`` before that commit (and of
        ``obs_window_start`` helpers that stamp a survey with the atmbc time just before it).
        Use the one your observations were stamped for.
    """

    n_iter: int = 4
    alpha_schedule: str = "constant"
    user_alphas: Optional[Sequence[float]] = None
    geometric_ratio: float = 0.7
    seed: Optional[int] = None
    update_initial_conditions: bool = False
    inflation_params: float = 1.0
    localisation_matrix: Optional[np.ndarray] = None
    cleanup_between_iterations: bool = False
    obs_time_convention: str = "state_time"

    def __post_init__(self):
        if self.obs_time_convention not in OBS_TIME_CONVENTIONS:
            raise ValueError(
                f"obs_time_convention must be one of {OBS_TIME_CONVENTIONS}, "
                f"got {self.obs_time_convention!r}."
            )
        # raises ValueError on any invalid n_iter / schedule / user_alphas
        self.alphas = compute_alphas(
            self.n_iter, self.alpha_schedule, self.user_alphas, self.geometric_ratio
        )
        if self.inflation_params <= 0:
            raise ValueError(f"inflation_params must be > 0, got {self.inflation_params}.")
        if self.localisation_matrix is not None:
            self.localisation_matrix = np.asarray(self.localisation_matrix)
            if self.localisation_matrix.size == 0:
                self.localisation_matrix = None


def _noise(data_cov, n_samples, rng):
    """Draw N(0, R) columns, shape (m, n_samples), from the generator ``rng``."""
    R = np.asarray(data_cov, dtype=float)
    if R.ndim == 1:
        return np.sqrt(R)[:, None] * rng.standard_normal((R.size, n_samples))
    if np.count_nonzero(R) == np.count_nonzero(np.diag(R)):
        std = np.sqrt(np.diag(R))
        return std[:, None] * rng.standard_normal((std.size, n_samples))
    return np.linalg.cholesky(R) @ rng.standard_normal((R.shape[0], n_samples))


def esmda_update(
    param,
    predicted,
    data,
    data_cov,
    alpha,
    rng=None,
    localisation=None,
    inflation=1.0,
):
    """One ES-MDA update of the parameter ensemble.

    Parameters
    ----------
    param : (n_param, n_ens) parameter ensemble.
    predicted : (n_obs, n_ens) forward-model predictions ``g(m_j)``.
    data : (n_obs,) or (n_obs, 1) observations.
    data_cov : (n_obs,) variances, or (n_obs, n_obs) covariance R.
    alpha : float, inflation coefficient of this iteration.
    rng : np.random.Generator or None
        Source of the observation perturbations (fresh draw per call and per member).
    localisation : (n_param, n_obs) taper or None.
    inflation : float
        Multiplicative inflation of the updated parameter anomalies.

    Returns
    -------
    (n_param, n_ens) updated parameter ensemble.
    """
    param = np.asarray(param, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    data = np.asarray(data, dtype=float).reshape(-1, 1)
    n_param, n_ens = param.shape
    n_obs = predicted.shape[0]

    if predicted.shape[1] != n_ens:
        raise ValueError(
            f"param has {n_ens} members but predicted has {predicted.shape[1]}."
        )
    if data.shape[0] != n_obs:
        raise ValueError(f"data has {data.shape[0]} values but predicted has {n_obs} rows.")
    if n_ens < 2:
        raise ValueError("At least 2 ensemble members are needed.")
    if alpha <= 0:
        raise ValueError(f"alpha must be > 0, got {alpha}.")
    R = np.asarray(data_cov, dtype=float)
    if R.shape not in ((n_obs,), (n_obs, n_obs)):
        raise ValueError(
            f"data_cov must have shape ({n_obs},) or ({n_obs}, {n_obs}), got {R.shape}."
        )

    rng = np.random.default_rng() if rng is None else rng

    # fresh perturbed observations: d_j = d + sqrt(alpha) * R^(1/2) z_j
    perturbed = data + np.sqrt(alpha) * _noise(R, n_ens, rng)
    innovations = perturbed - predicted

    param_anom = param - param.mean(axis=1, keepdims=True)
    obs_anom = predicted - predicted.mean(axis=1, keepdims=True)

    # (S S^T/(N-1) + alpha R)^-1 innovations, in the cheaper of obs / ensemble space
    weighted = enkf._solve_obs_space(obs_anom, alpha * R, innovations, n_ens)

    cross_cov = param_anom @ obs_anom.T / (n_ens - 1)  # (n_param, n_obs)
    if localisation is not None:
        localisation = np.asarray(localisation, dtype=float)
        if localisation.shape != cross_cov.shape:
            raise ValueError(
                f"localisation must have shape {cross_cov.shape}, got {localisation.shape}."
            )
        cross_cov = cross_cov * localisation

    updated = param + cross_cov @ weighted
    if inflation != 1.0:
        mean = updated.mean(axis=1, keepdims=True)
        updated = mean + inflation * (updated - mean)
    return updated


# ---------------------------------------------------------------------------
# Iteration driver (independent of CATHY: the forward model is a callable)
# ---------------------------------------------------------------------------
@dataclass
class ESMDAResult:
    """Output of ``run_esmda``.

    Attributes
    ----------
    param : (n_param, n_ens)
        Final parameter ensemble (columns of rejected members keep their last value).
    valid : (n_ens,) bool
        Members still active at the end.
    alphas : (n_iter,) inflation coefficients used.
    param_history : list of (n_param, n_ens)
        ``n_iter + 1`` arrays: the prior, then the ensemble after every update.
    prediction_history : list of (n_obs, n_ens)
        Predictions of every forward run (``n_iter`` runs, plus the final one when
        ``final_forward``); rows/columns that were not usable are NaN.
    misfit_history : (n_runs,) array
        Reduced chi-square of every forward run, ``mean_j r_j^T R^-1 r_j / n_obs``
        over the valid members and usable observations (about 1 at convergence).
    valid_history : list of (n_ens,) bool
        Active members at the start of every update.
    param_names : list of str or None
        Name of every row of ``param`` (set by the CATHY adapter).
    """

    param: np.ndarray
    valid: np.ndarray
    alphas: np.ndarray
    param_history: list
    prediction_history: list
    misfit_history: np.ndarray
    valid_history: list
    param_names: Optional[list] = None


def weighted_misfit(data, data_cov, prediction):
    """Reduced chi-square ``mean_j r_j^T R^-1 r_j / n_obs`` of an ensemble.

    ``r_j = data - prediction[:, j]``; ``data_cov`` is (n_obs,) or (n_obs, n_obs).
    """
    prediction = np.asarray(prediction, dtype=float)
    resid = np.asarray(data, dtype=float).reshape(-1, 1) - prediction
    R = np.asarray(data_cov, dtype=float)
    if R.ndim == 1:
        weighted = resid / R[:, None]
    else:
        weighted = np.linalg.solve(R, resid)
    return float(np.mean(np.sum(resid * weighted, axis=0)) / prediction.shape[0])


def _param_summary(param, max_parameters=4):
    """``p0 1.23+/-0.4; p1 ...`` for the first parameters of an ensemble (n_param, n_valid)."""
    rows = [f"p{j} {param[j].mean():.4g}+/-{param[j].std():.3g}" for j in range(min(len(param), max_parameters))]
    if len(param) > max_parameters:
        rows.append(f"... (+{len(param) - max_parameters} more)")
    return "; ".join(rows)


def _subset_cov(data_cov, rows):
    R = np.asarray(data_cov, dtype=float)
    return R[rows] if R.ndim == 1 else R[np.ix_(rows, rows)]


def run_esmda(
    param,
    forward,
    data,
    data_cov,
    config=None,
    project=None,
    callback=None,
    final_forward=True,
    min_members=2,
):
    """Run the ES-MDA iterations around a user supplied forward model.

    Parameters
    ----------
    param : (n_param, n_ens)
        Prior parameter ensemble, in the space where the update is done (apply any
        log transform before, and undo it inside ``forward``).
    forward : callable ``forward(param, valid) -> (prediction, ok)``
        Runs the model from t = 0 for the members flagged in ``valid`` and returns the
        predicted observations ``prediction`` (n_obs, n_ens) and a boolean ``ok``
        (n_ens,) telling which members ran successfully. Columns of members that were
        not run (or failed) are ignored, NaN is fine.
    data : (n_obs,) observations, all assimilation times stacked.
    data_cov : (n_obs,) variances or (n_obs, n_obs) covariance R.
    config : ESMDAConfig or None (defaults).
    project : callable ``project(param) -> param`` or None
        Applied after every update (bounds, reflection, back-transform checks).
    callback : callable ``callback(iteration, alpha, param, prediction, valid)`` or None
        Called after every forward run (``iteration`` = 0 .. n_iter; the last call has
        ``alpha=None`` when ``final_forward``). Use it to save outputs.
    final_forward : bool
        Run the model once more with the final parameters (posterior predictions).
    min_members : int
        Raise ``RuntimeError`` if fewer members remain valid.

    Notes
    -----
    Members that fail are dropped for the rest of the run (as ``ens_valid`` in the
    sequential DA). Observations whose prediction is non-finite for any valid member
    are left out of that update (data, R and localisation rows together).
    """
    config = ESMDAConfig() if config is None else config
    param = np.array(param, dtype=float)  # copy
    if param.ndim != 2:
        raise ValueError(f"param must be (n_param, n_ens), got shape {param.shape}.")
    n_param, n_ens = param.shape
    data = np.asarray(data, dtype=float).reshape(-1)
    R = np.asarray(data_cov, dtype=float)
    n_obs = data.size
    if R.shape not in ((n_obs,), (n_obs, n_obs)):
        raise ValueError(
            f"data_cov must have shape ({n_obs},) or ({n_obs}, {n_obs}), got {R.shape}."
        )
    loc = config.localisation_matrix
    if loc is not None and loc.shape != (n_param, n_obs):
        raise ValueError(
            f"localisation_matrix must have shape ({n_param}, {n_obs}), got {loc.shape}."
        )
    if min_members < 2:
        raise ValueError("min_members must be >= 2.")

    rng = np.random.default_rng(config.seed)
    valid = np.ones(n_ens, dtype=bool)
    param_history = [param.copy()]
    prediction_history, misfit_history, valid_history = [], [], []

    def _run(iteration, alpha):
        nonlocal valid
        pred, ok = forward(param.copy(), valid.copy())
        pred = np.array(pred, dtype=float)
        ok = np.asarray(ok, dtype=bool)
        if pred.shape != (n_obs, n_ens):
            raise ValueError(
                f"forward returned predictions of shape {pred.shape}, expected {(n_obs, n_ens)}."
            )
        if ok.shape != (n_ens,):
            raise ValueError(f"forward returned ok of shape {ok.shape}, expected {(n_ens,)}.")
        failed = valid & ~ok
        if failed.any():
            logger.warning(
                "ES-MDA iteration %d: member(s) %s failed and are dropped.",
                iteration,
                [int(i) + 1 for i in np.flatnonzero(failed)],
            )
        valid = valid & ok
        if valid.sum() < min_members:
            raise RuntimeError(
                f"Only {int(valid.sum())} valid member(s) left at iteration {iteration} "
                f"(min_members={min_members})."
            )
        pred[:, ~valid] = np.nan
        obs_ok = np.all(np.isfinite(pred[:, valid]), axis=1)
        if not obs_ok.any():
            raise RuntimeError(f"No usable observation at iteration {iteration}.")
        misfit = weighted_misfit(
            data[obs_ok], _subset_cov(R, obs_ok), pred[np.ix_(obs_ok, valid)]
        )
        prediction_history.append(pred)
        misfit_history.append(misfit)
        if alpha is None:
            logger.info(
                "ES-MDA final run (posterior): %d/%d members valid | misfit %.4g (reduced chi-square)",
                int(valid.sum()), n_ens, misfit,
            )
        else:
            logger.info(
                "ES-MDA iteration %d/%d: forward run done | %d/%d members valid | misfit %.4g "
                "(reduced chi-square) | updating with alpha=%.4g",
                iteration + 1, len(config.alphas), int(valid.sum()), n_ens, misfit, alpha,
            )
        if callback is not None:
            callback(iteration, alpha, param.copy(), pred.copy(), valid.copy())
        return pred, obs_ok

    for i, alpha in enumerate(config.alphas):
        pred, obs_ok = _run(i, float(alpha))
        valid_history.append(valid.copy())
        updated = esmda_update(
            param[:, valid],
            pred[np.ix_(obs_ok, valid)],
            data[obs_ok],
            _subset_cov(R, obs_ok),
            float(alpha),
            rng=rng,
            localisation=None if loc is None else loc[:, obs_ok],
            inflation=config.inflation_params,
        )
        param[:, valid] = updated
        if project is not None:
            projected = np.asarray(project(param.copy()), dtype=float)
            if projected.shape != param.shape:
                raise ValueError(
                    f"project returned shape {projected.shape}, expected {param.shape}."
                )
            param = projected
        param_history.append(param.copy())
        logger.info(
            "ES-MDA iteration %d/%d: updated parameters (mean +/- sd of the valid members): %s",
            i + 1, len(config.alphas), _param_summary(param[:, valid]),
        )

    if final_forward:
        _run(len(config.alphas), None)

    return ESMDAResult(
        param=param,
        valid=valid,
        alphas=config.alphas.copy(),
        param_history=param_history,
        prediction_history=prediction_history,
        misfit_history=np.asarray(misfit_history),
        valid_history=valid_history,
    )
