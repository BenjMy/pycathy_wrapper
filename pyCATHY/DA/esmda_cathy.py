"""ES-MDA on a CATHY ensemble: the glue between ``esmda.run_esmda`` and the ``DA`` object.

.. warning::
    EXPERIMENTAL. The orchestration is unit tested with a stub of the simulation object
    (``tests/test_esmda_cathy.py``) but it has **not** been run against the CATHY
    executable. Check a small case (few members, short window) before trusting results.

What it does (``run_esmda_on_cathy``)
-------------------------------------
Every iteration runs the whole simulation from t = 0 (the initial condition is never
updated as a *state*; an initial-condition *parameter* such as ``ic`` is updated like any
other parameter), reads the state at **every** observation time and stacks the predictions
in the order of the data. Supported observations: tensiometers (``tensio``), soil water
content (``swc``) and ERT (``ERT``, apparent resistivities through the Archie operator).
EM and ET observations raise ``NotImplementedError``.

For ERT the data go through the same steps as in the sequential analysis: unusable
quadrupoles are dropped (``DAConfig.mask_bad_ert``, from ``ERT_checks.quadrupole_quality_mask``),
then data and predictions are put in log10 and R is propagated to log10
(``DAConfig.use_log_transformed_obs`` / ``log_transform_cov``).

Which state is compared with an observation stamped ``t`` is set by
``ESMDAConfig.obs_time_convention`` (see there).

This module must not import ``cathy_DA`` (it gets the simulation as an argument).
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from pyCATHY.DA import esmda, mapper
from pyCATHY.DA.param_tables import _update_key
from pyCATHY.ERT import ERT_checks
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT

logger = logging.getLogger(__name__)

SUPPORTED_OBSERVATIONS = ("tensio", "swc", "ERT")
LOG_EPS = 1e-3  # offset added before log10: the value used by the sequential analysis
# reading failures that mean "this member's output is unusable" (truncated or missing files).
# KeyError / AttributeError are NOT here: they point to a configuration or code problem and
# must stop the run instead of silently dropping every member.
_MEMBER_READ_ERRORS = (OSError, IndexError, ValueError)


class ESMDATimeError(ValueError):
    """An observation time is not among the printed output times (configuration error)."""


# --------------------------------------------------------------------------- pure helpers
# CATHY prints the state at the FIRST TIME STEP >= a requested printing time (it does not shorten the step to land
# on it). In ONE run over the whole period (ES-MDA) a survey time is therefore printed up to DTMAX seconds late.
# PRINT_LAG_TOL [s] accepts a printed time in [t, t + PRINT_LAG_TOL] as the state at t. 0 (default) = exact match,
# the historical behaviour. Set it to a few DTMAX when the time steps are small (e.g. DTMAX = 60 s -> 120 s).
PRINT_LAG_TOL = 0.0


def state_rows_at_times(df, times, atol=1e-6, lag_tol=None):
    """Positions of the rows of ``df`` (index = output time) matching ``times``.

    A row matches if its time is within ``atol`` (relative to 1 + |t|) of ``t`` or, when ``lag_tol`` (default
    ``PRINT_LAG_TOL``) is positive, if it is the first row printed in ``[t, t + lag_tol]`` (CATHY prints at the
    first time step >= the requested time). Raises ``ESMDATimeError`` if a time has no such row.
    """
    lag_tol = PRINT_LAG_TOL if lag_tol is None else float(lag_tol)
    index = np.asarray(df.index, dtype=float)
    rows = []
    for t in times:
        hits = np.flatnonzero(np.abs(index - float(t)) <= atol * (1.0 + abs(float(t))))
        if hits.size == 0 and lag_tol > 0:
            late = np.flatnonzero((index >= float(t)) & (index <= float(t) + lag_tol))
            hits = late[:1]
        if hits.size == 0:
            raise ESMDATimeError(
                f"No output at t={t}: printed times are {index.tolist()}. "
                "Observation times must be CATHY printing times (or set esmda_cathy.PRINT_LAG_TOL)."
            )
        rows.append(int(hits[0]))
    return rows


def block_diag_cov(covs):
    """Block-diagonal observation covariance from a list of per-time covariances.

    Each block is a 1-D array of variances or a 2-D covariance. Returns a 1-D vector of
    variances when every block is diagonal, otherwise the dense block-diagonal matrix.
    """
    blocks = [np.atleast_1d(np.asarray(c, dtype=float)) for c in covs]
    if not blocks:
        raise ValueError("block_diag_cov needs at least one block.")
    for b in blocks:
        if b.ndim == 2 and b.shape[0] != b.shape[1]:
            raise ValueError(f"covariance block must be square, got {b.shape}.")
        if b.ndim > 2:
            raise ValueError(f"covariance block must be 1-D or 2-D, got {b.ndim}-D.")
    diagonal = []
    for b in blocks:
        if b.ndim == 1:
            diagonal.append(b)
        elif np.count_nonzero(b) == np.count_nonzero(np.diag(b)):
            diagonal.append(np.diag(b))
        else:
            diagonal = None
            break
    if diagonal is not None:
        return np.concatenate(diagonal)
    sizes = [b.shape[0] for b in blocks]
    out = np.zeros((sum(sizes), sum(sizes)))
    start = 0
    for b, n in zip(blocks, sizes):
        out[start : start + n, start : start + n] = np.diag(b) if b.ndim == 1 else b
        start += n
    return out


def log10_obs(x, eps=LOG_EPS):
    """``log10(max(x, 0) + eps)``, as ``DA._log_transform_obs`` (NaN stays NaN)."""
    return np.log10(np.clip(np.asarray(x, dtype=float), 0.0, None) + eps)


def log10_cov(data, cov, eps=LOG_EPS):
    """R propagated to log10 (first order), as ``cathy_DA._log10_cov``.

    ``cov`` is a 1-D vector of variances or a 2-D covariance of ``data``.
    """
    d = np.clip(np.asarray(data, dtype=float), 0.0, None) + eps
    jac = 1.0 / (np.log(10.0) * d)
    cov = np.asarray(cov, dtype=float)
    return cov * jac**2 if cov.ndim == 1 else cov * np.outer(jac, jac)


def _mask_cov(cov, keep):
    cov = np.asarray(cov, dtype=float)
    return cov[keep] if cov.ndim == 1 else cov[np.ix_(keep, keep)]


def resolve_parameter_names(dict_parm_pert, list_parm2update):
    """Keys of ``dict_parm_pert`` to update.

    Each entry of ``list_parm2update`` is an exact key (``"ic"``, ``"SATCELL_VG0"``) or the
    root of numbered keys (``"Ks"`` -> ``"Ks0"``, ``"Ks1"``, ...). It is NOT a substring
    match: updating ``"ic"`` must not pick ``"a_Archie"``. ``"St. var."`` is dropped (states
    are not updated by ES-MDA); ``"all"`` selects every key. An entry that matches nothing
    raises (it would silently not be updated).
    """
    if isinstance(list_parm2update, str):
        list_parm2update = [list_parm2update]
    wanted = [p for p in list_parm2update if "St. var." not in p]
    if "all" in wanted:
        names = list(dict_parm_pert)
    else:
        names, unmatched = [], []
        for w in wanted:
            if w in dict_parm_pert:
                matches = [w]
            else:
                matches = [k for k in dict_parm_pert if re.fullmatch(re.escape(w) + r"\d+", k)]
            if not matches:
                unmatched.append(w)
            names += [m for m in matches if m not in names]
        if unmatched:
            raise ValueError(
                f"list_parm2update entries {unmatched} are not perturbed parameters "
                f"(keys of dict_parm_pert: {list(dict_parm_pert)})."
            )
    if not names:
        raise ValueError(
            f"No parameter to update for list_parm2update={list(list_parm2update)}: "
            "ES-MDA updates parameters ('St. var.' alone is not enough)."
        )
    return names


def _entry_problem(key, entry):
    """Why an observation entry cannot be used, or None."""
    if not isinstance(entry, dict) or "data" not in entry:
        return "needs 'data'"
    if "ERT" in key:
        if "data_format" not in entry:
            return "needs 'data' and 'data_format' (build it with read_observations(data_type='ERT'))"
    elif "mesh_nodes" not in entry:
        return "needs 'data' and 'mesh_nodes' (the observed mesh nodes)"
    return None


def observation_plan(dict_obs, list_assimilated_obs):
    """``[(time stamp, [sensor keys], [obs dicts])]`` in time order, for the assimilated sensors.

    Same selection rule as ``DA._obs_key_select`` / ``_get_data2assimilate``. Raises
    ``NotImplementedError`` for observation types the smoother does not support and
    ``ValueError`` for incomplete entries.
    """
    if isinstance(list_assimilated_obs, str):
        list_assimilated_obs = [list_assimilated_obs]
    plan = []
    for time, sensors in dict_obs.items():
        keys, obs = [], []
        for sensor, entry in sensors.items():
            root = re.sub(r"\d", "", sensor)
            if "all" in list_assimilated_obs or root in list_assimilated_obs:
                keys.append(sensor)
                obs.append(entry)
        for key in keys:
            if not any(s in key for s in SUPPORTED_OBSERVATIONS):
                raise NotImplementedError(
                    f"ES-MDA on CATHY supports only {SUPPORTED_OBSERVATIONS} observations, "
                    f"got {key!r} at t={time}."
                )
        for key, entry in zip(keys, obs):
            problem = _entry_problem(key, entry)
            if problem:
                raise ValueError(f"Observation {key!r} at t={time} {problem}.")
        if keys:
            plan.append((float(time), keys, obs))
    if not plan:
        raise ValueError(f"No observation selected by list_assimilated_obs={list_assimilated_obs}.")
    return plan


def state_times(plan, atm_times, convention="state_time"):
    """Absolute time of the model state each observation of ``plan`` is compared with.

    ``"state_time"``: the stamp itself. ``"window_start"``: the atmbc time that follows the
    stamp (the stamp must be an atmbc time, and not the last one). See
    ``ESMDAConfig.obs_time_convention``.
    """
    stamps = [t for t, _, _ in plan]
    if convention == "state_time":
        return list(stamps)
    if convention != "window_start":
        raise ValueError(f"unknown obs_time_convention {convention!r}.")
    atm = np.asarray(atm_times, dtype=float)
    out = []
    for t in stamps:
        hit = np.flatnonzero(np.abs(atm - t) <= 1e-6 * (1.0 + abs(t)))
        if hit.size == 0:
            raise ValueError(
                f"obs_time_convention='window_start': t={t} is not an atmbc time "
                f"(atmbc times: {atm.tolist()})."
            )
        if hit[0] + 1 >= atm.size:
            raise ValueError(
                f"obs_time_convention='window_start': t={t} is the last atmbc time, "
                "no model window starts there."
            )
        out.append(float(atm[hit[0] + 1]))
    return out


def check_state_times(plan, atm_times, convention="state_time", acquisition_times=None, rtol=1e-6):
    """Problems (list of messages, empty = fine) pairing each observation with a model state.

    Every observation must be compared with a state that is an atmbc time, and, if
    ``acquisition_times`` is given (one per plan entry, same order), with the state at its
    acquisition time. Times are compared with a RELATIVE tolerance ``rtol`` (float32 times read
    back from the atmbc file differ from the float64 ones by up to ~1e-3 s at 1e4 s and
    ~0.06 s at 1e6 s): an absolute 1e-6 s tolerance would reject them.
    """
    atm = np.asarray(atm_times, dtype=float)
    try:
        states = state_times(plan, atm, convention)
    except ValueError as exc:
        return [str(exc)]
    problems = []
    for n, ((stamp, _, _), st) in enumerate(zip(plan, states)):
        tol = rtol * (1.0 + abs(st))
        if not (np.abs(atm - st) <= tol).any():
            nearest = float(atm[np.argmin(np.abs(atm - st))])
            problems.append(
                f"survey stamped {stamp:.3f} s ({stamp / 3600:.2f} h): the state at {st:.3f} s is not an "
                f"atmbc time (nearest {nearest:.3f} s). If the atmbc vector is a float32 read-back: "
                "DA.read_inputs('atmbc') overwrites sim.atmbc with it; keep a copy of sim.atmbc and "
                "restore it after the call."
            )
        if acquisition_times is not None and abs(st - float(acquisition_times[n])) > tol:
            problems.append(
                f"survey acquired at {float(acquisition_times[n]):.3f} s is compared with the state at "
                f"{st:.3f} s: the stamps do not follow obs_time_convention={convention!r}."
            )
    return problems


def uses_swc(plan):
    """True if the plan contains soil water content observations."""
    return any("swc" in key for _, keys, _ in plan for key in keys)


def uses_ert(plan):
    """True if the plan contains ERT observations."""
    return any("ERT" in key for _, keys, _ in plan for key in keys)


# --------------------------------------------------------------------------- observations
@dataclass
class ObsStack:
    """Data of all assimilation times as the filter sees them.

    ``data`` / ``cov`` are masked (usable quadrupoles) and log10-transformed when
    requested; ``keep`` is the mask over the RAW stacked layout (``layout``) and
    ``use_log`` tells the predictions to follow the same steps (``transform_prediction``).
    """

    data: np.ndarray
    cov: np.ndarray
    keep: np.ndarray
    layout: List[tuple]  # (plan index, sensor key, n raw data)
    use_log: bool

    def transform_prediction(self, pred_raw):
        pred = np.asarray(pred_raw, dtype=float)[self.keep]
        return log10_obs(pred) if self.use_log else pred


def _entry_values(key, entry):
    data = entry["data"]
    if "ERT" in key:
        col = "rhoa" if "pygimli" in str(entry.get("data_format")) else "resist"
        return np.asarray(data[col], dtype=float).ravel()
    return np.atleast_1d(np.asarray(data, dtype=float)).ravel()


def stack_data(sim, plan, list_assimilated_obs=None):
    """Observation vector, covariance and masks of all times of ``plan`` (see ``ObsStack``).

    The values are read from the observation entries (``_entry_values``, the same extraction as
    ``DA._get_data2assimilate``, checked against it in the tests where the library can do it):
    ``DA._get_data2assimilate`` cannot stack sensors of different lengths in one cycle (it
    squeezes a ragged list), the smoother needs the per-sensor segments anyway.
    ``list_assimilated_obs`` is not used (the selection is in ``plan``), kept for symmetry.
    """
    index_of = {float(t): i for i, t in enumerate(sim.dict_obs)}
    use_mask = bool(getattr(sim, "mask_bad_ert", True))
    use_log = bool(getattr(sim, "use_log_transformed_obs", False))
    log_cov = bool(getattr(sim, "log_transform_cov", True))

    layout, raw, keeps, covs_in = [], [], [], []
    for i, (t, keys, obs) in enumerate(plan):
        k = index_of[t]
        segments = [_entry_values(key, e) for key, e in zip(keys, obs)]
        sizes = [seg.size for seg in segments]
        values = np.concatenate(segments)
        cov = np.asarray(sim.stacked_data_cov[k], dtype=float)
        if cov.shape not in ((values.size,), (values.size, values.size)):
            raise ValueError(
                f"stacked_data_cov[{k}] has shape {cov.shape} for {values.size} data at t={t}."
            )
        keep_t = []
        for key, entry, n in zip(keys, obs, sizes):
            layout.append((i, key, n))
            k_seg = np.ones(n, dtype=bool)
            if use_mask and "ERT" in key and "rhoa" in getattr(entry["data"], "columns", []):
                k_mask, reasons = ERT_checks.quadrupole_quality_mask(entry["data"])
                if k_mask.size == n:
                    k_seg = np.asarray(k_mask, dtype=bool)
                    if not k_seg.all():
                        logger.warning(
                            "ERT %s at t=%s: %d/%d quadrupoles excluded %s",
                            key, t, int((~k_seg).sum()), n, reasons,
                        )
            keep_t.append(k_seg)
        keep_t = np.concatenate(keep_t)
        d_f, cov_f = values[keep_t], _mask_cov(cov, keep_t)
        if use_log and log_cov:
            cov_f = log10_cov(d_f, cov_f)
        raw.append(values)
        keeps.append(keep_t)
        covs_in.append(cov_f)

    keep = np.concatenate(keeps)
    data = np.concatenate(raw)[keep]
    if use_log:
        data = log10_obs(data)
    return ObsStack(data=data, cov=block_diag_cov(covs_in), keep=keep, layout=layout, use_log=use_log)


# --------------------------------------------------------------------------- CATHY readers
def _member_dir(sim, ens_nb):
    return os.path.join(sim.workdir, sim.project_name, f"DA_Ensemble/cathy_{ens_nb + 1}")


def predict_member_blocks(sim, ens_nb, plan, times, time_offset=0.0):
    """Point-sensor predictions (``tensio``, ``swc``) of one member: ``{(i, key): 1-D array}``.

    Same operators as ``DA._predict_member`` (tensiometer: the pressure head at the node;
    SWC: saturation times porosity) applied to the state of each of ``times`` (absolute
    times, one per plan entry) instead of the last one. ``time_offset`` is the absolute time
    of the first output. ERT is not handled here (``predict_ert_blocks``).
    """
    if not any(("tensio" in key or "swc" in key) for _, keys, _ in plan for key in keys):
        return {}
    path = _member_dir(sim, ens_nb)
    out_dir = os.path.join(path, sim.output_dirname)
    df_psi = sim.read_outputs(filename="psi", path=out_dir)
    df_sw, _ = sim.read_outputs(filename="sw", path=out_dir)

    rel = [t - time_offset for t in times]
    rows_psi = state_rows_at_times(df_psi, rel)
    rows_sw = state_rows_at_times(df_sw, rel)

    SPP_ensi = None
    if uses_swc(plan):
        SPP_ensi, _ = in_CT.read_soil(
            os.path.join(path, sim.input_dirname, "soil"),
            sim.dem_parameters,
            sim.cathyH["MAXVEG"],
        )

    blocks = {}
    for i, ((_, keys, obs), r_psi, r_sw) in enumerate(zip(plan, rows_psi, rows_sw)):
        state = [df_psi.iloc[r_psi], df_sw.iloc[r_sw].values]
        for key, entry in zip(keys, obs):
            if "tensio" in key:
                blocks[(i, key)] = np.atleast_1d(mapper.tensio_mapper(state, entry["mesh_nodes"]))
            elif "swc" in key:
                blocks[(i, key)] = np.atleast_1d(
                    mapper.swc_mapper(
                        state,
                        entry["mesh_nodes"],
                        sim.grid3d["mesh3d_nodes"],
                        sim.dem_parameters,
                        SPP_ensi,
                        sim.console,
                    )
                )
    return blocks


def predict_ert_blocks(sim, plan, times, ok_members, time_offset=0.0):
    """ERT predictions of the ensemble: ``{(i, key): (n_quadrupoles, len(ok_members)) array}``.

    The survey is picked by its position in ``dict_obs`` and the saturation of the row of the
    member ``output/sw`` file at the state time (``time_ass`` of ``mapper._map_ERT_parallel``).
    Columns follow ``ok_members`` (the members listed in ``sim.ens_valid``).
    """
    items = [(i, key) for i, (_, keys, _) in enumerate(plan) for key in keys if "ERT" in key]
    if not items:
        return {}
    index_of = {float(t): i for i, t in enumerate(sim.dict_obs)}
    out_dir = os.path.join(_member_dir(sim, int(ok_members[0])), sim.output_dirname)
    df_sw, _ = sim.read_outputs(filename="sw", path=out_dir)
    rows = state_rows_at_times(df_sw, [times[i] - time_offset for i, _ in items])

    paths, porosity = sim._ert_ensemble_inputs()  # paths of sim.ens_valid, soil read per member
    blocks = {}
    for (i, key), row in zip(items, rows):
        hx, _, _ = mapper._map_ERT_parallel(
            sim.dict_obs,
            index_of[plan[i][0]],
            sim.project_name,
            sim.Archie_parms,
            porosity,
            paths,
            list_assimilated_obs="all",
            default_state="psi",
            DA_cnb=index_of[plan[i][0]],
            savefig=False,
            verbose=False,
            time_ass=row,
        )
        hx = np.asarray(hx, dtype=float)
        if hx.ndim != 2 or hx.shape[1] != len(ok_members):
            raise ValueError(
                f"ERT forward returned shape {hx.shape} for {len(ok_members)} members at t={plan[i][0]}."
            )
        blocks[(i, key)] = hx
        reports = getattr(sim, "ert_check_reports", None)
        for rep in mapper.pop_ert_reports():
            if reports is not None:
                reports.append({"stage": "forward", "plan_index": i, **rep})
    return blocks


def _tol(t, rtol=1e-6):
    return rtol * (1.0 + abs(float(t)))


def build_run_summary(*, names, prior_physical, n_members, config, plan, times, print_times, atm_times,
                      n_data, n_data_raw, use_log, threshold_rejected, bounds_mode, parallel):
    """What a run will do and cost, as a plain ``dict`` (also the result of ``dry_run=True``).

    ``windows`` lists, for every assimilation time, ``(stamp, state time, previous state time,
    number of atmbc times in the window)``: the first window starts at the first atmbc time, the
    next ones at the previous state time. Unlike the sequential DA the model is NOT restarted at
    these times: they are the instants at which the state of ONE run from the start is compared
    with the data.
    """
    atm = np.asarray(atm_times, dtype=float)
    t_start, t_end = float(atm[0]), float(atm[-1])
    n_ert = sum(1 for _, keys, _ in plan for key in keys if "ERT" in key)
    runs = config.n_iter + 1  # n_iter updates + the final run with the posterior
    windows, previous = [], t_start
    for (stamp, _, _), t in zip(plan, times):
        n_steps = int(np.sum((atm > previous + _tol(previous)) & (atm <= float(t) + _tol(t))))
        windows.append((float(stamp), float(t), float(previous), n_steps))
        previous = float(t)
    sensors = sorted({re.sub(r"\d", "", key) for _, keys, _ in plan for key in keys})
    return {
        "parameters": list(names),
        "prior": {n: (float(np.mean(prior_physical[n])), float(np.std(prior_physical[n]))) for n in names},
        "n_members": int(n_members),
        "n_iter": int(config.n_iter),
        "alpha_schedule": config.alpha_schedule,
        "alphas": [float(a) for a in config.alphas],
        "inflation_params": float(config.inflation_params),
        "localisation": config.localisation_matrix is not None,
        "seed": config.seed,
        "obs_time_convention": config.obs_time_convention,
        "sensors": sensors,
        "n_windows": len(plan),
        "surveys": [(float(p[0]), float(t)) for p, t in zip(plan, times)],  # (stamp, state time)
        "windows": windows,
        "period": (t_start, t_end),
        "n_atmbc_times": int(atm.size),
        "print_times": [float(t) for t in print_times],
        "n_data": int(n_data),
        "n_data_raw": int(n_data_raw),
        "n_data_masked": int(n_data_raw - n_data),
        "log10": bool(use_log),
        "threshold_rejected": threshold_rejected,
        "param_bounds_mode": bounds_mode,
        "parallel": bool(parallel),
        "n_model_runs": int(runs),
        "n_member_simulations": int(runs * n_members),
        "n_ert_forwards": int(runs * n_ert * n_members),
    }


def format_run_plan(summary, max_parameters=12):
    """The recap printed (and logged) before the first model run: lines of text."""
    sm = summary
    hrs = lambda t: t / 3600.0  # noqa: E731
    t0, t1 = sm["period"]
    names = sm["parameters"]
    shown = [f"{n} {sm['prior'][n][0]:.4g}\u00b1{sm['prior'][n][1]:.3g}" for n in names[:max_parameters]]
    if len(names) > max_parameters:
        shown.append(f"... (+{len(names) - max_parameters} more)")
    alphas = ", ".join(f"{a:.3f}" for a in sm["alphas"])
    inv = sum(1.0 / a for a in sm["alphas"])
    lines = [
        "[ES-MDA] ============================== strategy ==============================",
        f"  parameters       : {len(names)} estimated, prior (mean\u00b1sd, physical): " + "; ".join(shown),
        f"  iterations       : {sm['n_iter']} update(s) + 1 final run | alpha ({sm['alpha_schedule']}) = {alphas}"
        f" | sum(1/alpha) = {inv:.3f} | inflation {sm['inflation_params']:g}"
        f" | localisation {'on' if sm['localisation'] else 'off'} | seed {sm['seed']}",
        f"  members          : {sm['n_members']} | a member is dropped for good when its water-balance error"
        f" exceeds {sm['threshold_rejected']} | bounds mode: {sm['param_bounds_mode']}",
        f"  observations     : {'+'.join(sm['sensors'])} | {sm['n_windows']} assimilation windows | {sm['n_data']} data"
        f" ({sm['n_data_masked']} masked out of {sm['n_data_raw']}) | log10: {sm['log10']}",
        f"  time convention  : {sm['obs_time_convention']}",
        f"  simulated period : t = {t0:.0f} -> {t1:.0f} s ({hrs(t0):.2f} -> {hrs(t1):.2f} h, {(t1 - t0) / 86400:.2f} d)"
        f" | {sm['n_atmbc_times']} atmbc times ({sm['n_atmbc_times'] - 1} atmbc windows), the whole period is run"
        " from t = start at every iteration",
        f"  cost             : {sm['n_model_runs']} model runs x {sm['n_members']} members = "
        f"{sm['n_member_simulations']} simulations" + (f", {sm['n_ert_forwards']} ERT forwards" if sm["n_ert_forwards"] else "")
        + f" | parallel: {sm['parallel']}",
        "[ES-MDA] assimilation windows: instants at which the state of ONE run is compared with the data",
        "     #   stamp [h]   state [h]  window from [h]  length [h]  atmbc times in window",
    ]
    for i, (stamp, state, previous, n_steps) in enumerate(sm["windows"], start=1):
        lines.append(
            f"  {i:4d}  {hrs(stamp):10.2f}  {hrs(state):10.2f}  {hrs(previous):15.2f}  {hrs(state - previous):10.2f}"
            f"  {n_steps:21d}"
        )
    lines.append("[ES-MDA] ======================================================================")
    return lines


def ensure_grid3d(sim):
    """Load the 3d mesh of the project if needed (same dry run as ``selec_atmbc_window``).

    ``DA.prepare_DA`` does not load it, and the SWC operator needs ``mesh3d_nodes``.
    """
    if len(sim.grid3d) == 0:
        sim.run_processor(IPRT1=3, DAFLAG=0, verbose=True)
        sim.grid3d = out_CT.read_grid3d(
            os.path.join(sim.workdir, sim.project_name, "output", "grid3d")
        )


def set_member_outputs(sim, print_times, t_end):
    """Print the state at ``print_times`` and run until ``t_end`` in every member folder."""
    for ens_nb in range(sim.NENS):
        sim.update_parm(
            TIMPRTi=list(print_times),
            TMAX=t_end,
            filename=os.path.join(_member_dir(sim, ens_nb), "input/parm"),
            IPRT1=2,
            backup=True,
            shellprint_update=False,
        )


# --------------------------------------------------------------------------- orchestration
# --------------------------------------------------------------------------- outputs (same files as the sequential DA)
def read_member_states(sim, ens_nb, times, time_offset=0.0):
    """Pressure head and saturation of ALL nodes of one member at ``times`` (absolute): two [n_times, n_nodes] arrays."""
    out_dir = os.path.join(_member_dir(sim, ens_nb), sim.output_dirname)
    df_psi = sim.read_outputs(filename="psi", path=out_dir)
    df_sw, _ = sim.read_outputs(filename="sw", path=out_dir)
    rel = [t - time_offset for t in times]
    psi = df_psi.iloc[state_rows_at_times(df_psi, rel)].values
    sw = df_sw.iloc[state_rows_at_times(df_sw, rel)].values
    return np.asarray(psi, dtype=float), np.asarray(sw, dtype=float)


def esmda_states_store(sim, result, times, t0, atm_times, da_config):
    """``StateStore`` of the POSTERIOR forward run (the last run, made with the final parameters), one chunk per state
    time, ``cycle`` = index of the matching atmbc time (as in the sequential DA); ``analysis`` is NaN (a smoother has
    no per-cycle analysis). Members that are invalid or whose output cannot be read are flagged ``rejected``."""
    from pyCATHY.DA.state_store import StateStore

    n_ens = int(sim.NENS)
    psi_all = sw_all = None
    bad = ~np.asarray(result.valid, dtype=bool)
    for j in np.flatnonzero(~bad):
        try:
            psi, sw = read_member_states(sim, int(j), times, t0)
        except Exception as exc:  # one unreadable member must not lose the others
            logger.warning("member %d: states of the posterior run not read (%s: %s)", j + 1, type(exc).__name__, exc)
            bad[j] = True
            continue
        if psi_all is None:
            psi_all = np.full((len(times), psi.shape[1], n_ens), np.nan, dtype=np.float32)
            sw_all = np.full_like(psi_all, np.nan)
        psi_all[:, :, j], sw_all[:, :, j] = psi, sw
    if psi_all is None:
        raise RuntimeError("no member output could be read")
    store = StateStore(n_ens, nodes=getattr(da_config, "store_nodes", None), full=getattr(da_config, "store_full", "all"))
    store.attrs.update(da_type="esmda", states_source="posterior_forward_run" if len(result.prediction_history) > len(result.alphas)
                       else "last_iteration_forward_run", n_iter=int(len(result.alphas)))
    atm = np.asarray(atm_times, dtype=float)
    for k, t in enumerate(times):
        cycle = int(np.argmin(np.abs(atm - float(t)))) if atm.size else k
        store.add(cycle, psi_all[k], sw_all[k], np.full_like(psi_all[k], np.nan), rejected=bad, is_analysis=True)
    return store


def esmda_params_dataset(result, to_physical, names):
    """``xarray.Dataset``: parameters of every iterate in PHYSICAL space [iteration, member(, dim)], predictions
    [run, obs, member], misfit [run], alpha [update]. Iteration 0 = prior, last = posterior."""
    import xarray as xr

    n_ens = result.param.shape[1]
    hist = [to_physical(h) for h in result.param_history]
    data = {}
    for n in names:
        a = np.stack([np.asarray(h[n], dtype=float) for h in hist])
        safe = re.sub(r"[^0-9A-Za-z_]", "_", str(n))
        if a.ndim == 2 and a.shape[1] == n_ens:
            data[safe] = (("iteration", "member"), a)
        elif a.ndim == 3 and a.shape[2] == n_ens:
            data[safe] = (("iteration", f"{safe}_dim", "member"), a)
        else:
            logger.warning("parameter %s: shape %s not stored in _params.nc", n, a.shape)
    if result.prediction_history:
        data["prediction"] = (("run", "obs", "member"), np.stack(result.prediction_history).astype(np.float32))
        data["misfit"] = (("run",), np.asarray(result.misfit_history, dtype=float))
    data["alpha"] = (("update",), np.asarray(result.alphas, dtype=float))
    data["valid"] = (("member",), np.asarray(result.valid, dtype=bool))
    ds = xr.Dataset(data)
    ds.attrs.update(da_type="esmda", n_iter=int(len(result.alphas)),
                    description="iteration 0 = prior, last = posterior; run n_iter = posterior forward run")
    return ds


def record_esmda_outputs(sim, result, names, to_physical, times, t0, atm_times, list_assimilated_obs, da_config):
    """Everything an ES-MDA run leaves next to the project, as the sequential DA does: issue log entries,
    ``_params.nc``, ``_states.nc`` (posterior run), and the standard backup (``dict_parm_pert``, ``df_obs``, ...).
    Each part is independent and guarded: a failure here must never lose a run that took hours."""
    from pyCATHY.DA.issues import esmda_issues
    from pyCATHY.DA.run_outputs import RunOutputs
    from pyCATHY.DA.state_store import save_dataset

    def guarded(label, fn):
        try:
            fn()
        except Exception as exc:
            logger.warning("ES-MDA output '%s' not written (%s: %s)", label, type(exc).__name__, exc)

    out = RunOutputs(sim.workdir, sim.project_name)
    os.makedirs(out.root, exist_ok=True)

    def issues():
        log = getattr(sim, "issue_log", None)
        if log is not None:
            for it, kind, level, msg in esmda_issues(result.param_history, result.valid, result.misfit_history, names):
                log.record(kind, msg, cycle=it, level=level)

    def params():
        save_dataset(esmda_params_dataset(result, to_physical, names), out.params_nc)

    def states():
        if not getattr(da_config, "write_states", True):
            return
        sim._state_store = esmda_states_store(sim, result, times, t0, atm_times, da_config)
        sim.config = da_config if getattr(sim, "config", None) is None else sim.config
        sim._finalize_state_backup()

    def backup():
        meta = {"listAssimilatedObs": list_assimilated_obs, "listUpdatedparm": list(names), "da_type": "esmda",
                "alphas": [float(a) for a in result.alphas], "misfit": [float(m) for m in result.misfit_history]}
        sim.backup_results_DA(meta)

    for label, fn in (("issues", issues), ("params", params), ("states", states), ("backup", backup)):
        guarded(label, fn)


def run_esmda_on_cathy(
    sim,
    dict_obs,
    list_parm2update,
    dict_parm_pert,
    list_assimilated_obs,
    config: Optional[esmda.ESMDAConfig] = None,
    parallel: bool = True,
    da_config=None,
    dry_run: bool = False,
    **kwargs,
):
    """Run ES-MDA on the CATHY ensemble of ``sim`` (a ``DA`` object). See the module docstring.

    ``dry_run=True`` does everything up to the first model run (ensemble folders, printing
    times, observation stacking, parameter names, checks) and returns a summary ``dict`` of
    what a real run would do and cost, without simulating anything: use it to find a setup
    mistake in seconds.

    The prior ensemble is ``dict_parm_pert[name]["ini_perturbation"]``. Each iteration
    overwrites it with the current iterate (back-transformed) so that the existing writers
    (``apply_ensemble_updates``) put it in the member folders; at the end it is restored,
    and the posterior is stored under ``dict_parm_pert[name]["update_nb1"]``. The
    ``ESMDAResult`` is returned and kept as ``sim.esmda_result``.
    """
    config = esmda.ESMDAConfig() if config is None else config
    if config.update_initial_conditions:
        raise NotImplementedError("update_initial_conditions=True is not supported yet.")
    if config.cleanup_between_iterations:
        raise NotImplementedError("cleanup_between_iterations=True is not supported yet.")

    plan = observation_plan(dict_obs, list_assimilated_obs)  # validates the sensor types
    if uses_ert(plan) and not parallel:
        raise NotImplementedError("ERT observations are only mapped with parallel=True.")

    callexe, all_atmbc_times, threshold_rejected, ens_times, _ = sim.prepare_DA(
        dict_obs,
        dict_parm_pert,
        list_parm2update,
        list_assimilated_obs,
        parallel,
        config=da_config,
        **kwargs,
    )
    names = resolve_parameter_names(sim.dict_parm_pert, list_parm2update)
    prior_key, posterior_key = _update_key(0), _update_key(1)
    for n in names:
        # the layered ic writer reads dict[n]["ini_<n>_withLayers"], not the key the smoother rewrites:
        # the iterates would silently never reach the member files
        if re.split(r"(\d+)", n)[0].casefold() == "ic" and sim.dict_parm_pert[n].get("pert_control_name") == "layers":
            raise NotImplementedError(
                f"{n!r}: a layered initial condition (pert_control_name='layers') is not supported by ES-MDA "
                "yet (DA._update_ic_member reads 'ini_ic_withLayers'); use a global ic."
            )
    if uses_ert(plan):
        archie = getattr(sim, "Archie_parms", None)
        if not isinstance(archie, dict) or "porosity" not in archie:
            raise ValueError("ERT observations need the Archie parameters (DA.set_Archie_parm).")

    if uses_swc(plan):
        ensure_grid3d(sim)  # before the member parm files are written, as in selec_atmbc_window

    t0 = float(all_atmbc_times[0])
    t_end = float(all_atmbc_times[-1]) - t0
    times = state_times(plan, all_atmbc_times, config.obs_time_convention)
    print_times = sorted({0.0, *[t - t0 for t in times]})
    set_member_outputs(sim, print_times, t_end)

    obs = stack_data(sim, plan, list_assimilated_obs)
    prior_physical = {n: np.array(sim.dict_parm_pert[n][prior_key], dtype=float) for n in names}

    summary = build_run_summary(
        names=names, prior_physical=prior_physical, n_members=sim.NENS, config=config, plan=plan, times=times,
        print_times=print_times, atm_times=all_atmbc_times, n_data=obs.data.size, n_data_raw=obs.keep.size,
        use_log=obs.use_log, threshold_rejected=threshold_rejected,
        bounds_mode=getattr(sim, "param_bounds_mode", "clip"), parallel=parallel,
    )
    logger.info("\n" + "\n".join(format_run_plan(summary)))  # console + log, before the first model run
    if dry_run:
        return summary

    param0 = np.asarray(sim.transform_parameters(names, update_key=prior_key), dtype=float)

    def to_physical(param_est):
        return {
            n: np.asarray(
                sim.transform_parameters([n], param=param_est[j], back=True)[0], dtype=float
            )
            for j, n in enumerate(names)
        }

    def forward(param_est, valid):
        n_ens = param_est.shape[1]
        for n, values in to_physical(param_est).items():
            sim.dict_parm_pert[n][prior_key] = values
        sim.ens_valid = [int(i) for i in np.flatnonzero(valid)]
        sim.apply_ensemble_updates(sim.dict_parm_pert, names, cycle_nb=0)
        sim._run_ensemble_hydrological_model(parallel, callexe)
        rejected = np.asarray(sim._check_before_analysis(sim.ens_valid, threshold_rejected), dtype=bool)
        ok = valid & ~rejected

        point = {}
        for j in np.flatnonzero(ok):
            try:
                point[int(j)] = predict_member_blocks(sim, int(j), plan, times, time_offset=t0)
            except ESMDATimeError:
                raise
            except _MEMBER_READ_ERRORS as exc:
                logger.warning("member %d: cannot read the output (%s: %s)", j + 1, type(exc).__name__, exc)
                ok[j] = False
        ok_members = [int(j) for j in np.flatnonzero(ok)]
        if not ok_members:
            raise RuntimeError("ES-MDA: no member produced a usable prediction.")
        sim.ens_valid = ok_members  # as _run_window: the ERT forward runs on the valid members

        ert = predict_ert_blocks(sim, plan, times, ok_members, time_offset=t0)

        rows = []
        for i, key, n in obs.layout:
            block = np.full((n, n_ens), np.nan)
            if (i, key) in ert:
                got = ert[(i, key)]
                if got.shape[0] != n:
                    raise ValueError(
                        f"{key} at t={plan[i][0]}: the forward gives {got.shape[0]} data, the "
                        f"observations hold {n} (different quadrupole sequence?)."
                    )
                block[:, ok_members] = got
            else:
                for j in ok_members:
                    col = point[j][(i, key)]
                    if col.size != n:
                        raise ValueError(
                            f"{key} at t={plan[i][0]}: member {j + 1} predicts {col.size} data, "
                            f"the observations hold {n}."
                        )
                    block[:, j] = col
            rows.append(block)
        return obs.transform_prediction(np.vstack(rows)), ok

    def project(param_est):
        """Bounds of the scenario (``clip_min`` / ``clip_max``) on the physical values."""
        physical = to_physical(param_est)
        out = param_est.copy()
        for j, n in enumerate(names):
            clipped = sim._clip_updated_parameter(n, physical[n])
            out[j] = np.asarray(sim.transform_parameters([n], param=clipped)[0], dtype=float)
        return out

    try:
        result = esmda.run_esmda(param0, forward, obs.data, obs.cov, config, project=project)
    finally:
        # keep dict_parm_pert meaningful whatever happens: prior under its own key
        for n in names:
            sim.dict_parm_pert[n][prior_key] = prior_physical[n]

    for n, values in to_physical(result.param).items():
        sim.dict_parm_pert[n][posterior_key] = values
    sim.ens_valid = [int(i) for i in np.flatnonzero(result.valid)]
    result.param_names = list(names)
    sim.esmda_result = result
    record_esmda_outputs(sim, result, names, to_physical, times, t0, all_atmbc_times, list_assimilated_obs,
                         da_config if da_config is not None else getattr(sim, "config", None))
    return result
