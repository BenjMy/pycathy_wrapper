"""Check how observation times line up with the model windows, before a run (B4).

How ``run_DA_sequential`` pairs observations and model states (read from the code):

* the model is run window by window over ``[t, t_next]``, for ``t`` in ``all_atmbc_times[:-1]``;
* after the window starting at ``t`` an analysis is made **if ``t`` is in ``ENS_times``**
  (exact comparison of floats);
* the data of the k-th analysis are the k-th entry of ``dict_obs`` (selected by counter,
  not by time);
* the predicted observation is computed from the **last time of the member output**, i.e.
  from the state at ``t_next``, the end of the window.

Consequences, reported by ``check_time_alignment``: an observation stamped ``t`` is compared
with the state at ``t_next``; an observation whose time is not exactly a window start is
never assimilated **and every later observation is then paired with the wrong analysis**;
the observation at the final time (TMAX) is never assimilated.
"""
from dataclasses import dataclass
from typing import List, Sequence

import numpy as np


@dataclass(frozen=True)
class Finding:
    level: str  # "error" (wrong pairing), "warning" (suspicious), "info"
    message: str


def check_time_alignment(ENS_times: Sequence[float], all_atmbc_times: Sequence[float]) -> List[Finding]:
    """Findings on the alignment of observation times with the model windows."""
    obs = np.asarray(list(ENS_times), dtype=float)
    atmbc = np.asarray(list(all_atmbc_times), dtype=float)
    findings: List[Finding] = []
    if obs.size == 0 or atmbc.size < 2:
        return findings

    starts = atmbc[:-1]  # windows start at every atmbc time except the last
    is_start = np.isin(obs, starts)  # exact comparison, as in the run loop
    final = obs == atmbc[-1]

    if np.any(np.diff(obs) < 0):
        findings.append(Finding("error", "ENS_times is not sorted: observations are paired with analyses by order"))
    if len(set(obs.tolist())) != obs.size:
        findings.append(Finding("warning", "ENS_times contains duplicated times"))

    lost = obs[~is_start]
    if lost.size:
        near = []
        for t in lost:
            scale = max(abs(t), 1.0)
            if np.any(np.abs(starts - t) <= 1e-6 * scale) and not final[list(obs).index(t)]:
                near.append(t)
        n_final = int(final.sum())
        reasons = []
        if n_final:
            reasons.append(f"{n_final} at the final time {atmbc[-1]:g} s (no window starts there)")
        if near:
            reasons.append(f"{len(near)} almost equal to a window start but not identical (floating point)")
        other = lost.size - n_final - len(near)
        if other:
            reasons.append(f"{other} not at any window start")
        findings.append(Finding(
            "error",
            f"{lost.size} of {obs.size} observation time(s) will never be assimilated "
            f"({'; '.join(reasons)}): {np.array2string(lost[:6], precision=6)}"
            f"{' ...' if lost.size > 6 else ''}. The k-th analysis uses the k-th observation, so every "
            "observation after the first skipped one is paired with the wrong analysis."))

    used = obs[is_start]
    if used.size:
        next_time = {s: e for s, e in zip(atmbc[:-1], atmbc[1:])}
        dt = np.array([next_time[t] - t for t in used])
        findings.append(Finding(
            "warning",
            f"{used.size} observation(s) are assimilated after the window that STARTS at their time: an "
            f"observation stamped t is compared with the model state at t + {dt.min():g} s"
            f"{'' if dt.min() == dt.max() else f' (up to t + {dt.max():g} s)'}, the last time of the "
            "member output. Check that the data were generated with the same convention."))
    return findings
