"""Light issue log of a DA run: what went wrong, once per distinct problem, not once per print.

``console_log.txt`` duplicates the whole console (tables, per-cycle prints) and grows without bound. The issue log
keeps only WARNING+ messages of the ``pyCATHY`` logger and the structured events recorded with ``IssueLog.record``.
Repeated messages (same text once the numbers are masked) are counted instead of repeated, so a 400-cycle run
that warns at every cycle still produces one line.

Line format (tab separated, greppable)::

    LEVEL  kind  count  cycles(first-last)  message

``IssueLog.summary()`` is a small dict for the run manifest.
"""
from __future__ import annotations

import logging
import os
import re
import threading

import numpy as np

_NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

# kinds of structured events (free text is allowed too: 'log' is used for captured logger messages)
KINDS = ("rejection", "collapse", "degeneracy", "bounds", "rank", "obs_perturbation", "ert", "clipping", "misfit", "run", "log")

# ensemble spread after / before the analysis (mean over nodes of the std over members)
SPREAD_COLLAPSE = 0.1      # below: the ensemble has (almost) collapsed -> ERROR level
SPREAD_WARN = 0.3          # below: strong loss of spread -> WARNING
MISFIT_WARN = 2.0          # ES-MDA: reduced chi-square of the last forward run above this = data not fitted
NEFF_FRACTION = 0.1        # particle filter: N_eff below this share of the members = degeneracy


class IssueLog:
    def __init__(self, path=None, max_distinct=200, cycle_getter=None):
        self.path = path
        self.cycle_getter = cycle_getter
        self.max_distinct = int(max_distinct)
        self._items = {}                       # key -> dict
        self._dropped = 0
        self._lock = threading.Lock()
        self._handler = None
        if path:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            open(path, "w").close()            # one log per run

    # ------------------------------------------------------------------ recording
    def record(self, kind, message, cycle=None, level="WARNING", **info):
        if cycle is None and self.cycle_getter is not None:
            try:
                cycle = self.cycle_getter()
            except Exception:
                cycle = None
        key = (kind, _NUM.sub("#", str(message)))
        with self._lock:
            it = self._items.get(key)
            if it is None:
                if len(self._items) >= self.max_distinct:
                    self._dropped += 1
                    return
                it = self._items[key] = dict(level=level, kind=kind, message=str(message), count=0,
                                             first=cycle, last=cycle, info=info)
            it["count"] += 1
            if cycle is not None:
                it["first"] = cycle if it["first"] is None else min(it["first"], cycle)
                it["last"] = cycle if it["last"] is None else max(it["last"], cycle)
            if level == "ERROR":
                it["level"] = "ERROR"
            self._flush()

    # ------------------------------------------------------------------ reading
    def summary(self):
        with self._lock:
            by_kind = {}
            for it in self._items.values():
                by_kind[it["kind"]] = by_kind.get(it["kind"], 0) + it["count"]
            return dict(n_distinct=len(self._items), n_events=sum(by_kind.values()), by_kind=by_kind,
                        n_errors=sum(1 for i in self._items.values() if i["level"] == "ERROR"),
                        dropped=self._dropped, file=self.path)

    def lines(self):
        out = []
        for it in sorted(self._items.values(), key=lambda i: (i["level"] != "ERROR", i["first"] is None, i["first"])):
            cyc = "-" if it["first"] is None else (str(it["first"]) if it["first"] == it["last"]
                                                  else f"{it['first']}-{it['last']}")
            msg = " ".join(it["message"].split())
            out.append(f"{it['level']}\t{it['kind']}\tx{it['count']}\t{cyc}\t{msg[:300]}")
        if self._dropped:
            out.append(f"WARNING\tlog\tx{self._dropped}\t-\t{self._dropped} further distinct issues not listed")
        return out

    def _flush(self):                          # tiny file, rewritten: a crash keeps the issues found so far
        if self.path:
            try:
                with open(self.path, "w", encoding="utf-8") as f:
                    f.write("\n".join(self.lines()) + "\n")
            except OSError:
                pass

    # ------------------------------------------------------------------ logger capture
    def attach(self, logger_name="pyCATHY", level=logging.WARNING):
        """Capture WARNING+ of the package logger (the 46 existing ``logger.warning`` calls) with no code change."""
        log = self

        class _H(logging.Handler):
            def emit(self, record):
                try:
                    log.record("log", record.getMessage(), level=record.levelname, logger=record.name)
                except Exception:              # diagnostics must never stop a run
                    pass

        self._handler = _H(level=level)
        logging.getLogger(logger_name).addHandler(self._handler)
        return self

    def detach(self, logger_name="pyCATHY"):
        if self._handler is not None:
            logging.getLogger(logger_name).removeHandler(self._handler)
            self._handler = None


# ---------------------------------------------------------------------- detectors (pure functions, no CATHY needed)
def spread(ensemble):
    """Mean over rows of the std over members of ``ensemble`` [n, n_members] (NaN members ignored)."""
    a = np.asarray(ensemble, dtype=float)
    if a.ndim != 2 or a.shape[1] < 2 or a.shape[0] == 0:
        return float("nan")
    return float(np.nanmean(np.nanstd(a, axis=1, ddof=1)))


def spread_ratio(before, after):
    """sd(after) / sd(before): 1 = unchanged, -> 0 = the ensemble collapsed (members become copies of each other)."""
    sb, sa = spread(before), spread(after)
    if not np.isfinite(sb) or not np.isfinite(sa):
        return float("nan")
    return float("inf") if sb == 0 and sa > 0 else (1.0 if sb == 0 else sa / sb)


def _members_last(a, n_members):
    """Orient a parameter array as [n_param, n_members] whatever the filter returned."""
    a = np.atleast_2d(np.asarray(a, dtype=float))
    if a.shape[1] != n_members and a.shape[0] == n_members:
        a = a.T
    return a


def analysis_issues(state_before=None, state_after=None, param_before=None, param_after=None, param_names=None,
                    pf_extra=None, n_members=None):
    """Problems visible in one analysis step: list of ``(kind, level, message)``.

    state_*: [n_nodes, n_members]; param_*: [n_param, n_members] or [n_members, n_param];
    pf_extra: ``AnalysisResult.extra`` of the particle filter (``n_eff``, ``resampled``).
    """
    out = []
    if state_before is not None and state_after is not None and np.shape(state_before) == np.shape(state_after):
        r = spread_ratio(state_before, state_after)
        if np.isfinite(r) and r < SPREAD_WARN:
            lvl = "ERROR" if r < SPREAD_COLLAPSE else "WARNING"
            out.append(("collapse", lvl, f"state ensemble spread x{r:.2f} after the analysis (members becoming identical)"))
    if param_before is not None and param_after is not None:
        n = int(n_members or np.atleast_2d(param_before).shape[-1])
        pb, pa = _members_last(param_before, n), _members_last(param_after, n)
        if pb.shape == pa.shape and pb.size:
            for j in range(pb.shape[0]):
                r = spread_ratio(pb[j:j + 1], pa[j:j + 1])
                if np.isfinite(r) and r < SPREAD_WARN:
                    name = param_names[j] if param_names is not None and j < len(param_names) else f"param{j}"
                    lvl = "ERROR" if r < SPREAD_COLLAPSE else "WARNING"
                    out.append(("collapse", lvl, f"parameter {name}: spread x{r:.2f} after the analysis"))
        if n and pb.shape[0] >= n:
            out.append(("rank", "WARNING", f"{pb.shape[0]} parameters >= {n} members: rank-deficient update"))
    if pf_extra and n_members:
        n_eff = pf_extra.get("n_eff")
        if n_eff is not None and np.isfinite(n_eff) and n_eff < NEFF_FRACTION * n_members:
            out.append(("degeneracy", "ERROR" if n_eff < 2 else "WARNING",
                        f"particle filter N_eff = {float(n_eff):.1f} of {n_members} members (degeneracy)"))
    return out


def esmda_issues(param_history, valid, misfit_history=None, names=None):
    """Problems visible in an ES-MDA run: list of ``(iteration, kind, level, message)``.

    ``param_history``: list of [n_param, n_members] (prior, then after every update); the spread change is checked
    between consecutive iterates on the valid members; ``misfit_history``: reduced chi-square of every forward run.
    """
    out = []
    valid = np.asarray(valid, dtype=bool)
    if valid.sum() >= 2:
        for k in range(len(param_history) - 1):
            b, a = np.asarray(param_history[k])[:, valid], np.asarray(param_history[k + 1])[:, valid]
            for kind, level, msg in analysis_issues(param_before=b, param_after=a, param_names=names,
                                                    n_members=int(valid.sum())):
                out.append((k + 1, kind, level, f"iteration {k + 1}: {msg}"))
    m = np.asarray(misfit_history if misfit_history is not None else [], dtype=float)
    if m.size and np.isfinite(m[-1]) and m[-1] > MISFIT_WARN:
        out.append((int(m.size), "misfit", "WARNING",
                    f"final reduced chi-square {m[-1]:.2f} > {MISFIT_WARN:g}: the data are not fitted"))
    if m.size > 1 and np.isfinite(m).all() and m[-1] > m[0]:
        out.append((int(m.size), "misfit", "WARNING", f"misfit increased from {m[0]:.2f} to {m[-1]:.2f} over the iterations"))
    if valid.size and (~valid).any():
        out.append((None, "rejection", "WARNING",
                    f"{int((~valid).sum())} of {valid.size} members dropped during the run: {[int(i) + 1 for i in np.flatnonzero(~valid)]}"))
    return out
