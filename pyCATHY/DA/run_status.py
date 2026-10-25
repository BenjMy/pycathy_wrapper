"""What did each DA run leave on disk? One line per run, whatever the DA type (EnKF, PF, ES-MDA).

    python -m pyCATHY.DA.run_status ../simulations --glob "SYN_NC_*"          # table
    python -m pyCATHY.DA.run_status ../simulations --glob "SYN_NC_*" --strict  # exit 1 if a run is not plottable

Verdicts: ``ok`` (states + parameters/obs backup) | ``partial`` (some outputs missing) | ``no states`` |
``no results`` (folder without outputs: not started, still running or crashed) | ``failed`` (the run's own json says so).
"""
from __future__ import annotations

import argparse
import fnmatch
import glob
import json
import os
import sys

BACKUP_FILES = ("dict_parm_pert", "df_obs", "df_performance", "ET_DA_xr", "df_Archie", "meta_DA")
REQUIRED = ("dict_parm_pert", "df_obs")                 # what a plot needs besides the states
_SKIP_JSON = ("_rejections.json", "_ERT_checks.json")


def _has(base, name):
    return any(os.path.exists(f"{base}_{name}{ext}") for ext in (".pkl", ".parquet"))


def _states_kind(base):
    if os.path.exists(base + "_states.nc"):
        return "nc"
    d = base + "_states"
    if os.path.isdir(d):
        n = len(glob.glob(os.path.join(d, "cycle_*.npz")))
        return f"chunks({n})" if n else "empty folder"
    return "none"


def _run_status(project_dir):
    """``status`` written by the run script in its own json (any *.json of the folder that has the key)."""
    for path in sorted(glob.glob(os.path.join(project_dir, "*.json"))):
        if path.endswith(_SKIP_JSON):
            continue
        try:
            with open(path, encoding="utf-8") as f:
                d = json.load(f)
            if isinstance(d, dict) and "status" in d:
                return str(d["status"])
        except (OSError, ValueError):
            continue
    name = os.path.basename(os.path.normpath(project_dir))      # fallback: the manifest written by the library itself
    try:
        with open(os.path.join(project_dir, f"{name}_run.json"), encoding="utf-8") as f:
            d = json.load(f)
        if isinstance(d, dict) and "library_status" in d:
            return str(d["library_status"])
    except (OSError, ValueError):
        pass
    return None


def _issue_counts(base):
    path = base + "_issues.log"
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]
    return {"errors": sum(ln.startswith("ERROR") for ln in lines), "warnings": sum(ln.startswith("WARNING") for ln in lines)}


def inspect_run(project_dir):
    """Dict describing one run folder (see module docstring for ``verdict``)."""
    project_dir = os.path.normpath(project_dir)
    name = os.path.basename(project_dir)
    base = os.path.join(project_dir, name)
    present = [n for n in BACKUP_FILES if _has(base, n)]
    states = _states_kind(base)
    params = os.path.exists(base + "_params.nc")
    status = _run_status(project_dir)
    results = bool(present) or states in ("nc",) or states.startswith("chunks")
    if status and status.lower() in ("failed", "error", "aborted", "crashed"):
        verdict = "failed"
    elif not results:
        verdict = "no results" + (" (status: finished!)" if status and status.lower() == "finished" else "")
    elif states in ("none", "empty folder"):
        verdict = "no states"
    elif all(r in present for r in REQUIRED):
        verdict = "ok"
    else:
        verdict = "partial"
    return dict(run=name, verdict=verdict, status=status, states=states, params=params,
                missing=[n for n in BACKUP_FILES if n not in present], issues=_issue_counts(base))


def inspect_runs(simulations_dir, pattern="*"):
    dirs = [d for d in sorted(glob.glob(os.path.join(simulations_dir, "*"))) if os.path.isdir(d)
            and fnmatch.fnmatch(os.path.basename(d), pattern)]
    return [inspect_run(d) for d in dirs]


def format_table(rows):
    if not rows:
        return "no run found"
    w = max(len(r["run"]) for r in rows)
    out = [f"{'run':<{w}}  {'verdict':<12} {'status':<9} {'states':<13} params  issues(E/W)  missing"]
    for r in rows:
        iss = "-" if r["issues"] is None else f"{r['issues']['errors']}/{r['issues']['warnings']}"
        miss = ",".join(m for m in r["missing"] if m != "meta_DA") or "-"
        out.append(f"{r['run']:<{w}}  {r['verdict']:<12} {str(r['status'] or '-'):<9} {r['states']:<13} "
                   f"{'yes' if r['params'] else '-':<6}  {iss:<12} {miss}")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("simulations_dir")
    ap.add_argument("--glob", default="*", help="run folder pattern, e.g. 'SYN_NC_*'")
    ap.add_argument("--strict", action="store_true", help="exit 1 if a run is not ok")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    rows = inspect_runs(a.simulations_dir, a.glob)
    print(json.dumps(rows, indent=1) if a.json else format_table(rows))
    return 1 if a.strict and any(r["verdict"] != "ok" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
