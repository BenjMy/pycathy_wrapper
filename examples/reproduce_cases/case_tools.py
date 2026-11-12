"""
Helpers to reproduce the study cases of CATHY_src_WithIrrigation (examples/Venezia..., examples/Salerno...)
through the pyCATHY wrapper.

  * get_repo()    : local clone of the source repo (shallow clone from GitHub if absent)
  * stage_case()  : trimmed copy of a case (sources + inputs, no outputs) to run in
  * compare_sw()  : numerical comparison of a nodal output file (sw / psi) with the reference run
"""
import os
import shutil
import subprocess

import numpy as np

REPO_URL = "https://github.com/CATHY-Org/CATHY_src_WithIrrigation.git"
CASE_DIRS = ("src", "input", "other-dem", "prepro-dem", "prepro", "R", "matlab")
CASE_FILES = ("cathy.fnames", "CATHY.H")


def get_repo(path=None):
    """Path of a local clone of CATHY_src_WithIrrigation (env CATHY_SRC_REPO or ./CATHY_src_WithIrrigation)."""
    path = os.path.abspath(path or os.environ.get("CATHY_SRC_REPO", "CATHY_src_WithIrrigation"))
    if not os.path.isdir(os.path.join(path, "examples")):
        subprocess.run(["git", "clone", "--depth", "1", REPO_URL, path], check=True)
    return path


def case_dir(repo, prefix):
    names = [d for d in os.listdir(os.path.join(repo, "examples")) if d.startswith(prefix) and not d.endswith(".py")]
    if not names:
        raise FileNotFoundError(f"no example starting with {prefix} in {repo}/examples")
    return os.path.join(repo, "examples", names[0])


def stage_case(src_case, work_root, subdirs=("src", "input", "other-dem", "prepro-dem", "prepro"), clean=True):
    """Copy sources and inputs of a case (never its outputs) to work_root/<case name>; returns the path."""
    dst = os.path.join(os.path.abspath(work_root), os.path.basename(src_case))
    if clean:
        shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst, exist_ok=True)
    for d in subdirs:
        if os.path.isdir(os.path.join(src_case, d)):
            shutil.copytree(os.path.join(src_case, d), os.path.join(dst, d), dirs_exist_ok=True)
    for f in CASE_FILES:
        if os.path.isfile(os.path.join(src_case, f)):
            shutil.copy(os.path.join(src_case, f), dst)
    for d in ("output", "vtk"):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    return dst


def compare_nodal(ref_path, new_path, reader, t_max_days=None):
    """
    Compare two nodal output files (same mesh) with ``reader`` (e.g. pyCATHY ``read_sw``).
    The new run is interpolated in time on the reference print times (the solver prints at the first
    step after each requested time, so the times differ slightly).
    Returns a dict of summary numbers and the aligned arrays (times in days).
    """
    ref, new = reader(ref_path), reader(new_path)
    ref = ref[0] if isinstance(ref, (tuple, list)) else ref       # read_sw returns (frame, ...) in some versions
    new = new[0] if isinstance(new, (tuple, list)) else new
    t_ref = ref.index.values.astype(float) / 86400.0
    t_new = new.index.values.astype(float) / 86400.0
    t_hi = min(t_ref.max(), t_new.max())
    if t_max_days is not None:
        t_hi = min(t_hi, t_max_days)
    sel = (t_ref <= t_hi) & (t_ref >= t_new.min())
    R = ref.values[sel]
    N = np.empty_like(R)
    for j in range(R.shape[1]):
        N[:, j] = np.interp(t_ref[sel], t_new, new.values[:, j])
    d = N - R
    out = {
        "n_times": int(sel.sum()), "t_window_days": float(t_hi),
        "max_abs_diff": float(np.abs(d).max()), "mean_abs_diff": float(np.abs(d).mean()),
        "rmse": float(np.sqrt((d ** 2).mean())), "ref_range": (float(R.min()), float(R.max())),
    }
    return out, t_ref[sel], R, N


def salerno_irrigation_from_tables(case):
    """
    (times, rates) of the Salerno irrigation_map rebuilt from the tables of the repository:
    R/Table_NEUMANN.csv (hourly irrigation), other-dem/veg-type_map (land use, classes 3 and 4 irrigated)
    and OPAC/map_kc_YYMM.csv (monthly crop coefficient). Identical to the shipped file (validated).
    """
    import pandas as pd
    from pyCATHY.irrigation import irrigation_map_from_series, select_irrigation_records

    tab = pd.read_csv(os.path.join(case, "R", "Table_NEUMANN.csv"), sep=";")
    land = np.loadtxt(os.path.join(case, "other-dem", "veg-type_map"), skiprows=6)
    kc = {f[7:11]: pd.read_csv(os.path.join(case, "OPAC", f), sep=";", header=None).values.astype(float)
          for f in sorted(os.listdir(os.path.join(case, "OPAC"))) if f.startswith("map_kc_") and f.endswith(".csv")}
    sub = tab.iloc[select_irrigation_records(tab["Irr_m3s"].values)]
    times = sub["Date_pos"].values - tab["Date_pos"].values[0]
    keys = [d[8:10] + d[3:5] for d in sub["Date"]]
    return irrigation_map_from_series(times, sub["Irr_m3s"].values, keys, land, kc)
