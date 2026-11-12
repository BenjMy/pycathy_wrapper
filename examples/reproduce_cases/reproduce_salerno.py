#!/usr/bin/env python3
"""
Salerno case (catchment with greenhouses and the irrigation module, ~611 d) through the pyCATHY wrapper,
on a TRUNCATED window (first N days), with the case's OWN source (`src/`, irrigation source of Sep 2024).

Stages (each can be run alone with --stage):
  1 stage   : trimmed copy of the case (sources + inputs, no outputs)
  2 check   : parm layout vs the wrapper table, irrigation_map header/size/nominal volume of the window
  3 atmbc   : input/atmbc is NOT in the repository (.gitignore, too heavy): give it with --atmbc FILE
              (--smoke writes a zero-flux homogeneous atmbc to test the mechanics ONLY: results are not comparable)
  4 run     : wrapper run of the window (mesh first, then IPRT1=2), original time stepping
  5 compare : hgatmsf / hgraph / hgirrigation of the run vs the reference outputs of the repository

usage:
    python reproduce_salerno.py --days 3 --atmbc /path/to/atmbc            # real comparison
    python reproduce_salerno.py --days 0.25 --smoke                         # mechanics only
Outputs: <work>/salerno_summary.json and <work>/salerno_compare.png
"""
import argparse
import json
import os
import shutil
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
from case_tools import case_dir, get_repo, salerno_irrigation_from_tables, stage_case
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.irrigation import read_irrigation_map
from pyCATHY.version_config import CATHY_DIR_LAYOUT

DAY = 86400.0


def last_per_step(arr):
    """a rejected step is written again: keep the last line of each step (first column = step number)"""
    last = {}
    for i, k in enumerate(arr[:, 0].astype(int)):
        last[k] = i
    return arr[sorted(last.values())]


def load_numeric(path, ncols=None):
    rows = []
    with open(path, errors="replace") as f:
        for line in f:
            if line[:1] == "#":
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                pass
    n = ncols or max(set(map(len, rows)), key=list(map(len, rows)).count)
    return np.array([r for r in rows if len(r) == n])


# ----------------------------------------------------------------------------- stages
def stage_check(proj, days):
    info = {}
    parm = in_CT.read_parm(os.path.join(proj, "input", "parm"))
    fmt = CATHY_DIR_LAYOUT["withIrr"]["parm_fmt"]
    info["parm_entries"], info["parm_fmt_sum"] = len(parm), int(sum(fmt))
    info["tmax_days"] = float(parm["TMAX"]) / DAY
    irr_path = os.path.join(proj, "other-dem", "irrigation_map")
    dem = np.loadtxt(os.path.join(proj, "prepro-dem", "dem"), skiprows=6)
    m = read_irrigation_map(irr_path, dem.shape[0], dem.shape[1])
    info.update(irr_flag=m["irr_flag"], irr_sw=m["irr_sw"], irr_hs=m["irr_hs"], irr_blocks=int(len(m["times"])))
    cell = float(parm.get("DELTA_X", 80.0)) if False else 80.0
    edges = np.append(m["times"], max(m["times"].max(), parm["TMAX"]))
    w = np.minimum(edges[1:], days * DAY) - np.minimum(edges[:-1], days * DAY)
    info["irr_nominal_window_m3"] = float(np.sum(m["rates"].sum(axis=(1, 2)) * cell * cell * w))
    print("check:", json.dumps(info, indent=1))
    assert info["parm_entries"] == info["parm_fmt_sum"], "parm layout differs from the wrapper table"
    return info


def stage_atmbc(proj, simu, atmbc, smoke, t_end):
    dst = simu._input_path("atmbc")
    if atmbc:
        shutil.copy(atmbc, dst)
        print("atmbc copied:", dst, round(os.path.getsize(dst) / 1e6, 1), "MB")
    elif smoke:
        simu.update_atmbc(HSPATM=1, IETO=1, time=[0.0, t_end], VALUE=[np.zeros(2), np.zeros(2)])
        print("SMOKE atmbc written (zero flux): results are NOT comparable with the reference")
    elif not os.path.isfile(dst):
        raise SystemExit("input/atmbc is not in the repository (.gitignore). Give it with --atmbc FILE "
                         "(or --smoke to test the mechanics only).")


def stage_run(proj, atmbc, smoke, days, regen_irrigation=False, src_case=None):
    simu = cathy_tools.CATHY(dirName=proj + "/", prj_name="", version="withIrr")
    t_end = days * DAY
    stage_atmbc(proj, simu, atmbc, smoke, t_end)
    if regen_irrigation:       # irrigation_map rebuilt with the wrapper tools from the tables of the repository
        times, rates = salerno_irrigation_from_tables(src_case)
        simu.update_irrigation_map(times=times, rates=rates, irr_sw=0.35, irr_flag=1)
        print("irrigation_map regenerated:", len(times), "blocks")
    # vegetation classes of the case (bounds of CATHY.H are never reduced below the case's own)
    veg = np.loadtxt(os.path.join(proj, "other-dem", "veg-type_map"), skiprows=6)
    simu.update_cathyH(MAXVEG=int(len(np.unique(veg))))
    # daily prints inside the window (the original print times are far later: hg* files are per step)
    t_print = list(np.arange(1.0, np.floor(days) + 1.0) * DAY) if days >= 1 else [t_end / 2, t_end]
    simu.run_processor(IPRT1=2, TRAFLAG=0, TMAX=t_end, TIMPRTi=t_print, verbose=False)
    return simu


def stage_compare(src_case, proj, work, days):
    ref_dir, new_dir = os.path.join(src_case, "output"), os.path.join(proj, "output")
    out = {}
    fig, ax = plt.subplots(2, 2, figsize=(14, 8))
    # hgatmsf: surface fluxes per step -> cumulative volumes
    R = last_per_step(load_numeric(os.path.join(ref_dir, "hgatmsf"), 10))
    N = last_per_step(load_numeric(os.path.join(new_dir, "hgatmsf"), 10))
    for name, col in (("POT", 3), ("ACT", 4), ("OVL", 5), ("SEEP", 7), ("REC", 8)):
        for tag, A in (("ref", R), ("new", N)):
            s = A[:, 2] <= days * DAY
            out[f"cum_{name}_{tag}_m3"] = float(np.sum(A[s, col] * A[s, 1]))
    t_r, t_n = R[:, 2] / DAY, N[:, 2] / DAY
    ax[0, 0].plot(t_r, R[:, 9], label="reference"); ax[0, 0].plot(t_n, N[:, 9], "--", label="wrapper")
    ax[0, 0].set_xlim(0, days); ax[0, 0].set_title("REC.VOL (hgatmsf)"); ax[0, 0].legend(); ax[0, 0].set_xlabel("days")
    for tag, A, ls in (("ref", R, "-"), ("new", N, "--")):
        ax[0, 1].plot(A[:, 2] / DAY, np.cumsum(A[:, 4] * A[:, 1]), ls, label=tag)
    ax[0, 1].set_xlim(0, days); ax[0, 1].set_title("cumulative actual surface flux (ACT) [m3]"); ax[0, 1].set_xlabel("days")
    # hgraph: time + outlet columns
    try:
        Hr = load_numeric(os.path.join(ref_dir, "hgraph"), 5); Hn = load_numeric(os.path.join(new_dir, "hgraph"), 5)
        sel = Hr[:, 0] <= days * DAY
        for j in range(1, 5):
            ax[1, 0].plot(Hr[:, 0] / DAY, Hr[:, j], label=f"ref col {j}"); ax[1, 0].plot(Hn[:, 0] / DAY, Hn[:, j], "--")
            out[f"hgraph_col{j}_maxdiff"] = float(np.max(np.abs(np.interp(Hr[sel, 0], Hn[:, 0], Hn[:, j]) - Hr[sel, j])))
        ax[1, 0].set_xlim(0, days); ax[1, 0].set_title("hgraph (outlet series), reference vs wrapper")
    except Exception as exc:
        out["hgraph_error"] = str(exc)
    # hgirrigation: dedupe steps, cumulative available / unused
    try:
        Gr = last_per_step(load_numeric(os.path.join(ref_dir, "hgirrigation")))
        Gn = last_per_step(load_numeric(os.path.join(new_dir, "hgirrigation")))
        for tag, G, ls in (("ref", Gr, "-"), ("new", Gn, "--")):
            s = G[:, 2] <= days * DAY
            out[f"irr_available_{tag}_m3"] = float(np.sum(G[s, 4] * G[s, 1]))
            out[f"irr_unused_{tag}_m3"] = float(np.sum(G[s, 6] * G[s, 1]))
            ax[1, 1].plot(G[:, 2] / DAY, np.cumsum(G[:, 4] * G[:, 1]), ls, label="available " + tag)
        ax[1, 1].set_xlim(0, days); ax[1, 1].set_title("irrigation available (hgirrigation) [m3]"); ax[1, 1].legend()
    except Exception as exc:
        out["hgirrigation_error"] = str(exc)
    fig.tight_layout(); fig.savefig(os.path.join(work, "salerno_compare.png"), dpi=110)
    print("compare:", json.dumps(out, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo"); ap.add_argument("--work", default="work_salerno")
    ap.add_argument("--days", type=float, default=3.0)
    ap.add_argument("--atmbc", help="path of the original (heavy) input/atmbc of the case")
    ap.add_argument("--smoke", action="store_true", help="zero-flux atmbc: mechanics only")
    ap.add_argument("--regen-irrigation", action="store_true", help="rebuild other-dem/irrigation_map from the tables of the repository")
    ap.add_argument("--stage", choices=["stage", "check", "run", "compare", "all"], default="all")
    a = ap.parse_args()
    repo = get_repo(a.repo)
    src_case = case_dir(repo, "Salerno")
    work = os.path.abspath(a.work); os.makedirs(work, exist_ok=True)
    proj = os.path.join(work, os.path.basename(src_case))
    summary = {}
    if a.stage in ("stage", "all"):
        proj = stage_case(src_case, work)
    if a.stage in ("check", "all"):
        summary["check"] = stage_check(proj, a.days)
    if a.stage in ("run", "all"):
        stage_run(proj, a.atmbc, a.smoke, a.days, a.regen_irrigation, src_case)
    if a.stage in ("compare", "all") and not a.smoke:
        summary["compare"] = stage_compare(src_case, proj, work, a.days)
    json.dump(summary, open(os.path.join(work, "salerno_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
