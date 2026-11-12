#!/usr/bin/env python3
"""
Venezia case (soil column, 30 d, SWC probes) through the pyCATHY wrapper.

Level A (this script): the case's OWN old source (`src/`, identical to dataset_oldCATHY/src) and its
original inputs are run with the wrapper (default version, layout input/ + prepro/) and the nodal
saturation `sw` is compared with the reference run shipped in the repository.

usage:
    python reproduce_venezia.py [--repo PATH_TO_CATHY_src_WithIrrigation] [--work work] [--days 3]
Writes <work>/venezia_levelA_summary.json and <work>/venezia_levelA_sw.png.

Notes
  * the old source does not write output/grid3d, hence run_processor(read_grid=False);
  * the original parm has a hand edited label (TMAX608400.086400.0): handled by the wrapper reader;
  * Level C (same inputs translated to the new layout and run with the new source) is not in this file yet.
"""
import argparse
import json
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
from case_tools import case_dir, compare_nodal, get_repo, stage_case
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo"); ap.add_argument("--work", default="work_venezia")
    ap.add_argument("--reuse", action="store_true", help="do not stage/run again if <work> already holds a finished run")
    ap.add_argument("--days", type=float, default=8.0, help="simulated days (the full case is 30; the first print of the case is at day 6.08)")
    a = ap.parse_args()

    repo = get_repo(a.repo)
    src_case = case_dir(repo, "Venezia")
    work = os.path.abspath(a.work); os.makedirs(work, exist_ok=True)
    existing = os.path.join(work, os.path.basename(src_case))
    reuse = a.reuse and os.path.isfile(os.path.join(existing, "output", "sw")) and os.path.getsize(os.path.join(existing, "output", "sw")) > 10000
    proj = existing if reuse else stage_case(src_case, work)

    # layers (depth of the node layers) from the case's own dem_parameters
    dem_par = in_CT.read_dem_parameters(os.path.join(proj, "input", "dem_parameters"))
    zkey = next(k for k in dem_par if k.startswith("zratio"))     # label spelling differs between files
    z = np.append(0.0, np.cumsum(dem_par[zkey]) * dem_par["base"])
    nrow_nodes, ncol_nodes = 11, 14                       # 10 x 13 cells -> 11 x 14 surface nodes

    t_end = a.days * 86400.0
    if not reuse:
        run_case(proj, t_end)
    compare_and_plot(a, src_case, proj, work, z, nrow_nodes, ncol_nodes)


def run_case(proj, t_end):
    simu = cathy_tools.CATHY(dirName=proj + "/", prj_name="", version="1.0.0")

    # original print times of the case, truncated to the window (TMAX alone is reset to the last print time)
    parm = in_CT.read_parm(os.path.join(proj, "input", "parm"))
    t_prt = [t for t in np.atleast_1d(parm["(TIMPRT(I),I=1,NPRT)"]) if t <= t_end]
    # number of vegetation classes of the case (MAXVEG of CATHY.H must not be smaller)
    nveg = len(np.unique(np.loadtxt(os.path.join(proj, "input", "root_map"), skiprows=6)))
    simu.update_cathyH(MAXVEG=max(nveg, 2))

    simu.run_preprocessor(verbose=False)
    simu.run_processor(IPRT1=2, TMAX=t_end, TIMPRTi=t_prt, read_grid=False, verbose=False)


def compare_and_plot(a, src_case, proj, work, z, nrow_nodes, ncol_nodes):

    ref_sw = os.path.join(src_case, "output", "sw")
    new_sw = os.path.join(proj, "output", "sw")
    stats, t, R, N = compare_nodal(ref_sw, new_sw, out_CT.read_sw, t_max_days=a.days)
    stats["days_requested"] = a.days
    print(json.dumps(stats, indent=1))
    json.dump(stats, open(os.path.join(work, "venezia_levelA_summary.json"), "w"), indent=1)

    # probes: column of nodes in the middle of the column experiment, depths 10, 30, 50 cm
    nnod2d = nrow_nodes * ncol_nodes
    centre = (nrow_nodes // 2) * ncol_nodes + ncol_nodes // 2
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for depth, col in ((0.10, "#1f77b4"), (0.30, "#2ca02c"), (0.50, "#d62728")):
        k = int(np.argmin(abs(z - depth)))
        j = k * nnod2d + centre
        ax[0].plot(t, R[:, j], color=col, lw=2, label="reference %.0f cm" % (z[k] * 100))
        ax[0].plot(t, N[:, j], color=col, ls="--", label="wrapper %.0f cm" % (z[k] * 100))
    ax[0].set_xlabel("days"); ax[0].set_ylabel("S_w [-]"); ax[0].legend(fontsize=7, ncol=2)
    ax[0].set_title("Venezia, level A (own old src): reference vs wrapper run")
    ax[1].hist(np.abs(N - R).ravel(), bins=60, log=True)
    ax[1].set_xlabel("|S_w difference|"); ax[1].set_title("all nodes, all times (max %.2g)" % stats["max_abs_diff"])
    fig.tight_layout(); fig.savefig(os.path.join(work, "venezia_levelA_sw.png"), dpi=120)


if __name__ == "__main__":
    main()
