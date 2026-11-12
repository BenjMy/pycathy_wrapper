#!/usr/bin/env python3
"""
Venezia case, level C: the ORIGINAL inputs of the case (old source, old layout `input/` + `prepro/`)
translated to the most updated source (`code/`, version "withIrr", layout `other-dem/`, `prepro-dem/`,
file names `atmbc-flow`, `ic-flow`, ...) and run with the wrapper; the nodal saturation is compared with the
reference run of the repository (which was made with the OLD source).

Translation (wrapper only, nothing edited by hand):
  * DEM, layers, zones, vegetation map, hap.in, dem_parameters : rebuilt with the wrapper API
    (update_dem_parameters / update_prepo_inputs / update_zone / update_veg_map) -> files in the new format
  * parm : read with the wrapper (old layout, glued label...) and written in the new layout (no NCOUT)
  * soil, ic, atmbc, nansfdirbc, nansfneubc, sfbc, retctab, nudging, base_map, effraininp,
    livelli_iniz_s, posizione_serb : same format in both layouts, copied to the location given by the
    new project's cathy.fnames (simu._input_path)

usage: python reproduce_venezia_levelC.py [--repo PATH] [--work work_C] [--days 8]
"""
import argparse
import json
import os
import shutil
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from case_tools import case_dir, compare_nodal, get_repo, stage_case
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT

COPY_AS_IS = ["soil", "ic", "atmbc", "nansfdirbc", "nansfneubc", "sfbc", "retctab", "nudging",
              "base_map", "effraininp", "livelli_iniz_s", "posizione_serb"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo"); ap.add_argument("--work", default="work_venezia_C")
    ap.add_argument("--days", type=float, default=8.0)
    a = ap.parse_args()
    repo = get_repo(a.repo)
    src_case = case_dir(repo, "Venezia")
    work = os.path.abspath(a.work); os.makedirs(work, exist_ok=True)
    old = stage_case(src_case, os.path.join(work, "old"))          # inputs of the case (never its outputs)
    t_end = a.days * 86400.0

    # ---- what the old inputs say
    dem_par = in_CT.read_dem_parameters(os.path.join(old, "input", "dem_parameters"))
    zkey = next(k for k in dem_par if k.startswith("zratio"))
    dem = np.loadtxt(os.path.join(old, "prepro", "dem"), skiprows=6)
    veg = np.loadtxt(os.path.join(old, "input", "root_map"), skiprows=6)
    parm = in_CT.read_parm(os.path.join(old, "input", "parm"))
    t_prt = [t for t in np.atleast_1d(parm["(TIMPRT(I),I=1,NPRT)"]) if t <= t_end]
    dx, dy = float(dem_par["delta_x"]), float(dem_par["delta_y"])

    # ---- new project with the most updated source (fetched from GitHub by the wrapper)
    simu = cathy_tools.CATHY(dirName=os.path.join(work, "new"), prj_name="venezia", version="withIrr")
    simu.update_dem_parameters(delta_x=dx, delta_y=dy, nstr=int(dem_par["nstr"]), base=float(dem_par["base"]),
                               zratio=list(np.atleast_1d(dem_par[zkey])), nzone=int(dem_par["nzone"]),
                               n1=int(dem_par["n1"]), isp=int(dem_par["isp"]), ivert=int(dem_par["ivert"]))
    simu.update_prepo_inputs(dem, delta_x=dx, delta_y=dy)
    simu.update_zone(np.ones(dem.shape, dtype=int))
    simu.update_veg_map(veg)

    # ---- physical inputs: same format, new names/locations
    for name in COPY_AS_IS:
        src = os.path.join(old, "input", name)
        if os.path.isfile(src):
            shutil.copy(src, simu._input_path(name))

    # ---- parm: values of the case in the new layout (NCOUT does not exist anymore)
    keep = {k: v for k, v in parm.items()
            if k in simu.parm and k not in ("NPRT", "(TIMPRT(I),I=1,NPRT)", "TMAX")}
    nveg = int(len(np.unique(veg)))
    simu.update_cathyH(MAXVEG=max(nveg, 2))
    simu.run_preprocessor(verbose=False)
    simu.run_processor(IPRT1=2, TRAFLAG=0, TMAX=t_end, TIMPRTi=t_prt, verbose=False, **{
        k: v for k, v in keep.items() if k not in ("IPRT1", "TRAFLAG")})

    ref_sw = os.path.join(src_case, "output", "sw")
    new_sw = os.path.join(simu.workdir, simu.project_name, "output", "sw")
    stats, t, R, N = compare_nodal(ref_sw, new_sw, out_CT.read_sw, t_max_days=a.days)
    print(json.dumps(stats, indent=1))
    json.dump(stats, open(os.path.join(work, "venezia_levelC_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
