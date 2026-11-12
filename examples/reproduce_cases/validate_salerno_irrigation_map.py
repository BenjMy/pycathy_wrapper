#!/usr/bin/env python3
"""
Regenerate the irrigation_map of the Salerno case with the wrapper tools from the tables of the repository
(R/Table_NEUMANN.csv hourly irrigation, other-dem/veg-type_map land use, OPAC/map_kc_*.csv monthly Kc)
and compare it block by block with the shipped other-dem/irrigation_map.

usage: python validate_salerno_irrigation_map.py [--repo PATH] [--write OUT_FILE]
"""
import argparse
import json
import os

import numpy as np

from case_tools import case_dir, get_repo, salerno_irrigation_from_tables
from pyCATHY.irrigation import (irrigation_map_from_series, read_irrigation_map,
                                select_irrigation_records, write_irrigation_map)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--repo"); ap.add_argument("--write")
    a = ap.parse_args()
    case = case_dir(get_repo(a.repo), "Salerno")

    t_gen, r_gen = salerno_irrigation_from_tables(case)
    land = np.loadtxt(os.path.join(case, "other-dem", "veg-type_map"), skiprows=6)

    ref = read_irrigation_map(os.path.join(case, "other-dem", "irrigation_map"), land.shape[0], land.shape[1])
    out = {"blocks_generated": int(len(t_gen)), "blocks_reference": int(len(ref["times"]))}
    n = min(len(t_gen), len(ref["times"]))
    out["times_equal"] = bool(len(t_gen) == len(ref["times"]) and np.array_equal(t_gen, ref["times"]))
    if len(t_gen) == len(ref["times"]):
        d = np.abs(r_gen - ref["rates"])
        scale = np.abs(ref["rates"]).max()
        out["max_abs_diff_m_s"] = float(d.max())
        out["max_rel_diff_to_peak"] = float(d.max() / scale)
        out["blocks_with_diff_gt_1e-6_of_peak"] = int((d.max(axis=(1, 2)) > 1e-6 * scale).sum())
        out["irrigated_cells_equal"] = bool(np.array_equal(r_gen > 0, ref["rates"] > 0))
        out["total_depth_generated_mm_per_unit_area"] = float((r_gen.sum(axis=(1, 2)) * np.diff(np.append(t_gen, t_gen[-1]))).sum() * 1000)
        out["total_depth_reference_mm_per_unit_area"] = float((ref["rates"].sum(axis=(1, 2)) * np.diff(np.append(ref["times"], ref["times"][-1]))).sum() * 1000)
    else:
        i = np.nonzero(t_gen[:n] != ref["times"][:n])[0]
        out["first_time_mismatch_index"] = int(i[0]) if len(i) else n
    print(json.dumps(out, indent=1))
    if a.write:
        write_irrigation_map(a.write, t_gen, r_gen, irr_sw=ref["irr_sw"], irr_flag=ref["irr_flag"])
        print("written", a.write)


if __name__ == "__main__":
    main()
