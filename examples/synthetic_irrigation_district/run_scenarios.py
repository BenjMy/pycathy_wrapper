#!/usr/bin/env python3
"""
Run the scenario matrix (module vs classic). Resumable: finished cases are skipped.

usage: python run_scenarios.py [scenario ...]            (default: all, in priority order)
env  : OUT=scen_runs  D0=90 D1=285 TEST=1 (6-day smoke test)  PRINT_EVERY=2
"""
import os
import sys
import time

from district_lib import run_case

OUT = os.path.abspath(os.environ.get("OUT", "scen_runs"))  # absolute: the wrapper changes the cwd
D0, D1 = int(os.environ.get("D0", 170)), int(os.environ.get("D1", 230))
if os.environ.get("TEST"):
    D0, D1 = 104, 110

SCENARIOS = {                      # priority order (runs are resumable)
    "S0_base":      dict(layout="blocks4", zroot=0.6),
    "S1_root0.3":   dict(layout="blocks4", zroot=0.3),
    "S3_night":     dict(layout="blocks4", zroot=0.6, timing="night"),
    "S4_strips":    dict(layout="strips", zroot=0.6),
    "S2_root1.0":   dict(layout="blocks4", zroot=1.0),
    "S5_checker":   dict(layout="checker", zroot=0.6),
    "S6_thresh0.5": dict(layout="blocks4", zroot=0.6, irr_sw=0.5, modes=("module",)),
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(SCENARIOS)
    max_cases = int(os.environ.get("MAX_CASES", 10**6))     # e.g. 1: one new case per call (sessions are short)
    n_new = 0
    for name in names:
        spec = dict(SCENARIOS[name])
        for mode in spec.pop("modes", ("module", "classic")):
            if os.path.exists(os.path.join(OUT, "results", f"{name}_{mode}.npz")):
                continue
            if n_new >= max_cases:
                raise SystemExit(0)
            n_new += 1
            t0 = time.time()
            f = run_case(name, mode, d0=D0, d1=D1, out_root=OUT, print_every=float(os.environ.get("PRINT_EVERY", 2)), **spec)
            print(f"DONE {name} {mode} -> {f}  [{time.time() - t0:.0f} s]", flush=True)
