#!/usr/bin/env python3
"""Compare the irrigation module with the classic approach (irrigation added to rain) over the scenarios.

usage: RES=scen_runs/results PLOTS=plots python plot_scenarios.py
"""
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RES = os.environ.get("RES", "scen_runs/results")
PLOTS = os.environ.get("PLOTS", "plots")
POROS = 0.55
COL = {"module": "#1f77b4", "classic": "#d62728"}


def load():
    out = {}
    for f in sorted(glob.glob(os.path.join(RES, "*.npz"))):
        scen, mode = re.match(r"(.*)_(module|classic)\.npz", os.path.basename(f)).groups()
        out.setdefault(scen, {})[mode] = np.load(f, allow_pickle=True)
    return out


def fields_of(d):
    return [str(n) for n in d["field_names"]]


def irrigated(names):
    return [n for n in names if "rainfed" not in n]


def metrics(d, other=None):
    """other = results of the other mode of the same scenario (used to estimate the irrigation that
    reached the soil in the classic mode: potential flux of classic minus module, same rain/ET)."""
    names = fields_of(d)
    irr = irrigated(names)
    m = {
        "rz_mean_irrigated": float(np.mean([d[f"rz_{n}"].mean() for n in irr])),
        "stress_days_irrigated": float(np.mean([d[f"stress_days_{n}"] for n in irr])),
        "deep_gain_mm_irrigated": float(np.mean([(d[f"deep_{n}"][-1] - d[f"deep_{n}"][0]) * 1000 for n in irr])),
        "rz_mean_rainfed": float(np.mean([d[f"rz_{n}"].mean() for n in names if "rainfed" in n])) if len(irr) < len(names) else np.nan,
        "irr_nominal_mm": float(d["irr_nominal_mm"]),
        # module: measured by the solver (hgirrigation: available - not used); classic: the irrigation
        # given in the atmbc (not measured: the solver has no separate record of it)
        "irr_applied_mm": float(d["irr_available_mm"] - d["irr_unused_mm"]) if "irr_available_mm" in d else float(d["irr_nominal_mm"]),
        "irr_unused_mm": float(d["irr_unused_mm"]) if "irr_unused_mm" in d else 0.0,
        "runoff_mm": float(d["cum_ovl_mm"]),
        "recharge_mm": float(d["cum_rec_mm"]),
        "net_flux_mm": float(d["cum_act_mm"]),
    }
    return m


def main():
    os.makedirs(PLOTS, exist_ok=True)
    runs = load()
    if not runs:
        raise SystemExit("no results in " + RES)
    table = {s: {m: metrics(d, v.get("module" if m == "classic" else "classic")) for m, d in v.items()} for s, v in runs.items()}
    json.dump(table, open(os.path.join(PLOTS, "scenario_metrics.json"), "w"), indent=1)

    # ---- figure 5: root-zone saturation module vs classic (every scenario that has both)
    both = [s for s, v in runs.items() if len(v) == 2]
    for s in both:
        d = runs[s]["module"]
        names = fields_of(d)
        fig, axs = plt.subplots(1, len(names), figsize=(4.6 * len(names), 3.8), sharey=True)
        axs = np.atleast_1d(axs)
        for a, n in zip(axs, names):
            for mode in ("module", "classic"):
                dd = runs[s][mode]
                a.plot(dd["t_days"] / 1.0, dd[f"rz_{n}"], color=COL[mode], ls="-" if mode == "module" else "--", label=mode)
            a.axhline(0.45, color="k", ls=":", lw=0.8)
            a.set_title(n, fontsize=10); a.set_xlabel("days from start of window")
        axs[0].set_ylabel("root-zone S_w [-]"); axs[0].legend()
        fig.suptitle(f"{s}: irrigation module (solid) vs classic, irrigation added to rain (dashed); dotted = stress threshold")
        fig.tight_layout(); fig.savefig(os.path.join(PLOTS, f"fig5_rootzone_{s}.png"), dpi=120); plt.close(fig)

    # ---- figure 6: summary over scenarios
    keys = [("rz_mean_irrigated", "mean root-zone S_w, irrigated fields [-]"),
            ("stress_days_irrigated", "stress days per irrigated field [d]"),
            ("deep_gain_mm_irrigated", "water gained below the roots [mm]"),
            ("irr_applied_mm", "irrigation reaching the soil [mm over district]\n(module: measured, hgirrigation; classic: as given in atmbc)"),
            ("runoff_mm", "cumulative overland flow (OVL) [mm]"),
            ("net_flux_mm", "cumulative net atmospheric exchange (ACT) [mm]")]
    scen = list(table)
    fig, axs = plt.subplots(2, 3, figsize=(17, 8))
    x = np.arange(len(scen))
    for a, (k, ttl) in zip(axs.ravel(), keys):
        for i, mode in enumerate(("module", "classic")):
            vals = [table[s].get(mode, {}).get(k, np.nan) for s in scen]
            a.bar(x + (i - 0.5) * 0.38, vals, 0.38, color=COL[mode], label=mode)
        a.set_xticks(x); a.set_xticklabels(scen, rotation=25, ha="right", fontsize=8); a.set_title(ttl, fontsize=10)
    axs[0, 0].legend()
    fig.suptitle("Scenarios: irrigation module vs classic (irrigation added to rain)")
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "fig6_scenarios_summary.png"), dpi=120); plt.close(fig)

    # ---- text table
    lines = ["| scenario | mode | rz Sw irr | stress d | deep gain mm | irrigation mm: module measured / classic as given (map nominal) | not used mm | OVL mm | ACT mm |", "|---|---|---|---|---|---|---|---|---|"]
    for s in scen:
        for mode in ("module", "classic"):
            if mode in table[s]:
                m = table[s][mode]
                lines.append(f"| {s} | {mode} | {m['rz_mean_irrigated']:.3f} | {m['stress_days_irrigated']:.1f} | {m['deep_gain_mm_irrigated']:.1f} | {m['irr_applied_mm']:.1f} ({m['irr_nominal_mm']:.1f}) | {m['irr_unused_mm']:.1f} | {m['runoff_mm']:.2f} | {m['net_flux_mm']:.1f} |")
    open(os.path.join(PLOTS, "scenario_table.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
