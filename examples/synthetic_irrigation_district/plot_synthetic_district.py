#!/usr/bin/env python3
"""Plots of the synthetic irrigation district (inputs, soil moisture, irrigation diagnostics, maps).

usage: OUT=synthetic_district PLOTS=plots python plot_synthetic_district.py
"""
import json
import os
import warnings
from datetime import datetime, timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

warnings.filterwarnings("ignore")
from pyCATHY import cathy_tools

OUT = os.environ.get("OUT", "synthetic_district")
PLOTS = os.environ.get("PLOTS", "plots")
DAY = 86400.0
T0 = datetime(2021, 1, 1)                     # synthetic calendar, 365-day years
COL = {"A": "#1f77b4", "B": "#2ca02c", "C": "#d62728", "D": "#7f7f7f"}
NAME = {"A": "A drip", "B": "B sprinkler", "C": "C furrow", "D": "D rainfed"}


def dates(t_s):
    return [T0 + timedelta(seconds=float(t)) for t in np.atleast_1d(t_s)]


def main():
    os.makedirs(PLOTS, exist_ok=True)
    prj = os.path.join(OUT, "district")
    F = np.load(os.path.join(prj, "forcing", "forcing.npz"))
    S = json.load(open(os.path.join(prj, "forcing", "setup.json")))
    nrow, ncol, dx, days = S["nrow"], S["ncol"], S["dx"], S["days"]
    area = nrow * ncol * dx * dx
    rain, et0 = F["rain_mm"], F["et0_mm"]
    et_pot = S["kc_pot"] * et0
    irr_t, irr_r = F["irr_times"], F["irr_rates"]
    masks = {k: F["mask_" + k] for k in "ABCD"}
    t_day = np.arange(days) * DAY
    years = int(round(days / 365))

    # irrigation depth per field per day [mm/d] from the piecewise constant map
    edges = np.append(irr_t, days * DAY)
    irr_day = {k: np.zeros(days) for k in "ABCD"}
    for i in range(len(irr_t)):
        for k in "ABC":
            rate = irr_r[i][masks[k]].mean() if masks[k].any() else 0.0
            if rate == 0:
                continue
            a, b = edges[i], edges[i + 1]
            for d in range(int(a // DAY), min(int(np.ceil(b / DAY)), days)):
                overlap = max(0.0, min(b, (d + 1) * DAY) - max(a, d * DAY))
                irr_day[k][d] += rate * overlap * 1000.0

    # ------------------------------------------------------------------ figure 1: inputs
    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    a = ax[0, 0]
    month = np.arange(days) // 30.4167
    nm = int(month.max()) + 1
    rain_m = np.array([rain[month == m].sum() for m in range(nm)])
    et_m = np.array([et_pot[month == m].sum() for m in range(nm)])
    x = np.arange(nm)
    a.bar(x - 0.2, rain_m, 0.4, color="#4c78a8", label="rain")
    a.bar(x + 0.2, et_m, 0.4, color="#f58518", label=f"potential ET (= {S['kc_pot']} ET0)")
    a.set_xticks(np.arange(0, nm, 6)); a.set_xticklabels([f"{T0.year + m // 12}-{m % 12 + 1:02d}" for m in range(0, nm, 6)], rotation=30)
    a.set_ylabel("mm / month"); a.set_title("Synthetic Mediterranean forcing (monthly)"); a.legend()
    a.text(0.01, 0.95, "annual mean: rain %.0f mm, ETpot %.0f mm" % (rain.sum() / years, et_pot.sum() / years), transform=a.transAxes, va="top", fontsize=9)

    a = ax[0, 1]
    layout = np.zeros((nrow, ncol))
    for i, k in enumerate("ABCD", 1):
        layout[masks[k]] = i
    cmap = ListedColormap(["#f0f0f0", COL["A"], COL["B"], COL["C"], COL["D"]])
    a.imshow(layout, cmap=cmap, vmin=0, vmax=4, extent=[0, ncol * dx, 0, nrow * dx], origin="upper")
    for k in "ABCD":
        r, c = np.argwhere(masks[k]).mean(0)
        a.text((c + 0.5) * dx, (nrow - r - 0.5) * dx, NAME[k], ha="center", va="center", color="white", fontsize=11, weight="bold")
    a.set_xlabel("x [m]"); a.set_ylabel("y [m]")
    a.set_title("District layout (16 x 16 cells of 25 m, quasi flat): fields and margins (grey)")

    a = ax[1, 0]
    for k in "ABC":
        a.fill_between(dates(t_day), irr_day[k], step="post", color=COL[k], alpha=0.8, label=NAME[k])
    a2 = a.twinx(); a2.bar(dates(t_day), rain, width=1.0, color="#4c78a8", alpha=0.35)
    a2.set_ylabel("rain [mm/d]", color="#4c78a8"); a2.invert_yaxis(); a2.set_ylim(60, 0)
    a.set_ylabel("nominal irrigation [mm/d]"); a.legend(loc="upper left", ncol=3)
    a.set_title("Irrigation pattern (irrigation_map) and rainfall (inverted axis)")
    a.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))

    a = ax[1, 1]
    labels = ["rain", "ET pot"] + [NAME[k] + "\nirrigation" for k in "ABC"]
    vals = [rain.sum() / years, et_pot.sum() / years] + [irr_day[k].sum() / years for k in "ABC"]
    a.bar(range(5), vals, color=["#4c78a8", "#f58518", COL["A"], COL["B"], COL["C"]])
    for i, v in enumerate(vals):
        a.text(i, v, "%.0f" % v, ha="center", va="bottom")
    a.set_xticks(range(5)); a.set_xticklabels(labels, fontsize=9); a.set_ylabel("mm / year (mean)")
    a.set_title("Annual totals")
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "fig1_forcing_and_layout.png"), dpi=130); plt.close(fig)

    # ------------------------------------------------------------------ outputs
    simu = cathy_tools.CATHY(dirName=OUT, prj_name="district", version="withIrr")
    sw = simu.read_outputs("sw")
    sw = sw[0] if isinstance(sw, tuple) else sw
    t_sw = sw.index.values.astype(float)
    nnod2d = (nrow + 1) * (ncol + 1)
    nodvp = S["nodvp"]; depth_layers = S["depth_layers"]
    cum_depth = np.cumsum(S["zratio"]) * S["base"]
    depth_of = {k: cum_depth[k - 1] for k in depth_layers}
    series = {}
    for i, name in enumerate("ABCD"):
        for j, k in enumerate(depth_layers):
            series[(name, k)] = sw.iloc[:, nodvp[i * len(depth_layers) + j] - 1].values

    # ------------------------------------------------------------------ figure 2: soil moisture
    fig, axs = plt.subplots(4, 1, figsize=(15, 12), sharex=True)
    shades = ["#9ecae1", "#3182bd", "#08306b"]
    for a, name in zip(axs, "ABCD"):
        for j, k in enumerate(depth_layers):
            a.plot(dates(t_sw), series[(name, k)], color=shades[j], lw=1.2, label="%.0f cm" % (depth_of[k] * 100))
        a.axhline(S["irr_sw"], color="k", ls=":", lw=0.8)
        for y in range(years):
            a.axvspan(T0 + timedelta(days=y * 365 + 104), T0 + timedelta(days=y * 365 + 258), color=COL[name] if name != "D" else "#bbbbbb", alpha=0.10)
        a.set_ylabel("S_w [-]"); a.set_title(NAME[name] + " (centre node)", loc="left", fontsize=10); a.set_ylim(0.1, 1.0)
        if name == "A":
            a.legend(ncol=3, loc="upper right", title="depth")
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    fig.suptitle("Simulated soil saturation; shaded = irrigation season; dotted = IRR_SW threshold %.2f" % S["irr_sw"])
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "fig2_soil_moisture.png"), dpi=130); plt.close(fig)

    # ------------------------------------------------------------------ figure 3: irrigation diagnostics (hgirrigation)
    h = np.loadtxt(os.path.join(prj, "output", "hgirrigation"), comments="#")
    th, atm, avail, unused = h[:, 2], h[:, 3], h[:, 4], h[:, 6]
    cum_atm, cum_avail, cum_unused = h[:, 7], h[:, 8], h[:, 10]
    to_mm = 1000.0 / area
    nominal = sum((irr_r[i].sum() * dx * dx) * (edges[i + 1] - edges[i]) for i in range(len(irr_t)))
    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    ax[0].plot(dates(th), cum_avail * to_mm, label="available irrigation (given in irrigation_map)")
    ax[0].plot(dates(th), (cum_avail - cum_unused) * to_mm, label="irrigation applied")
    ax[0].plot(dates(th), cum_unused * to_mm, label="irrigation not used (S_w >= IRR_SW)")
    ax[0].plot(dates(th), -cum_atm * to_mm, color="#f58518", label="cumulative potential ET - rain (ATMPOT)")
    ax[0].set_ylabel("mm over the whole district"); ax[0].legend(); ax[0].set_title("hgirrigation: cumulative volumes")
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    rel = (cum_avail[-1] - nominal) / nominal if nominal else np.nan
    ax[0].text(0.02, 0.60, "nominal (map) %.0f mm | solver available %.0f mm | diff %.2g %%" % (nominal * to_mm, cum_avail[-1] * to_mm, 100 * rel),
               transform=ax[0].transAxes, fontsize=9)
    frac = np.where(avail > 0, unused / np.where(avail > 0, avail, 1), np.nan)
    ax[1].plot(dates(th), frac, ".", ms=2, alpha=0.5)
    ax[1].set_ylabel("share of available irrigation not used"); ax[1].set_ylim(-0.02, 1.02)
    ax[1].set_title("Threshold effect (steps with irrigation only)"); ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "fig3_irrigation_diagnostics.png"), dpi=130); plt.close(fig)

    # ------------------------------------------------------------------ figure 4: maps and profiles
    pick = [np.argmin(abs(t_sw - d * DAY)) for d in (60, 200, 258, 330)]
    titles = ["day 60 (wet winter)", "day 200 (mid summer)", "day 258 (end of irrigation season)", "day 330 (autumn)"]
    fig, axs = plt.subplots(2, 4, figsize=(17, 8))
    for a, i, ttl in zip(axs[0], pick, titles):
        surf = sw.iloc[i, nnod2d * 1: nnod2d * 2].values.reshape(nrow + 1, ncol + 1)    # node layer 1 (6 cm)
        im = a.imshow(surf, vmin=0.2, vmax=1.0, cmap="YlGnBu", extent=[0, ncol * dx, 0, nrow * dx])
        a.set_title("S_w at 6 cm, " + ttl + "\n(t = %.0f d)" % (t_sw[i] / DAY), fontsize=9); a.set_xlabel("x [m]")
    fig.colorbar(im, ax=axs[0].tolist(), shrink=0.8, label="S_w [-]")
    centres = S["centres"]
    i_end = pick[2]
    for a, (name, (r, c)) in zip(axs[1], centres.items()):
        prof = [sw.iloc[i_end, k * nnod2d + r * (ncol + 1) + c] for k in range(len(S["zratio"]) + 1)]
        z = np.append(0, cum_depth)
        a.plot(prof, -z, "o-", color=COL[name]); a.set_ylim(-3.5, 0); a.set_xlim(0.1, 1.0)
        a.set_title(NAME[name] + ": profile at t = %.0f d" % (t_sw[i_end] / DAY), fontsize=9); a.set_xlabel("S_w [-]"); a.set_ylabel("depth [m]")
    fig.savefig(os.path.join(PLOTS, "fig4_maps_and_profiles.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # ------------------------------------------------------------------ numbers
    summary = {
        "years": years, "rain_mm_per_year": float(rain.sum() / years), "et_pot_mm_per_year": float(et_pot.sum() / years),
        "irrigation_nominal_mm_per_year_field": {NAME[k]: float(irr_day[k].sum() / years) for k in "ABC"},
        "irrigation_total_nominal_m3": float(nominal), "irrigation_available_solver_m3": float(cum_avail[-1]),
        "irrigation_not_used_m3": float(cum_unused[-1]),
        "irrigation_not_used_share": float(cum_unused[-1] / cum_avail[-1]) if cum_avail[-1] else None,
        "n_irrigation_map_times": int(len(irr_t)), "n_prints": int(len(t_sw)),
    }
    json.dump(summary, open(os.path.join(PLOTS, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
