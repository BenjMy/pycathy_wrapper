"""
Scenario library for the synthetic irrigation district (flat test case), wrapper only.

Two ways to represent the same irrigation are compared:
  * "module"  : irrigation module of the most updated source (version="withIrr"): `irrigation_map`
                (IRR_FLAG=1, IRR_SW threshold), homogeneous `atmbc` = rain - potential ET
  * "classic" : no irrigation module: irrigation is added to the rain in a spatially variable
                `atmbc` (HSPATM=0, per surface node): net = rain + irrigation - potential ET
"""
import json
import os
import shutil
import threading
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from pyCATHY import cathy_tools
from pyCATHY.irrigation import (build_irrigation_map, cell_raster_to_node_rates,
                                surface_node_areas, write_irrigation_map)
from run_synthetic_district import (BASE, DAY, DX, KC_POT, NCOL, NROW, NSTR, ZRATIO,
                                    synthetic_forcing)

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2 renamed trapz
POROS = 0.55                       # default porosity of the wrapper soil (water content = POROS * Sw)
SW_STRESS = 0.45                   # root-zone saturation below which a crop is considered stressed
SEASON = dict(start_doy=105, end_doy=258)
TEMPLATES = {
    "drip":      dict(interval_days=1,  start_hour=6,  duration_h=2.0, rate_mm_h=1.5),   # 3 mm/d
    "sprinkler": dict(interval_days=3,  start_hour=18, duration_h=4.0, rate_mm_h=5.0),   # 20 mm / 3 d
    "furrow":    dict(interval_days=10, start_hour=8,  duration_h=6.0, rate_mm_h=8.0),   # 48 mm / 10 d
}


def _mask(rows, cols):
    m = np.zeros((NROW, NCOL), bool)
    m[rows[0]:rows[1], cols[0]:cols[1]] = True
    return m


def make_layout(name):
    """Fields of a district: {field: dict(mask, schedule)} (schedule None = rainfed)."""
    if name == "blocks4":      # 4 blocks of 6 x 6 cells
        return {"A drip": dict(mask=_mask((1, 7), (1, 7)), schedule="drip"),
                "B sprinkler": dict(mask=_mask((1, 7), (9, 15)), schedule="sprinkler"),
                "C furrow": dict(mask=_mask((9, 15), (1, 7)), schedule="furrow"),
                "D rainfed": dict(mask=_mask((9, 15), (9, 15)), schedule=None)}
    if name == "strips":       # 3 long strips of 5 columns
        return {"A drip": dict(mask=_mask((1, 15), (1, 6)), schedule="drip"),
                "B sprinkler": dict(mask=_mask((1, 15), (6, 11)), schedule="sprinkler"),
                "D rainfed": dict(mask=_mask((1, 15), (11, 16)), schedule=None)}
    if name == "checker":      # checkerboard of 3 x 3 cell plots, sprinkler / rainfed
        irr, rain = np.zeros((NROW, NCOL), bool), np.zeros((NROW, NCOL), bool)
        for i in range(4):
            for j in range(4):
                blk = _mask((2 + 3 * i, 5 + 3 * i), (2 + 3 * j, 5 + 3 * j))
                (irr if (i + j) % 2 == 0 else rain).__ior__(blk)
        return {"E irrigated plots": dict(mask=irr, schedule="sprinkler"),
                "F rainfed plots": dict(mask=rain, schedule=None)}
    raise ValueError(name)


def irrigation_fields(layout, timing=None):
    out = []
    for fld in layout.values():
        if fld["schedule"] is None:
            continue
        spec = dict(TEMPLATES[fld["schedule"]], **SEASON)
        if timing == "night":
            spec["start_hour"] = 22
        elif timing == "midday":
            spec["start_hour"] = 12
        out.append(dict(mask=fld["mask"], **spec))
    return out


def shift_irrigation(times, rates, t0, t1):
    """Window [t0, t1) of a piecewise constant irrigation map, time origin moved to t0."""
    i0 = int(np.searchsorted(times, t0, side="right") - 1)
    keep = [i0] + [i for i in range(i0 + 1, len(times)) if times[i] < t1]
    t = np.maximum(times[keep] - t0, 0.0)
    t[0] = 0.0
    return t, rates[keep]


def _guard_hgflag(path, stop):
    """hgflag logs every air-dry surface node at every step (GBs in dry summers): keep it sparse."""
    while not stop.is_set():
        try:
            os.truncate(path, 0)
        except OSError:
            pass
        stop.wait(10)


def run_case(name, mode, layout="blocks4", zroot=0.6, timing=None, irr_sw=0.8, d0=90, d1=285,
             out_root="scen_runs", seed=42, print_every=2, keep_project=False):
    """Prepare (wrapper only) and run one case; returns the path of the small result file."""
    assert mode in ("module", "classic")
    out_root = os.path.abspath(out_root)          # the wrapper changes the working directory
    res_dir = os.path.join(out_root, "results"); os.makedirs(res_dir, exist_ok=True)
    os.chdir(out_root)                            # never stay inside a project that is deleted below
    res_file = os.path.join(res_dir, f"{name}_{mode}.npz")
    if os.path.exists(res_file):
        return res_file                                       # resumable
    nd = d1 - d0
    t_end = nd * DAY
    rain_mm, et0_mm = synthetic_forcing(365, seed)
    rain_mm, et0_mm = rain_mm[d0:d1], et0_mm[d0:d1]
    fields = make_layout(layout)
    irr_t, irr_r = build_irrigation_map(irrigation_fields(fields, timing), (NROW, NCOL), 365 * DAY)
    irr_t, irr_r = shift_irrigation(irr_t, irr_r, d0 * DAY, d1 * DAY)

    prj_root = os.path.join(out_root, f"{name}_{mode}")
    shutil.rmtree(prj_root, ignore_errors=True)
    simu = cathy_tools.CATHY(dirName=prj_root, prj_name="district", version="withIrr")
    ii, jj = np.meshgrid(np.arange(NROW), np.arange(NCOL), indexing="ij")
    DEM = 100.0 - 1e-4 * DX * (ii + jj) / 2.0
    simu.update_dem_parameters(delta_x=DX, delta_y=DX, nstr=NSTR, base=BASE, zratio=ZRATIO)
    simu.update_prepo_inputs(DEM, delta_x=DX, delta_y=DX)
    simu.update_zone(np.ones((NROW, NCOL), dtype=int))

    # one crop class per field (+ one for the margins), root depth = zroot for the crops
    veg = np.full((NROW, NCOL), len(fields) + 1.0)
    for k, fld in enumerate(fields.values(), 1):
        veg[fld["mask"]] = k
    simu.update_veg_map(veg)
    fp = simu.set_SOIL_defaults(FP_map_default=True)
    fp["ZROOT"] = float(zroot)
    fp.loc[len(fields) + 1, "ZROOT"] = 0.3
    simu.update_soil(FP_map=fp)

    # forcing
    et_ms = et0_mm * KC_POT / 1000.0 / DAY
    rain_ms = rain_mm / 1000.0 / DAY
    times = np.arange(nd) * DAY
    if mode == "module":
        simu.update_atmbc(HSPATM=1, IETO=1, time=np.append(times, t_end),
                          VALUE=[np.append(rain_ms, rain_ms[-1]), np.append(et_ms, et_ms[-1])])
        write_irrigation_map(simu._input_path("irrigation_map"), irr_t, irr_r, irr_sw=irr_sw, irr_flag=1)
    else:
        node_irr = cell_raster_to_node_rates(irr_r, DX, DX)               # (n_maps, n_nodes)
        t_all = np.unique(np.concatenate([times, irr_t, [t_end]]))
        net = []
        for t in t_all:
            d = min(int(t // DAY), nd - 1)
            k = int(np.searchsorted(irr_t, t, side="right") - 1)
            net.append(rain_ms[d] - et_ms[d] + node_irr[k])
        simu.update_atmbc(HSPATM=0, IETO=1, time=t_all, netValue=net)
    simu.update_ic(INDP=0, IPOND=0, pressure_head_ini=-1.5)

    nnod2d = (NROW + 1) * (NCOL + 1)
    t_print = np.arange(print_every, nd + 1e-9, print_every) * DAY
    simu.update_parm(TIMPRTi=list(t_print), VTKF=0, TMAX=t_end, IPRT=4,
                     NUMVP=1, **{"(NODVP(I),I=1,NUMVP)": [nnod2d + 1]})
    simu.run_preprocessor(verbose=False)

    out_dir = os.path.join(simu.workdir, simu.project_name, "output")
    stop = threading.Event()
    th = threading.Thread(target=_guard_hgflag, args=(os.path.join(out_dir, "hgflag"), stop), daemon=True)
    th.start()
    try:
        simu.run_processor(IPRT1=2, TRAFLAG=0, DELTAT=60.0, DTMIN=1e-4, DTMAX=3600.0, TMAX=t_end, verbose=False)
    finally:
        stop.set(); th.join(timeout=15)

    metrics = postprocess(simu, out_dir, fields, zroot, irr_t, irr_r, nd)
    np.savez(res_file, **metrics)
    os.chdir(out_root)
    if not keep_project:
        shutil.rmtree(prj_root, ignore_errors=True)
    return res_file


def last_per_step(arr):
    """Rows of a solver time-series file with repeated step numbers (a rejected step is redone):
    keep the last line of every step."""
    last = {}
    for i, step in enumerate(arr[:, 0].astype(int)):
        last[step] = i
    return arr[sorted(last.values())]


def postprocess(simu, out_dir, fields, zroot, irr_t, irr_r, nd):
    """Small set of derived series/numbers from the solver outputs (what the comparison plots use)."""
    nnod2d = (NROW + 1) * (NCOL + 1)
    z = np.append(0.0, np.cumsum(ZRATIO) * BASE)
    sw = simu.read_outputs("sw")
    sw = sw[0] if isinstance(sw, tuple) else sw
    t_days = sw.index.values.astype(float) / DAY
    out = {"t_days": t_days, "z": z, "field_names": np.array(list(fields))}
    zz = np.linspace(0.0, zroot, 60)
    for name, fld in fields.items():
        cells = np.argwhere(fld["mask"])
        nodes = sorted({(r + dr) * (NCOL + 1) + (c + dc) for r, c in cells for dr in (0, 1) for dc in (0, 1)})
        prof = np.stack([sw.iloc[:, [k * nnod2d + n for n in nodes]].values.mean(1) for k in range(len(z))], axis=1)  # (nt, nlayer)
        rz = np.array([_trapz(np.interp(zz, z, p), zz) / zroot for p in prof])
        zd = np.linspace(zroot, z[-1], 80)
        deep = np.array([POROS * _trapz(np.interp(zd, z, p), zd) for p in prof])          # m of water below the roots
        out[f"rz_{name}"] = rz
        out[f"deep_{name}"] = deep
        out[f"profile_end_{name}"] = prof[-1]
        out[f"stress_days_{name}"] = float(np.mean(rz < SW_STRESS) * nd)
    area = NROW * NCOL * DX * DX
    h = last_per_step(np.loadtxt(os.path.join(out_dir, "hgatmsf"), comments="#", ndmin=2))
    dt = h[:, 1]
    for col, key in ((3, "pot"), (4, "act"), (5, "ovl"), (7, "seep"), (8, "rec")):
        out[f"cum_{key}_mm"] = float(np.sum(h[:, col] * dt) * 1000.0 / area)
    out["rec_vol_mm"] = float(h[-1, 9] * 1000.0 / area)
    nominal = sum(irr_r[i].sum() * DX * DX * ((irr_t[i + 1] if i + 1 < len(irr_t) else nd * DAY) - irr_t[i])
                  for i in range(len(irr_t)))
    out["irr_nominal_mm"] = float(nominal * 1000.0 / area)
    hg_file = os.path.join(out_dir, "hgirrigation")
    if os.path.exists(hg_file) and os.path.getsize(hg_file) > 0:
        g = last_per_step(np.loadtxt(hg_file, comments="#", ndmin=2))
        if g.size:
            out["irr_available_mm"] = float(np.sum(g[:, 4] * g[:, 1]) * 1000.0 / area)
            out["irr_unused_mm"] = float(np.sum(g[:, 6] * g[:, 1]) * 1000.0 / area)
    return out
