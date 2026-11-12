#!/usr/bin/env python3
"""
Synthetic irrigation district (flat test case) prepared and run ONLY through the pyCATHY wrapper,
with the most updated irrigation source (version="withIrr" = CATHY_src_WithIrrigation/code/).

 * 16 x 16 cells of 25 m (400 m x 400 m = 16 ha), quasi flat terrain, homogeneous soil, 12 layers
 * four fields of 6 x 6 cells: A drip, B sprinkler, C furrow, D rainfed (control)
 * typical Mediterranean forcing (daily rain and ET0, several years) -> atmbc
 * irrigation pattern -> other-dem/irrigation_map (IRR_FLAG=1, IRR_SW, IRR_HS=0, TIMEIN + rasters)

usage:  YEARS=3 OUT=synthetic_district python run_synthetic_district.py
"""
import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from pyCATHY import cathy_tools
from pyCATHY.irrigation import build_irrigation_map, write_irrigation_map

YEARS = float(os.environ.get("YEARS", "3"))
DAYS = int(round(YEARS * 365))
OUT = os.environ.get("OUT", "synthetic_district")
NROW = NCOL = 16
DX = DY = 25.0
NSTR, BASE = 12, 6.0
ZRATIO = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.18, 0.24]  # sum = 1
IPRT = int(os.environ.get("IPRT", "4"))            # solver print option (4 also writes velocity fields: heavy)
PRINT_EVERY_DAYS = float(os.environ.get("PRINT_EVERY_DAYS", "1"))
IRR_SW = 0.80          # irrigation is applied at a surface node only if its saturation is below
KC_POT = 0.7           # potential evapotranspiration = KC_POT * ET0 (homogeneous atmbc)
DAY = 86400.0


# --------------------------------------------------------------------------- forcing
def synthetic_forcing(ndays, seed=42):
    """Typical Mediterranean daily forcing: wet winters, dry hot summers (~700 mm/y, ET0 ~1100 mm/y)."""
    rng = np.random.default_rng(seed)
    doy = np.arange(ndays) % 365
    wet_season = 0.5 * (1 + np.cos(2 * np.pi * (doy - 350) / 365))     # 1 mid-December, 0 mid-June
    sunny = 0.5 * (1 + np.cos(2 * np.pi * (doy - 196) / 365))         # 1 mid-July, 0 mid-January
    wet_day = rng.random(ndays) < (0.10 + 0.28 * wet_season)
    rain_mm = np.where(wet_day, rng.exponential(4.0 + 8.0 * wet_season), 0.0)
    et0_mm = (0.8 + 5.4 * sunny ** 1.2) * np.clip(1 + 0.12 * rng.standard_normal(ndays), 0.5, 1.6)
    et0_mm = np.where(rain_mm > 0, 0.5 * et0_mm, et0_mm)               # less demand on rainy days
    return rain_mm, et0_mm


# --------------------------------------------------------------------------- fields
def make_fields():
    def mask(r, c):
        m = np.zeros((NROW, NCOL), bool)
        m[r[0]:r[1], c[0]:c[1]] = True
        return m
    season = dict(start_doy=105, end_doy=258)        # 15 April - 15 September
    return {
        "A drip":      dict(mask=mask((1, 7), (1, 7)),   interval_days=1,  start_hour=6,  duration_h=2.0, rate_mm_h=1.5, **season),  # 3 mm/d
        "B sprinkler": dict(mask=mask((1, 7), (9, 15)),  interval_days=3,  start_hour=18, duration_h=4.0, rate_mm_h=5.0, **season),  # 20 mm / 3 d
        "C furrow":    dict(mask=mask((9, 15), (1, 7)),  interval_days=10, start_hour=8,  duration_h=6.0, rate_mm_h=8.0, **season),  # 48 mm / 10 d
        "D rainfed":   dict(mask=mask((9, 15), (9, 15)), interval_days=1,  start_hour=0,  duration_h=0.0, rate_mm_h=0.0, start_doy=1, end_doy=0),
    }


def main():
    rain_mm, et0_mm = synthetic_forcing(DAYS)
    fields = make_fields()
    t_end = DAYS * DAY
    irr_times, irr_rates = build_irrigation_map(list(fields.values()), (NROW, NCOL), t_end)

    simu = cathy_tools.CATHY(dirName=OUT, prj_name="district", version="withIrr")

    # quasi flat DEM: 0.01 % slope towards the south-east (a strictly flat DEM has no outlet)
    ii, jj = np.meshgrid(np.arange(NROW), np.arange(NCOL), indexing="ij")
    DEM = 100.0 - 1e-4 * DX * (ii + jj) / 2.0
    # layers first: update_prepo_inputs sizes CATHY.H (MAXSTR) from the number of layers
    simu.update_dem_parameters(delta_x=DX, delta_y=DY, nstr=NSTR, base=BASE, zratio=ZRATIO)
    simu.update_prepo_inputs(DEM, delta_x=DX, delta_y=DY)
    simu.update_zone(np.ones((NROW, NCOL), dtype=int))
    simu.update_veg_map(np.ones((NROW, NCOL)))
    simu.update_soil()

    # atmospheric forcing: net flux (rain - potential evapotranspiration), m/s, daily steps
    times = np.append(np.arange(DAYS) * DAY, t_end)
    rain_ms = np.append(rain_mm, rain_mm[-1]) / 1000.0 / DAY
    et_ms = np.append(et0_mm, et0_mm[-1]) * KC_POT / 1000.0 / DAY
    simu.update_atmbc(HSPATM=1, IETO=1, time=times, VALUE=[rain_ms, et_ms])
    simu.update_ic(INDP=0, IPOND=0, pressure_head_ini=-2.0)

    write_irrigation_map(simu._input_path("irrigation_map"), irr_times, irr_rates, irr_sw=IRR_SW, irr_flag=1)

    # vertical profiles (soil moisture time series) at the centre of each field, 3 depths
    nnod2d = (NROW + 1) * (NCOL + 1)
    centres = {"A": (4, 4), "B": (4, 12), "C": (12, 4), "D": (12, 12)}      # node (row, col)
    depth_layers = [1, 3, 5]                                               # 0.06, 0.36, 0.90 m
    nodvp = [k * nnod2d + r * (NCOL + 1) + c + 1 for (r, c) in centres.values() for k in depth_layers]
    t_print = np.arange(PRINT_EVERY_DAYS, DAYS + 1e-9, PRINT_EVERY_DAYS) * DAY
    simu.update_parm(TIMPRTi=list(t_print), VTKF=0, TMAX=t_end, IPRT=IPRT,
                     NUMVP=len(nodvp), **{"(NODVP(I),I=1,NUMVP)": nodvp})
    simu.run_preprocessor(verbose=False)
    simu.run_processor(IPRT1=2, TRAFLAG=0, DELTAT=60.0, DTMIN=1e-4, DTMAX=21600.0, TMAX=t_end, verbose=False)

    os.makedirs(os.path.join(OUT, "district", "forcing"), exist_ok=True)
    np.savez(os.path.join(OUT, "district", "forcing", "forcing.npz"), rain_mm=rain_mm, et0_mm=et0_mm,
             irr_times=irr_times, irr_rates=irr_rates, **{f"mask_{k.split()[0]}": v["mask"] for k, v in fields.items()})
    json.dump(dict(days=DAYS, nrow=NROW, ncol=NCOL, dx=DX, nstr=NSTR, base=BASE, zratio=ZRATIO,
                   irr_sw=IRR_SW, kc_pot=KC_POT, fields=list(fields), nodvp=nodvp,
                   centres={k: list(v) for k, v in centres.items()}, depth_layers=depth_layers),
              open(os.path.join(OUT, "district", "forcing", "setup.json"), "w"), indent=1)
    print("RUN DONE")


if __name__ == "__main__":
    main()
