"""Run the REAL perturbate()/perturbate_parm() of the cloned pyCATHY on each scenario
(stub simulation, no CATHY executable needed)."""
import sys, tomllib, warnings
import numpy as np, pandas as pd
from pyCATHY.DA import perturbate
from pyCATHY.DA.scenario import validate_scenario

IDX = pd.MultiIndex.from_product([[1], [0]], names=["zone", "layer"])
soil = pd.DataFrame({"PERMX":1e-4,"PERMY":1e-4,"PERMZ":1e-4,"POROS":0.4}, index=IDX)
feddes = pd.DataFrame({"ZROOT":[0.4]}, index=pd.Index([1]))
class Sim:
    soil_dict = {k:[v] for k,v in dict(PCWLT=-150.0, OMGC=1.0, PZ=1.0).items()}
    soil = soil_dict
    dem_parameters = {"base": 15.0}
    def read_inputs(self, key, **kw):
        return {"soil":(soil,feddes),"zone":(np.ones((2,2),int),{}),"root_map":(np.ones((2,2),int),{})}[key]

scen = tomllib.load(open("scenario.toml","rb"))["scenario"]
NENS = 64
for name, sc in scen.items():
    np.random.seed(42)
    print(f"\n=== {name}")
    errs, warns = validate_scenario(sc, nlayers=1, zone_raster=np.ones((2,2),int),
                                    veg_raster=np.ones((2,2),int), n_feddes_rows=1,
                                    raise_on_error=False)
    print("validator errors:", errs or "none")
    try:
        entries = perturbate.perturbate(Sim(), sc, NENS, validate_scenario=False)
        out = {}
        for e in entries:
            out = perturbate.perturbate_parm(out, parm=e, type_parm=e["type_parm"], mean=e["mean"],
                    sd=e["sd"], sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"],
                    per_type=e["per_type"], savefig=False)
        for k, v in out.items():
            x = np.asarray(v["ini_perturbation"], float)
            print(f"  {k:15s} mean={x.mean():9.4f} sd={x.std(ddof=1):8.4f} min={x.min():8.4f} max={x.max():8.4f}")
        upd = ["St. var."] + [d for d in out for l in sc["listUpdateParm"] if l in d]
        print("  updated keys:", upd[1:], "| noise:", sc.get("pert_sigma_Archie", "(script default)"))
    except Exception as exc:
        print("  FAILED:", type(exc).__name__, exc)
