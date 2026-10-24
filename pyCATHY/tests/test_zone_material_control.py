"""'zone_material' control (zone x group of layers) and Ks on the shared perturbation builder.

Why: 'layers' / 'material' write a value into EVERY zone of the layers they control, so the
compacted zone of a two-zone model (HC) is overwritten by the prior of zone 1. 'zone_material'
addresses [zone, [layers]] pairs only.
"""
import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.param_tables import _assign_soil_value
from pyCATHY.DA.scenario import validate_scenario

NL = 4
INDEX = pd.MultiIndex.from_product([[1, 2], range(NL)], names=["zone", "layer"])
ZGROUPS = [[1, [0, 1]], [1, [2, 3]], [2, [0, 1]]]      # zone 2 below layer 1 is not controlled
NEN = 8


def make_sim():
    soil = pd.DataFrame({"PERMX": 1e-6, "PERMY": 1e-6, "PERMZ": 1e-6, "POROS": 0.4}, index=INDEX)
    feddes = pd.DataFrame({"ZROOT": [0.4]}, index=pd.Index([1]))
    zone = np.ones((4, 4)); zone[:, 1:3] = 2

    class Sim:
        dem_parameters = {"base": 15.0}
        soil_dict = dict(PCWLT=[-150.0], OMGC=[1.0], PZ=[1.0])

        def read_inputs(self, key, **kw):
            return {"soil": (soil, feddes), "zone": (zone, {}), "root_map": (np.ones((4, 4)), {})}[key]

    Sim.soil = Sim.soil_dict
    return Sim()


def sc(name="Ks", control="zone_material", nom=(1e-6, 2e-6, 3e-7), sig=(1e-7, 2e-7, 3e-8), transf="log"):
    n = len(nom)
    return {
        "per_name": [name], "per_nom": [list(nom)], "per_mean": [list(nom)], "per_sigma": [list(sig)],
        "per_type": [["None"] * n], "sampling_type": ["normal"], "transf_type": [transf],
        "pert_control_name": [control], "zone_material_layers": ZGROUPS, "material_layers": [[0, 1], [2, 3]],
        "listUpdateParm": [f"{name}{k}" for k in range(n)],
    }


def test_assign_touches_only_its_zone_and_layers():
    df = pd.DataFrame(np.nan, index=INDEX, columns=["PERMX", "PERMY", "PERMZ"])
    _assign_soil_value(df, "zone_material", 2, ["PERMX", "PERMY", "PERMZ"], 5.0, layer_groups=ZGROUPS)
    hit = df.dropna(how="all").index.tolist()
    assert hit == [(2, 0), (2, 1)]                      # zone 1 untouched
    _assign_soil_value(df, "zone_material", 0, ["PERMX"], 7.0, layer_groups=ZGROUPS)
    assert df.loc[(1, 0), "PERMX"] == 7.0 and df.loc[(1, 1), "PERMX"] == 7.0 and np.isnan(df.loc[(1, 2), "PERMX"])
    assert df.loc[(2, 0), "PERMX"] == 5.0               # zone 2 keeps its own value
    # contrast: 'material' overwrites the same layers in BOTH zones
    df2 = pd.DataFrame(np.nan, index=INDEX, columns=["PERMX"])
    _assign_soil_value(df2, "material", 0, ["PERMX"], 7.0, layer_groups=[[0, 1], [2, 3]])
    assert df2.dropna().index.tolist() == [(1, 0), (1, 1), (2, 0), (2, 1)]
    with pytest.raises(ValueError):
        _assign_soil_value(df, "zone_material", 0, ["PERMX"], 1.0)


@pytest.mark.parametrize("name", ["Ks", "SATCELL_VG", "n_VG", "thetar_VG"])
def test_entries_one_per_zone_material(name):
    nom = (0.9, 0.8, 0.7) if name != "Ks" else (1e-6, 2e-6, 3e-7)
    entries = [e for e in perturbate.perturbate(make_sim(), sc(name, nom=nom, sig=tuple(0.1 * v for v in nom)),
                                                NEN, validate_scenario=False) if e["type_parm"].startswith(name)]
    assert [e["type_parm"] for e in entries] == [f"{name}0", f"{name}1", f"{name}2"]
    for e in entries:
        assert e["pert_control_name"] == "zone_material" and e["layer_groups"] == ZGROUPS
        assert e["transf_type"] == "log"


def test_ks_samples_and_builds_ensemble():
    entries = perturbate.perturbate(make_sim(), sc(), NEN, validate_scenario=False)
    out = {}
    np.random.seed(1)
    for e in entries:
        out = perturbate.perturbate_parm(out, parm=e, type_parm=e["type_parm"], mean=e["mean"], sd=e["sd"],
                                         sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"],
                                         per_type=e["per_type"], savefig=False)
    for k, nom in zip(("Ks0", "Ks1", "Ks2"), (1e-6, 2e-6, 3e-7)):
        v = np.ravel(out[k]["ini_perturbation"])
        assert v.size == NEN and abs(v.mean() / nom - 1) < 0.3 and v.std() > 0
        assert out[k]["pert_control_name"] == "zone_material" and out[k]["transf_type"] == "log"


def test_ks_legacy_scenarios_still_build():
    # scalar Ks without transf_type / pert_control_name (the former dedicated branch)
    s = {"per_name": ["Ks"], "per_nom": [1e-6], "per_mean": [1e-6], "per_sigma": [1e-7], "per_type": ["None"],
         "listUpdateParm": ["Ks"]}
    e = perturbate.perturbate(make_sim(), s, NEN, validate_scenario=False)
    assert [x["type_parm"] for x in e] == ["Ks"]
    # per layer, per zone, per material (new)
    s = sc(control="layers", nom=(1e-6,) * NL, sig=(1e-7,) * NL)
    s["listUpdateParm"] = [f"Ks{k}" for k in range(NL)]
    assert len(perturbate.perturbate(make_sim(), s, NEN, validate_scenario=False)) == NL
    s = sc(control="zone", nom=(1e-6, 3e-7), sig=(1e-7, 3e-8))
    s["listUpdateParm"] = ["Ks0", "Ks1"]
    assert len(perturbate.perturbate(make_sim(), s, NEN, validate_scenario=False)) == 2
    s = sc(control="material", nom=(1e-6, 2e-6), sig=(1e-7, 2e-7))
    s["listUpdateParm"] = ["Ks0", "Ks1"]
    e = perturbate.perturbate(make_sim(), s, NEN, validate_scenario=False)
    assert all(x["layer_groups"] == [[0, 1], [2, 3]] for x in e) and len(e) == 2   # groups shared, index selects


def test_validator():
    kw = dict(nlayers=NL, zone_raster=make_sim().read_inputs("zone")[0])
    errors, warns = validate_scenario(sc(), raise_on_error=False, **kw)
    assert not errors
    assert any("keep the base value" in w for w in warns)          # zone 2, layers 2-3 not controlled
    bad = sc(); bad["zone_material_layers"] = [[1, [0, 1]], [1, [1, 2]], [2, [0]]]
    assert any("several groups" in e for e in validate_scenario(bad, raise_on_error=False, **kw)[0])
    bad = sc(); bad["zone_material_layers"] = [[1, [0, 1]], [1, [2, 3]], [3, [0]]]
    assert any("zones must be" in e for e in validate_scenario(bad, raise_on_error=False, **kw)[0])
    bad = sc(); bad["zone_material_layers"] = [[1, [0, 1]], [1, [2, 9]], [2, [0]]]
    assert any("layers must be" in e for e in validate_scenario(bad, raise_on_error=False, **kw)[0])
    bad = sc(nom=(1e-6, 2e-6), sig=(1e-7, 2e-7))
    assert validate_scenario(bad, raise_on_error=False, **kw)[0]    # 2 entries, 3 groups
    bad = sc(); del bad["zone_material_layers"]
    assert any("zone_material_layers" in e for e in validate_scenario(bad, raise_on_error=False, **kw)[0])
    bad = sc("ic"); bad["per_nom"] = [[-1.0, -1.0, -1.0]]; bad["per_mean"] = bad["per_nom"]
    assert any("only supported" in e for e in validate_scenario(bad, raise_on_error=False, **kw)[0])
