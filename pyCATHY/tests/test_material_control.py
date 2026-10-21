"""'material' control: one value per group of layers (45 layer parameters -> a few materials)."""
import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.param_tables import _assign_soil_value
from pyCATHY.DA.scenario import validate_scenario

GROUPS = [[0, 1], [2], [3]]
INDEX = pd.MultiIndex.from_product([[1], [0, 1, 2, 3]], names=["zone", "layer"])
NEN = 8


def make_sim():
    soil = pd.DataFrame({"PERMX": 1e-4, "PERMY": 1e-4, "PERMZ": 1e-4, "POROS": 0.4}, index=INDEX)
    feddes = pd.DataFrame({"ZROOT": [0.4]}, index=pd.Index([1]))

    class Sim:
        dem_parameters = {"base": 15.0}
        soil_dict = dict(PCWLT=[-150.0], OMGC=[1.0], PZ=[1.0])

        def read_inputs(self, key, **kw):
            return {"soil": (soil, feddes), "zone": (np.ones((2, 2)), {}), "root_map": (np.ones((2, 2)), {})}[key]

    Sim.soil = Sim.soil_dict
    return Sim()


def sc(groups=GROUPS):
    return {
        "per_name": ["SATCELL_VG"], "per_nom": [[0.9, 0.8, 0.7]], "per_mean": [[0.9, 0.8, 0.7]],
        "per_sigma": [[0.09, 0.08, 0.07]], "per_type": [["None"] * 3],
        "sampling_type": ["normal"], "transf_type": ["None"],
        "pert_control_name": ["material"], "material_layers": groups,
        "listUpdateParm": ["SATCELL_VG0", "SATCELL_VG1", "SATCELL_VG2"],
    }


def test_entries_one_per_material():
    entries = perturbate.perturbate(make_sim(), sc(), NEN, validate_scenario=False)
    keys = [e["type_parm"] for e in entries if e["type_parm"].startswith("SATCELL_VG")]
    assert keys == ["SATCELL_VG0", "SATCELL_VG1", "SATCELL_VG2"]
    for e in entries:
        assert e["pert_control_name"] == "material"
        assert e["layer_groups"] == GROUPS


def test_assign_writes_every_layer_of_the_material():
    df = pd.DataFrame(np.nan, index=INDEX, columns=["VGPSATCELL"])
    _assign_soil_value(df, "material", 0, ["VGPSATCELL"], 5.0, layer_groups=GROUPS)
    _assign_soil_value(df, "material", 2, ["VGPSATCELL"], 7.0, layer_groups=GROUPS)
    assert df["VGPSATCELL"].tolist()[:2] == [5.0, 5.0]
    assert np.isnan(df["VGPSATCELL"].iloc[2]) and df["VGPSATCELL"].iloc[3] == 7.0


def test_assign_material_without_groups_fails():
    df = pd.DataFrame(np.nan, index=INDEX, columns=["VGPSATCELL"])
    with pytest.raises(ValueError):
        _assign_soil_value(df, "material", 0, ["VGPSATCELL"], 5.0)


def test_validator():
    assert not validate_scenario(sc(), nlayers=4, raise_on_error=False)[0]
    bad = sc([[0, 1], [2]])  # 2 groups, 3 entries
    assert validate_scenario(bad, nlayers=4, raise_on_error=False)[0]
    bad = sc([[0, 1], [1, 2], [3]])  # layer 1 twice
    assert validate_scenario(bad, nlayers=4, raise_on_error=False)[0]
    bad = sc([[0, 1], [2], [7]])  # layer out of range
    assert validate_scenario(bad, nlayers=4, raise_on_error=False)[0]


def test_multiplicative_with_bounds_warns():
    s = sc()
    s["per_type"] = [["multiplicative"] * 3]
    s["per_bounds"] = [[0.1, 2.0]]
    errors, warns = validate_scenario(s, nlayers=4, raise_on_error=False)
    assert any("multiplicative" in w for w in warns)
    s["per_type"] = [["None"] * 3]
    assert not any("multiplicative" in w for w in validate_scenario(s, nlayers=4, raise_on_error=False)[1])
