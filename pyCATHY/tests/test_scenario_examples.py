"""Worked examples E1 to E4 of the documentation, as tests: they build the perturbation
with the real ``perturbate`` on a stub simulation (no CATHY, no mesh) and assert which
rows of the soil table each parameter touches."""
import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.scenario import describe_ensemble

SOIL_INDEX = pd.MultiIndex.from_product([[1, 2], [0, 1]], names=["zone", "layer"])
ZONES = np.array([[1, 1, 2], [1, 2, 2]])
VEG = np.array([[1, 2, 3], [1, 2, 3]])  # class 2 spans both zones


class Sim:
    dem_parameters = {"base": 15.0}
    soil = {}

    def read_inputs(self, key, **kw):
        soil = pd.DataFrame(
            {"PERMX": 1e-4, "PERMY": 1e-4, "PERMZ": 1e-4, "POROS": 0.4}, index=SOIL_INDEX
        )
        feddes = pd.DataFrame({"ZROOT": [0.4, 0.6, 0.8]}, index=pd.Index([1, 2, 3]))
        return {"soil": (soil, feddes), "zone": (ZONES, {}), "root_map": (VEG, {})}[key]


def _scenario(name, control, nom, sd, update, bounds=None):
    listed = isinstance(nom, list)
    sc = {
        "per_name": [name], "per_type": [["None"] * len(nom) if listed else "None"],
        "per_nom": [nom], "per_mean": [nom], "per_sigma": [sd], "sampling_type": ["normal"],
        "transf_type": ["None"], "pert_control_name": [control], "listUpdateParm": update,
    }
    if bounds is not None:
        sc["per_bounds"] = [bounds]
    return sc


def _describe(sc, nens=8):
    entries = perturbate.perturbate(Sim(), sc, nens)
    return describe_ensemble(entries, sc["listUpdateParm"], SOIL_INDEX, ZONES, VEG)


def _rows(df, key):
    target = df.set_index("key").loc[key, "target"]
    inside = target.split(": ", 1)[1].split("  [")[0]
    return {tuple(int(v) for v in part.strip("() ").split(",")) for part in inside.split("), (")}


def test_E1_global_parameter_touches_every_row():
    df = _describe(_scenario("Ks", None, 1e-4, 1e-5, ["Ks"]))
    assert _rows(df, "Ks") == {(1, 0), (1, 1), (2, 0), (2, 1)}  # a scalar Ks is "Ks", not "Ks0"


def test_E2_layers():
    df = _describe(_scenario("Ks", "layers", [1e-4, 1e-5], [1e-5, 1e-6], ["St. var.", "Ks0"]))
    assert _rows(df, "Ks0") == {(1, 0), (2, 0)}  # layer 0 in every zone
    assert _rows(df, "Ks1") == {(1, 1), (2, 1)}
    assert df.set_index("key").loc["Ks0", "updated"] and not df.set_index("key").loc["Ks1", "updated"]
    assert df.attrs["warnings"] == []


def test_E3_zones_are_one_based_while_the_index_is_zero_based():
    df = _describe(_scenario("Ks", "zone", [1e-4, 1e-5], [1e-5, 1e-6], ["Ks1"]))
    assert _rows(df, "Ks0") == {(1, 0), (1, 1)}  # Ks0 -> zone 1
    assert _rows(df, "Ks1") == {(2, 0), (2, 1)}  # Ks1 -> zone 2
    assert df.attrs["warnings"] == []


def test_E4_root_map_overlap_is_reported_and_feddes_uses_row_position():
    sc = _scenario("Ks", "root_map", [1e-4, 1e-5, 1e-6], [1e-5] * 3, ["Ks2"])
    df = _describe(sc)
    assert _rows(df, "Ks0") == {(1, 0), (1, 1)}
    assert _rows(df, "Ks1") == {(1, 0), (1, 1), (2, 0), (2, 1)}  # class 2 spans both zones
    assert _rows(df, "Ks2") == {(2, 0), (2, 1)}
    # Ks1 is written over Ks0 in zone 1: the user must be told
    assert any("Ks0, Ks1" in w and "last one written wins" in w for w in df.attrs["warnings"])

    feddes = describe_ensemble(
        [{"type_parm": "ZROOT1", "mean": 0.4, "sd": 0.005, "clip_min": 0, "clip_max": 1.5}], ["ZROOT1"]
    )
    assert "Feddes table row 1" in feddes.loc[0, "target"]
    assert feddes.loc[0, "updated"]


def test_a_write_that_would_fail_in_the_run_is_shown_as_an_error():
    df = describe_ensemble(
        [{"type_parm": "Ks5", "pert_control_name": "zone", "mean": 1, "sd": 1}], [], SOIL_INDEX, ZONES, VEG
    )
    assert df.loc[0, "target"].startswith("ERROR")


def test_root_map_without_maps_is_an_error_not_a_guess():
    df = describe_ensemble(
        [{"type_parm": "Ks0", "pert_control_name": "root_map", "mean": 1, "sd": 1}], [], SOIL_INDEX
    )
    assert "needs the zone and vegetation maps" in df.loc[0, "target"]


def test_update_name_must_be_a_key_that_was_built():
    """A scalar Ks is built as 'Ks': asking to update 'Ks0' used to be silently dropped."""
    from pyCATHY.DA.scenario import ScenarioError

    with pytest.raises(ScenarioError, match=r"'Ks0'.*did you mean 'Ks'"):
        perturbate.perturbate(Sim(), _scenario("Ks", None, 1e-4, 1e-5, ["Ks0"]), 8)
    # the state update is not a parameter
    perturbate.perturbate(Sim(), _scenario("Ks", None, 1e-4, 1e-5, ["St. var.", "Ks"]), 8)


def test_B21_nodata_is_not_a_zone_in_perturbate(monkeypatch):
    """A zone map with nodata must ask for one entry per real zone, not one more."""
    zones_nodata = np.array([[-9999, 1, 1], [-9999, 2, 2]])

    class SimNodata(Sim):
        def read_inputs(self, key, **kw):
            return (zones_nodata, {}) if key == "zone" else super().read_inputs(key, **kw)

    sc = _scenario("Ks", "zone", [1e-4, 1e-5], [1e-5, 1e-6], ["Ks1"])
    entries = perturbate.perturbate(SimNodata(), sc, 8)
    assert [e["type_parm"] for e in entries] == ["Ks0", "Ks1"]  # was 3 values -> IndexError


def test_n_classes_helper():
    assert perturbate._n_classes(np.array([[1, 2], [2, 3]])) == 3
    assert perturbate._n_classes(np.array([[-9999, 1], [1, 2]])) == 2


# ---- the layered sampling path through perturbate_parm (never exercised before) ---------
def _sample(sc, nens=64):
    entries = perturbate.perturbate(Sim(), sc, nens)
    out = {}
    for e in entries:
        out = perturbate.perturbate_parm(
            out, parm=e, type_parm=e["type_parm"], mean=e["mean"], sd=e["sd"],
            sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"],
            per_type=e["per_type"], savefig=False,
        )
    return entries, out


def test_layers_control_samples_through_perturbate_parm():
    """The standard call (no nlayers keyword) used to fail with a TypeError on None."""
    sc = _scenario("Ks", "layers", [1e-4, 1e-5], [1e-5, 1e-6], ["Ks0"])
    entries, out = _sample(sc)
    assert sorted(out) == ["Ks0", "Ks1"]
    assert all(e["nlayers"] == 2 for e in entries)
    assert out["Ks0"]["nlayers"] == 2


def test_layers_without_nlayers_fails_with_a_clear_message():
    entry = {"type_parm": "Ks0", "pert_control_name": "layers", "nominal": 1e-4, "mean": 1e-4,
             "sd": 1e-5, "sampling_type": "normal", "per_type": "None", "ensemble_size": 8,
             "clip_min": None, "clip_max": None, "transf_type": None}
    with pytest.raises(ValueError, match="number of layers"):
        perturbate.perturbate_parm({}, parm=entry, type_parm="Ks0", mean=1e-4, sd=1e-5,
                                   sampling_type="normal", ensemble_size=8, per_type="None",
                                   savefig=False)


def test_DA_describes_its_own_perturbation():
    """describe_perturbation reads dict_parm_pert, the soil table and the maps of the project."""
    from types import SimpleNamespace

    from pyCATHY.DA.cathy_DA import DA

    sc = _scenario("Ks", "root_map", [1e-4, 1e-5, 1e-6], [1e-5] * 3, ["Ks2"])
    _, dict_parm_pert = _sample(sc, nens=8)
    sim = Sim()
    da = SimpleNamespace(dict_parm_pert=dict_parm_pert, read_inputs=sim.read_inputs)
    table = DA.describe_perturbation(da, ["Ks2"])
    assert table["key"].tolist() == ["Ks0", "Ks1", "Ks2"]
    assert table.set_index("key").loc["Ks2", "updated"]
    assert any("last one written wins" in w for w in table.attrs["warnings"])  # class 2 spans zones

    da.read_inputs = lambda *a, **k: (_ for _ in ()).throw(KeyError("no project"))
    DA.describe_perturbation(da, [])  # missing inputs: described without rows, no exception
