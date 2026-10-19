"""Which scenario layouts really work, measured on the real ``perturbate`` + ``perturbate_parm``
with a stub simulation, and the check that the validator reports an error exactly when the real
code fails. This table is also the source of the support matrix in the documentation
(doc/content/DA_guides/perturbation.rst): if the code changes, this test fails first."""
import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.scenario import validate_scenario

SOIL_INDEX = pd.MultiIndex.from_product([[1, 2], [0, 1]], names=["zone", "layer"])
ZONES = np.array([[1, 1, 2], [1, 2, 2]])
VEG3 = np.array([[1, 2, 3], [1, 2, 3]])
NLAYERS = 2
NEN = 8


def make_sim(n_feddes_rows):
    soil = pd.DataFrame({"PERMX": 1e-4, "PERMY": 1e-4, "PERMZ": 1e-4, "POROS": 0.4}, index=SOIL_INDEX)
    feddes = pd.DataFrame({"ZROOT": [0.4] * n_feddes_rows}, index=pd.Index(range(1, n_feddes_rows + 1)))

    class Sim:
        dem_parameters = {"base": 15.0}
        # the Feddes parameters PCWLT, OMGC and PZ are read from here
        soil_dict = {k: [v] * n_feddes_rows for k, v in dict(PCWLT=-150.0, OMGC=1.0, PZ=1.0).items()}

        def read_inputs(self, key, **kw):
            return {"soil": (soil, feddes), "zone": (ZONES, {}), "root_map": (VEG3, {})}[key]

    Sim.soil = Sim.soil_dict
    return Sim()


NOMINAL = {"Ks": 1e-4, "porosity": 0.4, "n_VG": 2.0, "alpha_VG": 2.0, "thetar_VG": 0.1,
           "SATCELL_VG": 0.9, "ic": -1.0, "WTPOSITION": 10.0, "ZROOT": 0.4, "PCREF": -4.0,
           "PCWLT": -150.0, "OMGC": 1.0, "PZ": 1.0}
COUNT = {"layers": NLAYERS, "zone": 2, "root_map": 3}


def scenario(name, control, listed, n_entries=None):
    v = NOMINAL[name]
    n = n_entries or COUNT.get(control, 1)
    wrap = (lambda x: [x * (1 + 0.1 * i) for i in range(n)]) if listed else (lambda x: x)
    feddes = name in ("ZROOT", "PCREF", "PCWLT", "OMGC", "PZ")
    sc = {
        "per_name": [name], "per_nom": [wrap(v)], "per_mean": [wrap(v)],
        "per_sigma": [[abs(v) * 0.1] * n if listed else abs(v) * 0.1],
        "per_type": ["None" if (feddes or not listed) else ["None"] * n],
        "sampling_type": ["normal"], "transf_type": ["None"], "listUpdateParm": ["St. var."],
    }
    if control is not None:
        sc["pert_control_name"] = [control]
    return sc


def real_outcome(sc, sim):
    """'ok' if the real build and sampling work, else the exception name."""
    try:
        entries = perturbate.perturbate(sim, sc, NEN, validate_scenario=False)
        out = {}
        for e in entries:
            out = perturbate.perturbate_parm(
                out, parm=e, type_parm=e["type_parm"], mean=e["mean"], sd=e["sd"],
                sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"],
                per_type=e["per_type"], savefig=False)
        return "ok" if entries else "nothing built"
    except Exception as exc:  # noqa: BLE001 - the point is to classify any failure
        return type(exc).__name__


def validator_errors(sc, n_rows):
    return validate_scenario(sc, nlayers=NLAYERS, zone_raster=ZONES, veg_raster=VEG3,
                             n_feddes_rows=n_rows, raise_on_error=False)[0]


# (name, control, listed, n_feddes_rows, works?)
SOIL = ["Ks", "n_VG", "thetar_VG", "SATCELL_VG"]
CASES = (
    [(n, None, False, 3, True) for n in SOIL + ["porosity", "ic", "WTPOSITION"]]
    + [(n, c, True, 3, True) for n in SOIL for c in ("layers", "zone", "root_map")]
    + [(n, c, True, 3, True) for n in ("porosity", "ic") for c in ("layers", "zone")]
    + [(n, "root_map", True, 3, False) for n in ("porosity", "ic")]
    + [("WTPOSITION", c, True, 3, False) for c in ("layers", "zone", "root_map")]
    + [(n, None, True, 3, False) for n in SOIL + ["porosity", "ic"]]  # a list needs a control
    # Feddes: ZROOT and PCREF scalar for 1 row, one entry per row for several
    + [(n, None, False, 1, True) for n in ("ZROOT", "PCREF")]
    + [(n, None, True, 1, False) for n in ("ZROOT", "PCREF")]
    + [(n, None, True, 3, True) for n in ("ZROOT", "PCREF")]
    + [(n, None, False, 3, False) for n in ("ZROOT", "PCREF")]
    # PCWLT, OMGC, PZ: scalar only
    + [(n, None, False, k, True) for n in ("PCWLT", "OMGC", "PZ") for k in (1, 3)]
    + [(n, None, True, k, False) for n in ("PCWLT", "OMGC", "PZ") for k in (1, 3)]
)


@pytest.mark.parametrize("name, control, listed, rows, works", CASES,
                         ids=[f"{c[0]}-{c[1]}-{'list' if c[2] else 'scalar'}-feddes{c[3]}" for c in CASES])
def test_layout_is_measured_and_the_validator_agrees(name, control, listed, rows, works):
    sc = scenario(name, control, listed, n_entries=rows if name in ("ZROOT", "PCREF") and listed else None)
    outcome = real_outcome(sc, make_sim(rows))
    errors = validator_errors(sc, rows)
    if works:
        assert outcome == "ok", f"real code: {outcome}"
        assert errors == [], f"validator rejects a layout that works: {errors}"
    else:
        assert outcome != "ok", "this layout was expected to fail in the real code"
        assert errors, f"validator accepts a layout that fails with {outcome}"


def test_VGP_is_silently_ignored_by_perturbate_and_rejected_by_the_validator():
    sc = {"per_name": ["VGP"], "per_nom": [1.0], "per_mean": [1.0], "per_sigma": [0.1],
          "per_type": ["None"], "sampling_type": ["normal"], "transf_type": ["None"],
          "listUpdateParm": ["St. var."]}
    assert real_outcome(sc, make_sim(1)) == "nothing built"  # silently ignored
    assert any("VGP" in e and "commented out" in e for e in validator_errors(sc, 1))


def test_unknown_per_type_is_a_clear_error_not_an_unbound_local():
    with pytest.raises(ValueError, match="unknown per_type"):
        perturbate.perturbate_dist({"nominal": 1.0, "type_parm": "ZROOT0"}, ["None", "None"],
                                   np.ones(4), 4)
    assert perturbate.perturbate_dist({"nominal": 1.0}, "None", np.ones(4), 4).tolist() == [1.0] * 4


def test_alpha_VG_is_built_and_sampled_but_cannot_be_written_to_the_members():
    """The build/sample test above cannot see this: the failure is in the member write."""
    from pyCATHY.DA.param_tables import _SOIL_PARM_COLUMNS

    sc = scenario("n_VG", None, False)
    sc.update(per_name=["alpha_VG"])
    assert real_outcome(sc, make_sim(3)) == "ok"                 # builds and samples
    assert "alpha_VG" not in _SOIL_PARM_COLUMNS                  # ... but has no column to write to
    assert any("alpha_VG" in e and "Unknown parameter" in e for e in validator_errors(sc, 3))


def test_every_soil_parameter_the_validator_accepts_has_a_column_to_be_written_to():
    """Dispatch -> _update_soil_params_from_perturbations needs a column for each of them."""
    from pyCATHY.DA.param_tables import _SOIL_PARM_COLUMNS
    from pyCATHY.DA.scenario import KNOWN_PARAMETERS, NOT_IMPLEMENTED

    for name in SOIL + ["porosity"]:
        assert name in _SOIL_PARM_COLUMNS and name not in NOT_IMPLEMENTED, name
    soil_like = {"Ks", "porosity", "n_VG", "alpha_VG", "thetar_VG", "SATCELL_VG"} & set(KNOWN_PARAMETERS)
    assert {n for n in soil_like if n not in _SOIL_PARM_COLUMNS} <= set(NOT_IMPLEMENTED)
