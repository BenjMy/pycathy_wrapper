import copy
import pathlib
import tomllib

import numpy as np
import pytest

from pyCATHY.DA.scenario import ScenarioError, validate_scenario

HERE = pathlib.Path(__file__).parent


def _scenario(name):
    with open(HERE / "scenario_ert_example.toml", "rb") as f:
        return tomllib.load(f)["scenario"][name]


@pytest.mark.parametrize("name", ["ic", "ic_ZROOT_NOupdate", "ic_ZROOT_upd_ZROOT"])
def test_project_scenarios_are_valid(name):
    errors, warnings = validate_scenario(_scenario(name), raise_on_error=False)
    assert errors == []


def _errs(sc, **kw):
    return validate_scenario(sc, raise_on_error=False, **kw)[0]


def test_unknown_parameter_name_suggests_the_right_one():
    sc = copy.deepcopy(_scenario("ic_ZROOT_upd_ZROOT"))
    sc["per_name"][1] = "Zroot"
    assert any("unknown parameter 'Zroot'" in e and "ZROOT" in e for e in _errs(sc))


def test_misaligned_lists():
    sc = copy.deepcopy(_scenario("ic_ZROOT_upd_ZROOT"))
    sc["per_sigma"] = [1.75]
    errs = _errs(sc)
    assert len(errs) == 1 and "per_sigma" in errs[0]


def test_update_of_a_parameter_that_is_not_perturbed():
    sc = copy.deepcopy(_scenario("ic_ZROOT_NOupdate"))
    sc["listUpdateParm"] = ["St. var.", "Ks0"]
    assert any("Ks0" in e and "not perturbed" in e for e in _errs(sc))


def test_bad_bounds_and_mean_outside_bounds():
    sc = copy.deepcopy(_scenario("ic_ZROOT_upd_ZROOT"))
    sc["per_bounds"][1] = {"mn": 0, "mx": 1.5}
    assert any("unknown bound key" in e for e in _errs(sc))
    sc["per_bounds"][1] = {"min": 0.5, "max": 1.5}  # mean 0.4 is outside
    assert any("outside the bounds" in e for e in _errs(sc))


def _layered(nom, control="layers"):
    return {
        "per_name": ["Ks"], "per_type": [["None"] * len(nom)], "per_nom": [nom],
        "per_mean": [nom], "per_sigma": [[0.1] * len(nom)], "sampling_type": ["normal"],
        "transf_type": ["None"], "pert_control_name": [control], "listUpdateParm": ["Ks0"],
    }


def test_layer_count_must_match_the_model():
    assert _errs(_layered([1e-4, 1e-5]), nlayers=2) == []
    assert any("3 layers" in e for e in _errs(_layered([1e-4, 1e-5]), nlayers=3))


def test_control_without_list_is_an_error():
    sc = _layered([1e-4])
    sc["per_nom"], sc["per_mean"], sc["per_sigma"] = [1e-4], [1e-4], [0.1]
    assert any("needs one entry per layer" in e for e in _errs(sc))


def test_zone_and_class_numbering_must_be_contiguous():
    """B20/B21: entry k is written to class k + 1; a gap or nodata silently misplaces it."""
    zones = np.array([[1, 1, 3], [1, 3, 3]])  # class 2 missing
    errs = _errs(_layered([1e-4, 1e-5], "zone"), zone_raster=zones)
    assert any("numbered 1..2" in e for e in errs)
    ok = np.array([[1, 1, 2], [1, 2, -9999]])  # nodata is not a zone
    assert _errs(_layered([1e-4, 1e-5], "zone"), zone_raster=ok) == []


def test_all_errors_are_reported_together():
    sc = copy.deepcopy(_scenario("ic_ZROOT_upd_ZROOT"))
    sc["per_name"][1] = "Zroot"
    sc["sampling_type"][0] = "gauss"
    with pytest.raises(ScenarioError) as exc:
        validate_scenario(sc)
    assert str(exc.value).count("\n  - ") >= 3  # name, sampling, and the update of ZROOT0


def test_update_index_beyond_a_listed_parameter():
    sc = _layered([1e-4, 1e-5])
    sc["listUpdateParm"] = ["Ks0", "Ks1"]
    assert _errs(sc, nlayers=2) == []
    sc["listUpdateParm"] = ["Ks2"]
    assert any("Ks2" in e and "2 entries" in e for e in _errs(sc, nlayers=2))


# ---- B22: zones are paired with their parameter by name, not by position -----------------
from pyCATHY.DA.scenario import select_zone_parameters  # noqa: E402


def test_pairing_follows_the_name_not_the_order():
    lst = ["St. var.", "ZROOT1", "ZROOT0"]  # the positional rule gave zone 0 the parameter ZROOT1
    assert select_zone_parameters(lst, 0) == ["St. var.", "ZROOT0"]
    assert select_zone_parameters(lst, 1) == ["St. var.", "ZROOT1"]


def test_several_parameters_per_zone():
    lst = ["St. var.", "Ks0", "ZROOT0", "Ks1", "ZROOT1"]
    assert select_zone_parameters(lst, 0) == ["St. var.", "Ks0", "ZROOT0"]
    assert select_zone_parameters(lst, 1) == ["St. var.", "Ks1", "ZROOT1"]


def test_zone_without_a_parameter_updates_the_state_only():
    """The positional rule raised IndexError here."""
    assert select_zone_parameters(["St. var.", "ZROOT0"], 2) == ["St. var."]


def test_state_only_and_parameter_only_lists():
    assert select_zone_parameters(["St. var."], 3) == ["St. var."]
    assert select_zone_parameters(["ZROOT0", "ZROOT1"], 1) == ["ZROOT1"]


def test_legacy_positional_rule_is_kept_for_names_without_suffix():
    assert select_zone_parameters(["St. var.", "Ks"], 0) == ["St. var.", "Ks"]
    assert select_zone_parameters(["St. var."], 0) == ["St. var."]


def test_DA_uses_the_name_pairing_only_for_zone_localisations():
    from types import SimpleNamespace

    from pyCATHY.DA.cathy_DA import DA

    lst = ["St. var.", "ZROOT1", "ZROOT0"]
    for loc, expected in (("veg_map", ["St. var.", "ZROOT0"]), ("zones", ["St. var.", "ZROOT0"]),
                          ("nodes", lst), (None, lst)):
        da = SimpleNamespace(localisation=loc)
        assert DA._zone_parameters(da, lst, 0) == expected
