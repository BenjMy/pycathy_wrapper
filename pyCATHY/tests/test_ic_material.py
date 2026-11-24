"""Per-group initial pressure head: ``ic`` with pert_control_name='material'.

One value per group of layers (scenario key ``material_layers``), so that the prior mean of every
group can be the truth of a layered initial condition (e.g. [-2.4, -1.3, -1.15]) instead of one
scalar for all the layers.
"""
import numpy as np
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.scenario import validate_scenario

GROUPS = [[0, 1, 2], [3, 4], [5, 6, 7]]
IC_TRUE = [-2.4, -1.3, -1.15]
NEN = 6


def sc(control="material", nom=IC_TRUE, sig=(0.5, 0.5, 0.5), bounds=None):
    n = len(nom)
    d = {
        "per_name": ["ic"], "per_nom": [list(nom)], "per_mean": [list(nom)], "per_sigma": [list(sig)],
        "per_type": [["None"] * n], "sampling_type": ["normal"], "transf_type": [["None"] * n],
        "pert_control_name": [control], "material_layers": GROUPS, "listUpdateParm": ["St. var."],
    }
    if bounds is not None:
        d["per_bounds"] = [bounds]
    return d


class Sim:
    dem_parameters = {"base": 2.0}

    def read_inputs(self, key, **kw):
        import pandas as pd
        idx = pd.MultiIndex.from_product([[1], range(8)], names=["zone", "layer"])
        soil = pd.DataFrame({"PERMX": 1e-6}, index=idx)
        return {"soil": (soil, pd.DataFrame({"ZROOT": [0.4]}, index=[1])),
                "zone": (np.ones((2, 2)), {}), "root_map": (np.ones((2, 2)), {})}[key]

    soil = dict(PCWLT=[-150.0], OMGC=[1.0], PZ=[1.0])


def build(**kw):
    return perturbate.perturbate(Sim(), sc(**kw), NEN, check_scenario=False)


def test_one_entry_per_group_with_its_own_prior_mean():
    ent = build()
    ics = [e for e in ent if e["type_parm"].startswith("ic")]
    assert [e["type_parm"] for e in ics] == ["ic0", "ic1", "ic2"]
    assert [e["mean"] for e in ics] == IC_TRUE
    assert all(e["pert_control_name"] == "material" and e["layer_groups"] == GROUPS for e in ics)


def test_ic_has_its_own_groups_when_ic_layers_is_given():
    """The layers of a layered initial condition need not be the soil materials (here 2 materials, 3 ic groups)."""
    s = sc()
    s["material_layers"] = [[0, 1, 2, 3], [4, 5, 6, 7]]
    s["ic_layers"] = GROUPS
    assert not validate_scenario(s, nlayers=8, raise_on_error=False)[0]
    ics = [e for e in perturbate.perturbate(Sim(), s, NEN, check_scenario=False) if e["type_parm"].startswith("ic")]
    assert all(e["layer_groups"] == GROUPS and e["nlayers"] == 3 for e in ics)
    s["ic_layers"] = [[0, 1], [3, 4], [5, 6, 7]]                 # layer 2 uncovered
    assert any("cover" in e for e in validate_scenario(s, nlayers=8, raise_on_error=False)[0])


def test_bounds_per_group_are_applied():
    ent = build(bounds={"min": -5.0, "max": -0.05})
    ics = [e for e in ent if e["type_parm"].startswith("ic")]
    assert all((e["clip_min"], e["clip_max"]) == (-5.0, -0.05) for e in ics)


def test_validator_accepts_material_and_still_rejects_zone_material():
    errs, _ = validate_scenario(sc(), nlayers=8, raise_on_error=False)
    assert not errs, errs
    bad = sc(control="zone_material")
    bad["zone_material_layers"] = [[1, g] for g in GROUPS]
    assert any("only supported" in e for e in validate_scenario(bad, nlayers=8, raise_on_error=False)[0])


def test_missing_material_layers_is_an_error():
    s = sc()
    del s["material_layers"]
    with pytest.raises(ValueError, match="material_layers"):
        perturbate.perturbate(Sim(), s, NEN, check_scenario=False)


def test_group_count_mismatch_is_an_error():
    with pytest.raises(ValueError, match="material_layers"):
        perturbate.perturbate(Sim(), sc(nom=IC_TRUE[:2], sig=(0.5, 0.5)), NEN, check_scenario=False)


def test_sampled_prior_is_centred_on_each_group_truth():
    np.random.seed(1)
    n = 4000
    for e, truth in zip([e for e in build() if e["type_parm"].startswith("ic")], IC_TRUE):
        out = perturbate.perturbate_parm({}, e, e["type_parm"], e["mean"], e["sd"], e["per_type"], "normal",
                                         ensemble_size=n, nlayers=e["nlayers"])
        v = out[e["type_parm"]]["ini_perturbation"]
        assert v.shape == (n,)
        assert abs(v.mean() - truth) < 0.05 and abs(v.std() - 0.5) < 0.05     # offset ~ 0 sigma


# ---- member writer ------------------------------------------------------------------------
def stub(n_members=NEN):
    d = DA.__new__(DA)
    d.project_name = "p"
    d.calls, d.mapped = [], []
    d.map_prop_2mesh_markers = lambda name, vals, **k: (d.mapped.append(np.array(vals)), np.array(vals))[1]
    d.update_ic = lambda **k: d.calls.append(k)
    return d


def dpp():
    rng = np.random.default_rng(0)
    return {f"ic{g}": {"pert_control_name": "material", "layer_groups": GROUPS,
                       "ini_perturbation": rng.normal(IC_TRUE[g], 0.5, NEN)} for g in range(3)}


def test_writer_expands_group_values_to_every_layer_of_the_group():
    d, p = stub(), dpp()
    d._update_ic_member(0, p, 2, "ic1", False)
    assert len(d.calls) == 1
    v = d.mapped[0]
    assert v.shape == (8,)
    for g, layers in enumerate(GROUPS):
        assert np.allclose(v[layers], p[f"ic{g}"]["ini_perturbation"][2])
    assert d.calls[0]["INDP"] == 1 and d.calls[0]["IPOND"] == 0


def test_writer_gives_the_same_vector_whatever_group_key_triggers_it():
    d, p = stub(), dpp()
    for k in ("ic0", "ic1", "ic2"):
        d._update_ic_member(0, p, 3, k, False)
    assert np.allclose(d.mapped[0], d.mapped[1]) and np.allclose(d.mapped[1], d.mapped[2])


def test_writer_keeps_the_analysed_state_after_cycle_0():
    d = stub()
    d._update_ic_member(3, dpp(), 0, "ic0", False)
    assert d.calls == []


def test_writer_rejects_groups_that_leave_a_layer_uncovered():
    d, p = stub(), dpp()
    for k in p:
        p[k]["layer_groups"] = [[0, 2], [3, 4], [5, 7]]          # layers 1 and 6 have no value
    with pytest.raises(ValueError, match="cover every layer"):
        d._update_ic_member(0, p, 0, "ic0", False)
