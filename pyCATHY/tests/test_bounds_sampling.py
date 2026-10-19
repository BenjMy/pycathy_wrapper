"""per_bounds must reach the prior sample of every control (layers, zone, root_map), B31."""
import importlib.util
import pathlib

import numpy as np
import pytest

from pyCATHY.DA import perturbate

_spec = importlib.util.spec_from_file_location(
    "examples_helpers", pathlib.Path(__file__).with_name("test_scenario_examples.py"))
H = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(H)


def _sampled(control, nom, bounds, nens=128):
    sc = H._scenario("Ks", control, nom, [1e-5] * len(nom), ["St. var."], bounds=bounds)
    entries = perturbate.perturbate(H.Sim(), sc, nens)
    out = {}
    for e in entries:
        out = perturbate.perturbate_parm(
            out, parm=e, type_parm=e["type_parm"], mean=e["mean"], sd=e["sd"],
            sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"],
            per_type=e["per_type"], savefig=False)
    return entries, out


@pytest.mark.parametrize("control, nom", [("layers", [1e-4] * 2), ("zone", [1e-4] * 2), ("root_map", [1e-4] * 3)])
def test_bounds_are_stored_and_respected_by_the_prior(control, nom):
    """Bounds far narrower than the prior (sd 1e-5): the sample must stay inside them."""
    entries, out = _sampled(control, nom, [9e-5, 1.1e-4])
    for key, val in out.items():
        assert (val["clip_min"], val["clip_max"]) == (9e-5, 1.1e-4), f"{control} {key}: bounds dropped"
        sample = np.asarray(val["sampling"]).ravel()
        assert sample.min() >= 9e-5 and sample.max() <= 1.1e-4, f"{control} {key}"


def test_max_only_bound_is_not_ignored():
    _, out = _sampled("root_map", [1e-4] * 3, {"max": 1.05e-4})
    for val in out.values():
        assert val["clip_min"] is None and val["clip_max"] == 1.05e-4
        assert np.asarray(val["sampling"]).max() <= 1.05e-4


def test_per_entry_bounds_use_the_entry_own_bounds():
    spec = [{"min": 9.5e-5, "max": 1.05e-4}, {"min": 9e-5, "max": 1.1e-4}]
    _, out = _sampled("zone", [1e-4, 1e-4], spec)
    assert (out["Ks0"]["clip_min"], out["Ks0"]["clip_max"]) == (9.5e-5, 1.05e-4)
    assert (out["Ks1"]["clip_min"], out["Ks1"]["clip_max"]) == (9e-5, 1.1e-4)


def test_an_open_side_does_not_inherit_the_previous_entry():
    """Each entry starts from its own defaults (the loop used to feed the previous result back)."""
    spec = [{"min": 9e-5, "max": 1.1e-4}, {"max": 1.2e-4}]
    entries, _ = _sampled("layers", [1e-4, 1e-4], spec)
    assert (entries[1]["clip_min"], entries[1]["clip_max"]) == (None, 1.2e-4)


def test_unbounded_scenarios_are_unchanged():
    entries, out = _sampled("layers", [1e-4, 1e-4], None)
    assert all(e["clip_min"] is None and e["clip_max"] is None for e in entries)


def test_ic_and_zroot_keep_their_defaults_when_the_scenario_is_silent():
    sc = H._scenario("Ks", None, 1e-4, 1e-5, ["Ks"])
    assert all(e["clip_min"] is None for e in perturbate.perturbate(H.Sim(), sc, 8))


# ---- the truncated sampler ---------------------------------------------------------------
def test_truncated_sampler_accepts_open_sides_and_zero_sd():
    x = perturbate.sampling_dist_trunc(None, 1.0, 500, loc=0.0, scale=1.0)
    assert x.max() <= 1.0 and x.min() < -1.0
    y = perturbate.sampling_dist_trunc(-1.0, None, 500, loc=0.0, scale=1.0)
    assert y.min() >= -1.0 and y.max() > 1.0
    z = perturbate.sampling_dist_trunc(0.0, 1.0, 5, loc=2.0, scale=0.0)  # no spread: mean, inside
    np.testing.assert_allclose(z, 1.0)
