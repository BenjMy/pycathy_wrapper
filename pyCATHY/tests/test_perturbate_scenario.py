"""perturbate() end to end on a stub simulation (no CATHY run, no mesh)."""
import pathlib
import tomllib

import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.scenario import ScenarioError

HERE = pathlib.Path(__file__).parent


def _scenario(name="ic_ZROOT_upd_ZROOT"):
    with open(HERE / "scenario_ert_example.toml", "rb") as f:
        return tomllib.load(f)["scenario"][name]


class _Sim:
    dem_parameters = {"base": 15.0}
    soil = {}

    def read_inputs(self, key, **kw):
        idx = pd.MultiIndex.from_product([[1], [0]], names=["zone", "layer"])
        return {
            "soil": (pd.DataFrame({"PERMX": [1e-4], "POROS": [0.4]}, index=idx),
                     pd.DataFrame({"ZROOT": [0.4]}, index=pd.Index([1]))),
            "zone": (np.ones((3, 3)), {}),
            "root_map": (np.ones((3, 3)), {}),
        }[key]


def test_bounds_of_the_project_scenario_reach_the_entries():
    entries = {e["type_parm"]: e for e in perturbate.perturbate(_Sim(), _scenario(), 32)}
    assert (entries["ic"]["clip_min"], entries["ic"]["clip_max"]) == (-5.0, -0.05)
    assert (entries["ZROOT0"]["clip_min"], entries["ZROOT0"]["clip_max"]) == (0, 1.5)


def test_misspelled_parameter_is_reported_before_building():
    bad = dict(_scenario())
    bad["per_name"] = ["ic", "Zroot"]
    with pytest.raises(ScenarioError, match="did you mean 'ZROOT'"):
        perturbate.perturbate(_Sim(), bad, 32)


def test_validation_can_be_skipped():
    bad = dict(_scenario())
    bad["listUpdateParm"] = ["St. var.", "Ks0"]  # not perturbed: normally an error
    perturbate.perturbate(_Sim(), bad, 32, validate_scenario=False)
