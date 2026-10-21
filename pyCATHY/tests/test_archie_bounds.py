import numpy as np
import pytest
from pyCATHY.DA.perturbate import Archie_pert_rules, sampling_dist


def rules(name, mean, sd, **bounds):
    parm = {"nominal": mean, "type_parm": name, **bounds}
    return Archie_pert_rules(parm, name, 400, mean, sd, None, "normal")[0]


def test_scenario_bounds_are_used():
    np.random.seed(0)
    x = rules("a_Archie", 0.5, 0.5, clip_min=0.1, clip_max=1.0)
    assert x.min() >= 0.1 and x.max() <= 1.0


def test_defaults_kept_without_bounds():
    np.random.seed(0)
    x = rules("m_Archie", 2.0, 5.0)
    assert x.min() >= 1.3 and x.max() <= 3.0
    x = rules("n_Archie", 2.0, 5.0)
    assert x.min() >= 1.3 and x.max() <= 3.0
    assert rules("rFluid_Archie", 1.0, 5.0).min() >= 0


def test_one_sided_bound_overrides_only_that_side():
    np.random.seed(0)
    x = rules("n_Archie", 2.0, 5.0, clip_min=1.0)   # max keeps the default 3
    assert x.min() >= 1.0 and x.max() <= 3.0 and x.min() < 1.3


def test_uniform_without_minmax_uses_mean_pm_sqrt3_sd():
    np.random.seed(0)
    x = sampling_dist("uniform", 2.0, 0.1, 2000)
    h = np.sqrt(3) * 0.1
    assert x.min() >= 2.0 - h and x.max() <= 2.0 + h and abs(x.std() - 0.1) < 0.01
