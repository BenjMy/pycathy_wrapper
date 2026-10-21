import numpy as np
from pyCATHY.DA.cathy_DA import DA


def _stub(mode):
    d = DA.__new__(DA)
    d.dict_parm_pert = {"n_VG0": {"clip_min": 1.0, "clip_max": 3.0}}
    d.count_DA_cycle = 2
    d.param_bounds_mode = mode
    return d


def test_clip_prints_collapse(capsys):
    out = DA._clip_updated_parameter(_stub("clip"), "n_VG0", np.array([0.1, 0.2, 0.3]))
    txt = capsys.readouterr().out
    assert "[DA-PARAM]" in txt and "CLIP" in txt and "SPREAD COLLAPSE" in txt
    assert set(out) == {1.0}


def test_reflect_prints_mode_no_collapse(capsys):
    DA._clip_updated_parameter(_stub("reflect"), "n_VG0", np.array([0.1, 0.2, 0.3]))
    txt = capsys.readouterr().out
    assert "REFLECT" in txt and "SPREAD COLLAPSE" not in txt


def test_nothing_printed_when_inside(capsys):
    DA._clip_updated_parameter(_stub("clip"), "n_VG0", np.array([1.5, 2.0, 2.5]))
    assert capsys.readouterr().out == ""


def test_layout_banner_lists_material_mapping_and_rank_warning(capsys):
    import pandas as pd
    d = DA.__new__(DA)
    d.NENS = 2
    d.param_bounds_mode = "reflect"
    d.dict_parm_pert = {"SATCELL_VG1": {"pert_control_name": "material", "layer_groups": [[0, 1], [2, 3]]}}
    cols = ["key", "control", "target", "mean", "sd", "min", "max", "updated"]
    d.describe_perturbation = lambda l: pd.DataFrame([dict.fromkeys(cols, 1)])
    DA._log_perturbation_layout(d, ["SATCELL_VG0", "SATCELL_VG1"])
    txt = capsys.readouterr().out
    assert "REFLECT" in txt and "SATCELL_VG1 -> layers [2, 3]" in txt and "rank-deficient" in txt
