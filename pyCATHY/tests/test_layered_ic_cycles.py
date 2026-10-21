"""A layered ic must set the prior state at cycle 0 only: at later cycles the analysed
pressure head (written just before by _update_states_from_analysis) must survive."""
import numpy as np
from pyCATHY.DA.cathy_DA import DA


def stub():
    d = DA.__new__(DA)
    d.project_name = "p"
    d.calls = []
    d.map_prop_2mesh_markers = lambda *a, **k: ("nodes",)[0]
    d.update_ic = lambda **k: d.calls.append(k)
    return d


ENTRY = {"ic0": {"pert_control_name": "layers", "ini_ic0_withLayers": np.zeros((3, 4))}}


def test_cycle0_applies_prior_layers():
    d = stub()
    d._update_ic_member(0, ENTRY, 1, "ic0", False)
    assert len(d.calls) == 1


def test_later_cycles_do_not_overwrite_the_analysed_state(capsys):
    d = stub()
    d._update_ic_member(2, ENTRY, 1, "ic0", False)
    assert d.calls == []
    assert "[DA-PARAM]" in capsys.readouterr().out
