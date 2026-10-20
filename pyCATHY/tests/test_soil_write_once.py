"""P14: the soil file of a member is written once, also for a partial list of layers.

Before the fix the file was written when no NaN was left in the table: a partial update of
the layers (``listUpdateParm`` with a few layers only) wrote nothing, silently, and a full
update rewrote the file at every later key (31 writes per member for 45 van Genuchten keys).
Synthetic: the real ``_update_soil_params_from_perturbations`` and ``apply_ensemble_updates``
run on a stub that records the writes.
"""
import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA.cathy_DA import DA

NENS, NL = 3, 4
IDX = pd.MultiIndex.from_product([[1], range(NL)], names=["zone", "layer"])
BASE = pd.DataFrame(
    {"VGPSATCELL": 0.4, "VGNCELL": 1.4, "VGRMCCELL": 0.05, "PERMX": 1.0}, index=IDX
)


class Stub(DA):
    def __init__(self, tmp_path):
        self.NENS = NENS
        self.tmp = tmp_path
        self.writes = []  # (member, table)
        self.verbose = False
        for i in range(NENS):
            (tmp_path / f"m{i}").mkdir()

    def _member_dir(self, ens_nb):
        return str(self.tmp / f"m{ens_nb}")

    def _update_states_from_analysis(self, *a, **k):
        pass

    def _soil_tables_to_fill(self):
        return [BASE.copy() for _ in range(NENS)]

    def _write_member_soil(self, df, ctl, shell, FP_map=None):
        self.writes.append((int(self.cwd_member()), df.copy(), ctl))

    def cwd_member(self):
        import os
        return os.path.basename(os.getcwd())[1:]


def _pert(names):
    pert = {}
    for name in names:
        for k in range(NL):
            pert[f"{name}{k}"] = {
                "pert_control_name": "layers",
                "update_nb1": np.arange(NENS) + 10.0 * k + (100.0 if name == "n_VG" else 0.0),
            }
    return pert


def _run(tmp_path, keys, names=("SATCELL_VG", "n_VG")):
    s = Stub(tmp_path)
    s.apply_ensemble_updates(_pert(names), ["St. var."] + keys, cycle_nb=1, analysis=np.zeros((1, NENS)))
    return s


def test_partial_layers_are_written_once(tmp_path):
    s = _run(tmp_path, ["SATCELL_VG0", "SATCELL_VG1"])
    assert len(s.writes) == NENS  # one write per member (was 0)
    for member, df, ctl in s.writes:
        v = df["VGPSATCELL"].to_numpy()
        assert v[0] == member and v[1] == member + 10.0  # updated layers
        assert (v[2:] == 0.4).all()  # others untouched, no NaN
        assert not df.isna().any().any() and ctl == "layers"


def test_two_families_one_write(tmp_path):
    keys = [f"SATCELL_VG{k}" for k in range(NL)] + [f"n_VG{k}" for k in range(NL)]
    s = _run(tmp_path, keys)
    assert len(s.writes) == NENS  # was 5 per member
    for member, df, _ in s.writes:
        assert np.allclose(df["VGPSATCELL"], member + 10.0 * np.arange(NL))
        assert np.allclose(df["VGNCELL"], member + 100.0 + 10.0 * np.arange(NL))
        assert (df["PERMX"] == 1.0).all() and (df["VGRMCCELL"] == 0.05).all()


def test_no_soil_key_no_write(tmp_path):
    s = Stub(tmp_path)
    s.apply_ensemble_updates({}, ["St. var."], cycle_nb=1, analysis=np.zeros((1, NENS)))
    assert s.writes == []
