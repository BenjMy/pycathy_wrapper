"""B5 / B10 / B4: diagnostics that say who was rejected, why, and what the update did."""
import json
import logging
from types import SimpleNamespace

import numpy as np
import pytest

from pyCATHY.DA.cathy_DA import DA, _check_positive_pressure_heads


def _da(tmp_path, nens=10):
    (tmp_path / "prj").mkdir()
    da = SimpleNamespace(
        workdir=str(tmp_path), project_name="prj", NENS=nens, count_DA_cycle=0,
        rejection_log=[], dict_parm_pert={},
    )
    da._record_rejection = lambda *a: DA._record_rejection(da, *a)
    da.rejection_summary = lambda: DA.rejection_summary(da)
    return da


def test_ledger_records_stage_member_and_reason(tmp_path):
    da = _da(tmp_path)
    da.count_DA_cycle = 1
    DA._record_rejection(da, "after_run", [4], "df_mbeconv['CUM.']True")
    DA._record_rejection(da, "after_analysis", [0, 6], "unfeasible root depth ZROOT0 = -5.56")
    assert DA.rejection_summary(da) == {"after_run": [5], "after_analysis": [1, 7]}  # 1-based
    saved = json.loads((tmp_path / "prj" / "prj_rejections.json").read_text())
    assert saved[0] == {"cycle": 1, "stage": "after_run", "member": 5, "reason": "df_mbeconv['CUM.']True"}
    assert len(saved) == 3


def test_error_message_names_stages_and_members(tmp_path):
    da = _da(tmp_path)
    DA._record_rejection(da, "after_run", [1, 2], "len(df_mbeconv)['TIME']0")
    DA._record_rejection(da, "after_analysis", [3], "unfeasible root depth ZROOT0 = -1")
    with pytest.raises(ValueError) as exc:
        DA._raise_if_too_many_rejected(da, [True] * 4 + [False] * 6, threshold_rejected=30)
    msg = str(exc.value)
    assert msg.startswith("% number of rejected ensemble is too high:40.0")  # old prefix kept
    assert "'after_run': [2, 3]" in msg and "'after_analysis': [4]" in msg
    assert "prj_rejections.json" in msg


def test_below_threshold_does_not_raise(tmp_path):
    DA._raise_if_too_many_rejected(_da(tmp_path), [True] + [False] * 9, threshold_rejected=30)


def test_positive_psi_reports_members_and_share(caplog):
    psi = np.full((10, 4), -1.0)
    psi[:5, 1] = 0.5          # member 2: half of the nodes
    psi[:, 3] = np.nan        # invalid member: ignored
    with caplog.at_level(logging.WARNING):
        found = _check_positive_pressure_heads(psi, "analysed")
    assert found == {2: 0.5}
    assert "analysed ensemble" in caplog.text and "[2]" in caplog.text


def test_no_positive_psi_is_silent(caplog):
    with caplog.at_level(logging.WARNING):
        assert _check_positive_pressure_heads(np.full((5, 3), -2.0)) == {}
    assert caplog.text == ""


def test_assimilation_times_are_logged(caplog):
    da = SimpleNamespace(ENS_times=[0.0, 3600.0], count_DA_cycle=1,
                         count_atmbc_cycle=1, all_atmbc_times=[0.0, 3600.0, 7200.0])
    with caplog.at_level(logging.INFO):
        DA._log_assimilation_times(da, 3600.0)
    assert "observation stamped t=3600.0" in caplog.text and "[3600.0, 7200.0]" in caplog.text
    assert "differs from the end of the model window" in caplog.text  # 3600 != 7200


def _bind(da):
    for name in ("_record_rejection", "rejection_summary", "_raise_if_too_many_rejected"):
        setattr(da, name, (lambda n: lambda *a, **k: getattr(DA, n)(da, *a, **k))(name))
    return da


def test_check_before_analysis_feeds_the_ledger(tmp_path):
    da = _bind(_da(tmp_path))
    reasons = {1: "df_mbeconv['CUM.']True", 3: None, 4: "df_mbeconv['time']0.02"}
    da._member_rejection_reason = lambda n: reasons.get(n)
    rejected = DA._check_before_analysis(da, ens_valid=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], threshold_rejected=90)
    assert [i + 1 for i, r in enumerate(rejected) if r] == [2, 5]
    assert DA.rejection_summary(da) == {"after_run": [2, 5]}
    assert {r["reason"] for r in da.rejection_log} == {"df_mbeconv['CUM.']True", "df_mbeconv['time']0.02"}


def test_check_after_analysis_records_the_negative_root_depths(tmp_path):
    """The run of the report: updated ZROOT0 negative for some members."""
    da = _bind(_da(tmp_path, nens=5))
    da.count_DA_cycle = 0
    da.dict_parm_pert = {"ZROOT0": {"update_nb1": [0.4, -0.2, 17.5, -5.56, 0.5]}}
    da._max_root_depth = lambda: 15.0
    valid = DA._check_after_analysis(da, "update_nb1", ["St. var.", "ZROOT0"])
    assert sorted(valid) == [0, 4]
    assert DA.rejection_summary(da) == {"after_analysis": [2, 3, 4]}
    reasons = {r["member"]: r["reason"] for r in da.rejection_log}
    assert "-5.56" in reasons[4] and "ZROOT0" in reasons[4]
