import numpy as np
import pytest

from pyCATHY.DA.time_alignment import check_time_alignment


def _levels(findings):
    return [f.level for f in findings]


ATMBC = [0.0, 3600.0, 7200.0, 10800.0]


def test_observations_at_window_starts_get_only_the_offset_warning():
    findings = check_time_alignment([0.0, 3600.0, 7200.0], ATMBC)
    assert _levels(findings) == ["warning"]
    assert "t + 3600 s" in findings[0].message and "3 observation(s)" in findings[0].message


def test_final_time_is_never_assimilated_and_shifts_nothing_before_it():
    findings = check_time_alignment([0.0, 3600.0, 10800.0], ATMBC)
    err = [f for f in findings if f.level == "error"][0]
    assert "1 of 3" in err.message and "final time 10800" in err.message


def test_a_rounded_time_is_detected_as_almost_equal():
    """The script comment: 'round elements to the nearest hundred ... a bug in the simulation'."""
    findings = check_time_alignment([0.0, 3600.0000001, 7200.0], ATMBC)
    err = [f for f in findings if f.level == "error"][0]
    assert "almost equal to a window start" in err.message
    assert "paired with the wrong analysis" in err.message


def test_unsorted_and_duplicated_times():
    assert "error" in _levels(check_time_alignment([3600.0, 0.0], ATMBC))
    assert any("duplicated" in f.message for f in check_time_alignment([0.0, 0.0, 3600.0], ATMBC))


def test_irregular_windows_report_the_range():
    findings = check_time_alignment([0.0, 600.0], [0.0, 600.0, 3600.0])
    msg = [f for f in findings if f.level == "warning"][0].message
    assert "t + 600 s" in msg and "up to t + 3000 s" in msg


def test_nothing_to_check():
    assert check_time_alignment([], ATMBC) == []
    assert check_time_alignment([0.0], [0.0]) == []


def test_the_pairing_rule_matches_the_run_loop():
    """Mirror of the loop in _run_DA_sequential_impl/_run_window: which observation each
    analysis gets, and which ones are skipped (the checker must agree with it)."""
    obs, atmbc = [0.0, 3599.5, 7200.0, 10800.0], ATMBC
    analyses = []  # (window start, observation used)
    for t in atmbc[:-1]:
        if t in obs:
            analyses.append((t, obs[len(analyses)]))  # data selected by counter
    assert analyses == [(0.0, 0.0), (7200.0, 3599.5)]  # 3599.5 is paired with the window at 7200
    err = [f for f in check_time_alignment(obs, atmbc) if f.level == "error"][0]
    assert "2 of 4" in err.message


def test_the_run_logs_the_findings_without_stopping(caplog):
    import logging
    from types import SimpleNamespace

    from pyCATHY.DA.cathy_DA import DA

    da = SimpleNamespace(ENS_times=[0.0, 3599.5, 10800.0], all_atmbc_times=ATMBC)
    with caplog.at_level(logging.WARNING):
        DA._log_time_alignment(da)
    assert any(r.levelno == logging.ERROR and "never be assimilated" in r.message for r in caplog.records)

    caplog.clear()
    broken = SimpleNamespace(ENS_times=None, all_atmbc_times=ATMBC)
    DA._log_time_alignment(broken)  # a diagnostic must never raise
    assert caplog.text == ""
