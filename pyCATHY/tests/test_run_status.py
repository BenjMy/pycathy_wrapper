import json

import pytest

from pyCATHY.DA import run_status as rs


def _run(sim, name, files=(), states=None, status=None, issues=None):
    d = sim / name
    d.mkdir(parents=True)
    for f in files:
        (d / f"{name}_{f}").write_text("x")
    if states == "empty":
        (d / f"{name}_states").mkdir()
    elif states == "chunks":
        (d / f"{name}_states").mkdir()
        (d / f"{name}_states" / "cycle_0000.npz").write_text("x")
    if status:
        (d / f"{name}_run.json").write_text(json.dumps({"status": status}))
    if issues:
        (d / f"{name}_issues.log").write_text(issues)
    return d


def test_every_case_of_the_plots_log_gets_a_verdict(tmp_path):
    ok = ["dict_parm_pert.pkl", "df_obs.parquet", "states.nc"]
    _run(tmp_path, "SYN_NC_ic_enkf", ok, status="finished")
    _run(tmp_path, "SYN_NC_ic_esmda", ["dict_parm_pert.pkl", "states.nc"])               # no df_obs: as in the log
    _run(tmp_path, "SYN_NC_ic_VGP_esmda", states="empty")                                # nothing at all
    _run(tmp_path, "SYN_NC_x_pf", ["dict_parm_pert.pkl", "df_obs.parquet"])
    _run(tmp_path, "SYN_NC_y_enkf", states="chunks", status="failed")
    v = {r["run"]: r["verdict"] for r in rs.inspect_runs(str(tmp_path), "SYN_NC_*")}
    assert v == {"SYN_NC_ic_enkf": "ok", "SYN_NC_ic_esmda": "partial", "SYN_NC_ic_VGP_esmda": "no results",
                 "SYN_NC_x_pf": "no states", "SYN_NC_y_enkf": "failed"}


def test_finished_status_without_outputs_is_called_out(tmp_path):
    _run(tmp_path, "r", status="finished")
    assert "finished!" in rs.inspect_run(str(tmp_path / "r"))["verdict"]


def test_issue_counts_and_table(tmp_path):
    _run(tmp_path, "r1", ["dict_parm_pert.pkl", "df_obs.parquet", "states.nc"],
         issues="ERROR\tcollapse\tx2\t3-4\tm\nWARNING\tlog\tx1\t-\tm\n")
    rows = rs.inspect_runs(str(tmp_path))
    assert rows[0]["issues"] == {"errors": 1, "warnings": 1} and "1/1" in rs.format_table(rows)


def test_strict_exit_code(tmp_path, capsys):
    _run(tmp_path, "bad", states="empty")
    assert rs.main([str(tmp_path)]) == 0 and rs.main([str(tmp_path), "--strict"]) == 1
