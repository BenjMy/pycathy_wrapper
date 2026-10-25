"""Homogeneous DA outputs (run_outputs) and the light issue log (issues). No CATHY executable needed."""
import logging

import numpy as np
import pytest

from pyCATHY.DA import esmda_cathy
from pyCATHY.DA.analysis import AnalysisResult, DAConfig
from pyCATHY.DA.issues import IssueLog, analysis_issues, spread_ratio
from pyCATHY.DA.run_outputs import RunOutputs
from pyCATHY.DA.state_store import StateStore, load_states
from pyCATHY.tests.test_esmda_cathy import _da_shell


# ---- layout
def test_every_output_lives_in_the_project_folder_with_the_project_prefix(tmp_path):
    o = RunOutputs(str(tmp_path), "prj")
    for path in (o.states_nc, o.states_chunks, o.issues_log, o.rejections_json, o.ert_checks_json, o.run_json):
        assert path.startswith(str(tmp_path / "prj" / "prj_"))
    assert o.console_log.endswith("console_log.txt")


def test_registry_lists_only_what_exists(tmp_path):
    o = RunOutputs(str(tmp_path), "prj")
    (tmp_path / "prj").mkdir()
    (tmp_path / "prj" / "prj_issues.log").write_text("x")
    assert o.registry() == {"issues": "prj_issues.log"}


# ---- issue log
def test_repeated_issue_is_one_line_with_a_count_and_the_cycle_range(tmp_path):
    log = IssueLog(str(tmp_path / "i.log"))
    for c in range(1, 41):
        log.record("collapse", f"spread x0.0{c % 9} at cycle {c}", cycle=c)
    lines = (tmp_path / "i.log").read_text().strip().splitlines()
    assert len(lines) == 1 and "x40" in lines[0] and "1-40" in lines[0]
    assert log.summary()["by_kind"] == {"collapse": 40}


def test_the_log_stays_small_whatever_happens(tmp_path):
    log = IssueLog(str(tmp_path / "i.log"), max_distinct=5)
    for i in range(50):
        log.record("log", f"distinct problem number word{'x' * i}")
    assert len((tmp_path / "i.log").read_text().strip().splitlines()) <= 6
    assert log.summary()["dropped"] == 45


def test_package_warnings_are_captured_and_info_is_not(tmp_path):
    log = IssueLog(str(tmp_path / "i.log"), cycle_getter=lambda: 7).attach()
    lg = logging.getLogger("pyCATHY.DA.somewhere")
    lg.setLevel(logging.INFO)
    lg.info("progress, not an issue")
    lg.warning("member %d rejected", 3)
    lg.error("something broke")
    log.detach()
    text = (tmp_path / "i.log").read_text()
    assert "progress" not in text and "member 3 rejected" in text and "ERROR" in text and "\t7\t" in text
    lg.warning("after detach")                                  # no longer captured
    assert "after detach" not in (tmp_path / "i.log").read_text()


# ---- detectors
def test_collapse_of_states_and_parameters_is_detected():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(30, 64))
    p = rng.normal(size=(2, 64))
    kinds = analysis_issues(x, x * 0.05, p, np.vstack([p[0], p[1] * 0.01]), ["ks", "ic"], n_members=64)
    assert ("collapse", "ERROR") in {(k, l) for k, l, _ in kinds}
    assert any("ic" in m for _, _, m in kinds) and not any("ks" in m for _, _, m in kinds)


def test_a_healthy_analysis_reports_nothing():
    x = np.random.default_rng(1).normal(size=(30, 64))
    assert analysis_issues(x, x * 0.8, n_members=64) == []
    assert spread_ratio(x, x) == pytest.approx(1.0)


def test_pf_degeneracy_and_rank_deficiency_are_reported():
    kinds = {k for k, _, _ in analysis_issues(pf_extra={"n_eff": 1.2}, n_members=64)}
    assert kinds == {"degeneracy"}
    p = np.random.default_rng(2).normal(size=(70, 64))
    assert {k for k, _, _ in analysis_issues(param_before=p, param_after=p, n_members=64)} == {"rank"}


# ---- wiring in DA (no model run)
def test_default_run_writes_the_light_issue_log_and_no_console_log(tmp_path, monkeypatch):
    def fake_run(sim, *a, **k):
        print("heavy console output")
        logging.getLogger("pyCATHY.DA.x").warning("a real problem at member %d", 2)
        return "done"

    monkeypatch.setattr(esmda_cathy, "run_esmda_on_cathy", fake_run)
    sim = _da_shell(tmp_path)
    sim.count_DA_cycle = 0
    assert sim.run_DA_esmda({}, ["Ks"], {}, ["swc"]) == "done"
    assert not (tmp_path / "prj" / "console_log.txt").exists()
    assert "a real problem at member 2" in (tmp_path / "prj" / "prj_issues.log").read_text()


def test_a_failing_run_leaves_the_error_in_the_issue_log(tmp_path, monkeypatch):
    monkeypatch.setattr(esmda_cathy, "run_esmda_on_cathy", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    sim = _da_shell(tmp_path)
    with pytest.raises(RuntimeError):
        sim.run_DA_esmda({}, ["Ks"], {}, ["swc"])
    text = (tmp_path / "prj" / "prj_issues.log").read_text()
    assert "ERROR" in text and "boom" in text


def test_analysis_hook_feeds_the_issue_log(tmp_path):
    sim = _da_shell(tmp_path)
    (tmp_path / "prj").mkdir()
    sim.issue_log = IssueLog(str(tmp_path / "prj" / "i.log"))
    x = np.random.default_rng(3).normal(size=(20, 16))
    res = AnalysisResult(analysis=x * 0.02, analysis_param=np.zeros((16, 0)), extra={"n_eff": 1.0})
    sim._record_analysis_issues(res, x, np.zeros((0, 16)), ["St. var."], 16)
    assert set(sim.issue_log.summary()["by_kind"]) == {"collapse", "degeneracy"}


# ---- states are written for EnKF/PF whatever save_results says
def _store(n=4):
    st = StateStore(n)
    rng = np.random.default_rng(0)
    for c in range(3):
        a = rng.normal(size=(5, n))
        st.add(c, a, a, a * 1.1)
    return st


@pytest.mark.parametrize("save_results", [False, True])
def test_states_nc_is_written_without_save_results(tmp_path, save_results):
    sim = _da_shell(tmp_path)
    (tmp_path / "prj").mkdir()
    sim.config, sim._state_store = DAConfig(save_results=save_results), _store()
    sim._finalize_state_backup()
    ds = load_states(RunOutputs(str(tmp_path), "prj").states_nc)
    assert ds.sizes["cycle"] == 3 and ds.sizes["member"] == 4


def test_write_states_false_disables_the_file(tmp_path):
    sim = _da_shell(tmp_path)
    (tmp_path / "prj").mkdir()
    sim.config, sim._state_store = DAConfig(write_states=False), _store()
    sim._finalize_state_backup()
    assert not (tmp_path / "prj" / "prj_states.nc").exists()


def test_console_log_option_is_validated():
    with pytest.raises(ValueError, match="console_log"):
        DAConfig(console_log="verbose")


# ---- run manifest (<project>_run.json), written by the library for every DA type
def _manifest_shell(tmp_path, **attrs):
    sim = _da_shell(tmp_path)
    sim.NENS = 8
    sim.count_DA_cycle = 3
    sim.rejection_log = [dict(cycle=1, stage="run", member=2, reason="x"), dict(cycle=2, stage="run", member=2, reason="x"),
                         dict(cycle=2, stage="analysis", member=5, reason="y")]
    for k, v in attrs.items():
        setattr(sim, k, v)
    return sim


def test_manifest_records_options_counts_and_outputs(tmp_path):
    import json
    from pyCATHY.DA.analysis import AnalysisConfig

    sim = _manifest_shell(tmp_path, config=DAConfig(inflation=1.1, seed=1, analysis=AnalysisConfig(perturb_obs=True)))
    with sim._run_logging(sim.config) as log:
        sim._run_da_type, sim._run_params_updated = "enkf_Evensen2009", []     # what run_DA_sequential sets
        log.record("collapse", "spread x0.1", cycle=2)
    data = json.load(open(sim._outputs().run_json))
    assert data["library_status"] == "finished" and data["n_cycles"] == 3
    assert data["n_issues"] == 1 and data["n_rejected_members"] == 2          # member 2 counted once
    assert data["wall_time_s"] is not None and data["versions"]["numpy"]
    opt = data["options"]
    assert opt["perturb_obs"] is True and opt["sakov"] is True and opt["inflation"] == 1.1 and opt["seed"] == 1
    assert opt["state_only"] is True and opt["n_updated_parameters"] == 0
    assert opt["param_bounds_mode"] is None and opt["vgp_corr"] is None          # irrelevant without parameter update
    assert data["outputs"]["issues"] == "prj_issues.log" and data["outputs"]["run"] == "prj_run.json"


def test_manifest_flags_a_failed_run_and_keeps_the_original_error(tmp_path):
    import json

    sim = _manifest_shell(tmp_path, config=DAConfig())
    with pytest.raises(RuntimeError, match="boom"):
        with sim._run_logging(sim.config):
            raise RuntimeError("boom")
    data = json.load(open(sim._outputs().run_json))
    assert data["library_status"] == "failed" and data["n_issues"] == 1


def test_manifest_with_parameters_updated_is_not_state_only(tmp_path):
    sim = _manifest_shell(tmp_path, config=DAConfig(), _run_da_type="esmda", _run_params_updated=["Ks", "porosity"])
    opt = sim.run_summary()["options"]
    assert opt["state_only"] is False and opt["updated_parameters"] == ["Ks", "porosity"] and "param_bounds_mode" not in opt


def test_manifest_does_not_shadow_the_status_of_the_script_manifest(tmp_path):
    from pyCATHY.DA.run_status import _run_status
    import json

    sim = _manifest_shell(tmp_path, config=DAConfig())
    sim.write_run_summary(status="finished")
    folder = tmp_path / "prj"
    assert _run_status(str(folder)) == "finished"                                # fallback to library_status
    (folder / "prj_nb7.json").write_text(json.dumps({"status": "failed"}))
    assert _run_status(str(folder)) == "failed"                                  # the script's own status wins
