"""ES-MDA leaves the same kind of files as the sequential DA (states, params, backup, issues). No CATHY needed."""
import os

import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import esmda_cathy
from pyCATHY.DA.analysis import DAConfig
from pyCATHY.DA.esmda import ESMDAResult
from pyCATHY.DA.issues import IssueLog, esmda_issues
from pyCATHY.DA.run_outputs import RunOutputs
from pyCATHY.DA.state_store import load_states
from pyCATHY.tests.test_esmda_cathy import _da_shell

N_ENS, N_NODES = 4, 6
TIMES = [3600.0, 7200.0]            # absolute state times
ATM = [0.0, 3600.0, 7200.0]


def _member_of(path):
    return int(path.split("cathy_")[1].split(os.sep)[0])


def _fake_sim(tmp_path, unreadable=()):
    sim = _da_shell(tmp_path)
    sim.NENS, sim.output_dirname, sim.dict_obs, sim.dict_parm_pert = N_ENS, "output", {}, {"Ks": {"x": 1}}
    sim.count_DA_cycle = 0
    sim.config = DAConfig()

    def read_outputs(filename, path):
        m = _member_of(path)
        if m in unreadable:
            raise OSError("missing psi")
        idx = [0.0, 3600.0, 7200.0]
        base = np.array([[-m - 0.01 * n - 0.001 * t for n in range(N_NODES)] for t in range(3)])
        if filename == "psi":
            return pd.DataFrame(base, index=idx)
        return pd.DataFrame(-base / 100, index=idx), None

    sim.read_outputs = read_outputs
    sim.backup_results_DA = lambda meta: setattr(sim, "_backup_meta", meta)
    return sim


def _result(valid=(1, 1, 1, 1), misfit=(8.0, 5.0, 3.0, 1.2)):
    prior = np.vstack([np.linspace(1, 2, N_ENS), np.linspace(0, 1, N_ENS)])
    hist = [prior, prior * 0.9, prior * 0.8, prior * 0.7]
    pred = [np.random.default_rng(i).normal(size=(5, N_ENS)) for i in range(4)]
    return ESMDAResult(param=hist[-1], valid=np.array(valid, bool), alphas=np.array([3.0, 3.0, 3.0]),
                       param_history=hist, prediction_history=pred, misfit_history=np.array(misfit),
                       valid_history=[np.array(valid, bool)] * 3, param_names=["Ks", "ic"])


def _phys(p):
    return {"Ks": 10 ** p[0], "ic": p[1]}


def test_states_of_the_posterior_run_in_the_common_store(tmp_path):
    sim = _fake_sim(tmp_path)
    st = esmda_cathy.esmda_states_store(sim, _result(valid=(1, 1, 0, 1)), TIMES, 0.0, ATM, DAConfig())
    cycles, members, nodes, psi = st.arrays("psi")
    assert list(cycles) == [1, 2] and psi.shape == (2, N_ENS, N_NODES)
    assert np.isnan(psi[:, 2]).all() and not np.isnan(psi[:, [0, 1, 3]]).any()      # invalid member -> NaN
    assert psi[0, 1, 0] == pytest.approx(-2 - 0.001 * 1, abs=1e-5)                   # member 2, node 0, t = 3600
    assert st.attrs["da_type"] == "esmda" and st.attrs["states_source"] == "posterior_forward_run"
    assert np.isnan(st.arrays("analysis")[3]).all()                                   # a smoother has no analysis


def test_an_unreadable_member_is_flagged_not_fatal(tmp_path, caplog):
    sim = _fake_sim(tmp_path, unreadable=(3,))
    st = esmda_cathy.esmda_states_store(sim, _result(), TIMES, 0.0, ATM, DAConfig())
    assert st.rejected[0].tolist() == [0, 0, 1, 0] and "member 3" in caplog.text


def test_params_dataset_has_every_iterate_in_physical_space(tmp_path):
    ds = esmda_cathy.esmda_params_dataset(_result(), _phys, ["Ks", "ic"])
    assert ds["Ks"].dims == ("iteration", "member") and ds.sizes["iteration"] == 4
    assert ds["prediction"].dims == ("run", "obs", "member") and ds["misfit"].size == 4
    assert float(ds["Ks"][0, 0]) == pytest.approx(10.0)


def test_record_outputs_writes_everything_a_sequential_run_leaves(tmp_path):
    sim = _fake_sim(tmp_path)
    sim.issue_log = IssueLog(RunOutputs(str(tmp_path), "prj").issues_log)
    esmda_cathy.record_esmda_outputs(sim, _result(misfit=(9.0, 7.0, 5.0, 4.0)), ["Ks", "ic"], _phys, TIMES, 0.0, ATM,
                                     ["swc"], sim.config)
    o = RunOutputs(str(tmp_path), "prj")
    assert load_states(o.states_nc).attrs["da_type"] == "esmda"
    import xarray as xr
    assert xr.load_dataset(o.params_nc).attrs["da_type"] == "esmda"
    assert sim._backup_meta["da_type"] == "esmda" and sim._backup_meta["listUpdatedparm"] == ["Ks", "ic"]
    assert "misfit" in sim.issue_log.summary()["by_kind"]                          # 4.0 > 2: data not fitted
    assert set(o.registry()) >= {"states", "params", "issues"}


def test_a_failing_part_never_stops_the_others(tmp_path):
    sim = _fake_sim(tmp_path, unreadable=(1, 2, 3, 4))                             # no state can be read
    esmda_cathy.record_esmda_outputs(sim, _result(), ["Ks", "ic"], _phys, TIMES, 0.0, ATM, ["swc"], sim.config)
    o = RunOutputs(str(tmp_path), "prj")
    assert not os.path.exists(o.states_nc) and os.path.exists(o.params_nc) and hasattr(sim, "_backup_meta")


def test_write_states_false_is_honoured(tmp_path):
    sim = _fake_sim(tmp_path)
    cfg = DAConfig(write_states=False)
    esmda_cathy.record_esmda_outputs(sim, _result(), ["Ks"], _phys, TIMES, 0.0, ATM, ["swc"], cfg)
    assert not os.path.exists(RunOutputs(str(tmp_path), "prj").states_nc)


def test_esmda_issues():
    r = _result(valid=(1, 1, 1, 0), misfit=(2.0, 3.0, 4.0, 5.0))
    kinds = {k for _, k, _, _ in esmda_issues(r.param_history, r.valid, r.misfit_history, ["Ks", "ic"])}
    assert kinds == {"misfit", "rejection"}
    collapsed = [r.param_history[0], r.param_history[0] * 0.0 + r.param_history[0].mean(axis=1, keepdims=True)]
    assert any(k == "collapse" for _, k, _, _ in esmda_issues(collapsed, np.ones(N_ENS, bool), None, ["Ks", "ic"]))


def test_load_backup_reads_states_params_and_says_clearly_when_missing(tmp_path, capsys, monkeypatch):
    import builtins

    # cathy_tools prints with rich; test_enkf.py stubs rich.console globally: use the builtin to stay order independent
    monkeypatch.setattr("pyCATHY.cathy_tools.print", builtins.print, raising=False)
    sim = _fake_sim(tmp_path)
    esmda_cathy.record_esmda_outputs(sim, _result(), ["Ks"], _phys, TIMES, 0.0, ATM, ["swc"], sim.config)
    out = sim.load_backup()
    assert "states" in out and "params" in out
    capsys.readouterr()
    os.remove(RunOutputs(str(tmp_path), "prj").states_nc)
    os.makedirs(RunOutputs(str(tmp_path), "prj").states_chunks)                      # empty folder, as in the plots log
    sim.load_backup()
    out = capsys.readouterr().out
    assert out.count("\n") == 1 and "missing:" in out and "states" in out        # ONE line, not a stack of errors


# ---- chunks (<project>_states/) and the kind of states, same rules as the sequential DA
def test_esmda_writes_chunks_only_with_save_results_and_they_reload_like_the_netcdf(tmp_path):
    from pyCATHY.DA.state_store import StateStore

    o = RunOutputs(str(tmp_path), "prj")
    sim = _fake_sim(tmp_path)
    esmda_cathy.esmda_states_store(sim, _result(), TIMES, 0.0, ATM, DAConfig())             # save_results False
    assert not os.path.isdir(o.states_chunks)

    st = esmda_cathy.esmda_states_store(sim, _result(), TIMES, 0.0, ATM, DAConfig(save_results=True))
    assert sorted(os.listdir(o.states_chunks)) == ["cycle_00001.npz", "cycle_00002.npz"]
    back = StateStore.from_chunks(o.states_chunks)
    assert back.attrs["da_type"] == "esmda" and back.attrs["states_kind"] == "smoothed_posterior"   # attrs survive
    a, b = st.to_xarray(), back.to_xarray()
    assert a.attrs["states_kind"] == b.attrs["states_kind"] == "smoothed_posterior"
    np.testing.assert_array_equal(a["psi"].values, b["psi"].values)
    assert "ALL the data" in a.attrs["description"]


def test_stale_chunks_of_a_previous_run_are_removed(tmp_path):
    o = RunOutputs(str(tmp_path), "prj")
    os.makedirs(o.states_chunks)
    open(os.path.join(o.states_chunks, "cycle_00099.npz"), "wb").close()
    sim = _fake_sim(tmp_path)
    esmda_cathy.esmda_states_store(sim, _result(), TIMES, 0.0, ATM, DAConfig(save_results=True))
    assert "cycle_00099.npz" not in os.listdir(o.states_chunks)


def test_record_outputs_with_save_results_leaves_chunks_and_netcdf(tmp_path):
    sim = _fake_sim(tmp_path)
    cfg = DAConfig(save_results=True)
    sim.config = cfg
    esmda_cathy.record_esmda_outputs(sim, _result(), ["Ks"], _phys, TIMES, 0.0, ATM, ["swc"], cfg)
    o = RunOutputs(str(tmp_path), "prj")
    assert os.path.exists(o.states_nc) and len(os.listdir(o.states_chunks)) == 2
    assert load_states(o.states_nc).attrs["states_kind"] == "smoothed_posterior"


def test_sequential_store_is_tagged_as_forecast(tmp_path):
    from pyCATHY.DA.state_store import StateStore

    st = StateStore(2)
    st.add(1, np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((3, 2)))
    assert st.to_xarray().attrs["states_kind"] == "forecast_before_update"                  # default of the sequential path


def test_states_without_the_final_forward_run_are_not_called_posterior(tmp_path):
    sim = _fake_sim(tmp_path)
    r = _result()
    r.prediction_history = r.prediction_history[:3]                  # final_forward=False: one run per iteration only
    st = esmda_cathy.esmda_states_store(sim, r, TIMES, 0.0, ATM, DAConfig())
    ds = st.to_xarray()
    assert ds.attrs["states_kind"] == "last_iteration_forward_run" and "NOT the posterior" in ds.attrs["description"]
    assert ds.attrs["states_source"] == "last_iteration_forward_run"
    assert esmda_cathy.esmda_states_store(sim, _result(), TIMES, 0.0, ATM, DAConfig()).to_xarray().attrs[
        "states_kind"] == "smoothed_posterior"
