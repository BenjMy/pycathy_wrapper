"""Tests of the ES-MDA / CATHY glue with a STUB of the simulation object.

No CATHY executable is involved: the stub reproduces the methods the adapter calls and
writes a deterministic forward model through them. This checks the orchestration (what is
written, run, read, stacked and restored), not CATHY's file formats.
"""
import logging
import re

import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import esmda_cathy
from pyCATHY.DA.esmda import ESMDAConfig
from pyCATHY.DA.esmda_cathy import (
    ESMDATimeError,
    block_diag_cov,
    observation_plan,
    predict_member_blocks,
    resolve_parameter_names,
    run_esmda_on_cathy,
    stack_data,
    state_rows_at_times,
    state_times,
)

TIMES = [3600.0, 7200.0]
ATMBC_TIMES = np.array([0.0, 3600.0, 7200.0, 10800.0])
NODES = [0, 1]


# --------------------------------------------------------------------------- helpers
def test_state_rows_at_times():
    df = pd.DataFrame(np.zeros((3, 2)), index=[0.0, 3600.0, 7200.0])
    assert state_rows_at_times(df, [7200.0, 0.0]) == [2, 0]
    assert state_rows_at_times(df, [3600.0000001]) == [1]  # within tolerance
    with pytest.raises(ESMDATimeError, match="No output at t=5000"):
        state_rows_at_times(df, [5000.0])


def test_block_diag_cov_diagonal_returns_vector():
    out = block_diag_cov([np.array([1.0, 2.0]), np.diag([3.0, 4.0]), [5.0]])
    assert out.ndim == 1 and np.allclose(out, [1, 2, 3, 4, 5])


def test_block_diag_cov_dense_when_a_block_is_full():
    full = np.array([[2.0, 0.5], [0.5, 2.0]])
    out = block_diag_cov([np.array([1.0]), full])
    assert out.shape == (3, 3)
    assert np.allclose(out, [[1, 0, 0], [0, 2, 0.5], [0, 0.5, 2]])


def test_block_diag_cov_validation():
    with pytest.raises(ValueError, match="at least one"):
        block_diag_cov([])
    with pytest.raises(ValueError, match="square"):
        block_diag_cov([np.ones((2, 3))])
    with pytest.raises(ValueError, match="1-D or 2-D"):
        block_diag_cov([np.ones((2, 2, 2))])


def test_resolve_parameter_names():
    d = {"Ks": {}, "porosity": {}, "ZROOT": {}}
    assert resolve_parameter_names(d, ["St. var.", "Ks", "porosity"]) == ["Ks", "porosity"]
    assert resolve_parameter_names(d, ["all"]) == ["Ks", "porosity", "ZROOT"]
    assert resolve_parameter_names(d, "Ks") == ["Ks"]
    with pytest.raises(ValueError, match="No parameter"):
        resolve_parameter_names(d, ["St. var."])


def test_resolve_parameter_names_is_not_a_substring_match():
    """Updating 'ic' must not pick 'a_Archie' ('ic' is a substring of 'Archie')."""
    d = {"ic": {}, "a_Archie": {}, "m_Archie": {}, "SATCELL_VG0": {}, "SATCELL_VG1": {}, "n_VG0": {}}
    assert resolve_parameter_names(d, ["St. var.", "ic"]) == ["ic"]
    # the root of numbered keys selects all of them, and only them
    assert resolve_parameter_names(d, ["SATCELL_VG"]) == ["SATCELL_VG0", "SATCELL_VG1"]
    # exact keys keep working, no duplicates
    assert resolve_parameter_names(d, ["n_VG0", "n_VG0", "ic"]) == ["n_VG0", "ic"]


def test_resolve_parameter_names_unknown_entry_raises():
    d = {"Ks": {}, "porosity": {}}
    with pytest.raises(ValueError, match=r"\['poro'\] are not perturbed"):
        resolve_parameter_names(d, ["Ks", "poro"])


def _dict_obs():
    return {
        3600.0: {"tensio1": {"data": np.array([1.0, 2.0]), "mesh_nodes": NODES},
                 "swc1": {"data": np.array([0.3]), "mesh_nodes": [4]}},
        7200.0: {"tensio1": {"data": np.array([3.0, 4.0]), "mesh_nodes": NODES}},
    }


def test_observation_plan_selection_and_order():
    plan = observation_plan(_dict_obs(), ["tensio"])
    assert [(t, k) for t, k, _ in plan] == [(3600.0, ["tensio1"]), (7200.0, ["tensio1"])]
    plan_all = observation_plan(_dict_obs(), ["all"])
    assert plan_all[0][1] == ["tensio1", "swc1"]


def test_observation_plan_rejects_unsupported_and_empty():
    obs = _dict_obs()
    obs[3600.0]["EM1"] = {"data": {}}
    with pytest.raises(NotImplementedError, match="EM1"):
        observation_plan(obs, ["all"])
    assert observation_plan(obs, ["tensio"])  # not selected -> fine
    with pytest.raises(ValueError, match="No observation selected"):
        observation_plan(_dict_obs(), ["ETact"])


# --------------------------------------------------------------------------- stub simulation
class StubSim:
    """Linear forward model psi(t, node) = A[t, node] * Ks + B[t, node] * porosity."""

    workdir = "/nowhere"
    project_name = "stub"
    output_dirname = "output"
    input_dirname = "input"
    dem_parameters = {}
    cathyH = {"MAXVEG": 1}
    grid3d = {"mesh3d_nodes": np.zeros((5, 3))}
    console = None

    def __init__(self, n_ens=40, seed=0, reject=None, unreadable=(), drop_time=None):
        rng = np.random.default_rng(seed)
        self.A = rng.uniform(0.5, 1.5, size=(3, 3))
        self.B = rng.uniform(0.5, 1.5, size=(3, 3))
        self.printed = [0.0, 3600.0, 7200.0]
        if drop_time is not None:
            self.printed = [t for t in self.printed if t != drop_time]
        self.reject = reject or {}  # run number -> members rejected
        self.unreadable = set(unreadable)
        self.n_ens = n_ens
        self.n_runs = 0
        self.prepare_called = False
        self.written = []  # ini_perturbation seen by apply_ensemble_updates, per call
        self.parm_updates = []
        self.member_params = {}
        self.ens_valid = None

    # ---- the interface used by the adapter
    def prepare_DA(self, dict_obs, dict_parm_pert, list_parm2update, list_assimilated_obs,
                   parallel, config=None, **kwargs):
        self.prepare_called = True
        self.dict_obs = dict_obs
        self.dict_parm_pert = dict_parm_pert
        self.NENS = self.n_ens
        self.ens_valid = list(range(self.NENS))
        return "./exe", ATMBC_TIMES, 100, [float(t) for t in dict_obs], "ini_perturbation"

    def transform_parameters(self, names, update_key=None, back=False, param=None):
        out = []
        for n in names:
            if update_key:
                param = self.dict_parm_pert[n][update_key]
            if self.dict_parm_pert[n]["transf_type"] == "log":
                out.append(np.exp(param) if back else np.log(param))
            else:
                out.append(np.asarray(param, dtype=float))
        return np.array(out)

    def _clip_updated_parameter(self, name, values):
        e = self.dict_parm_pert[name]
        return np.clip(values, e.get("clip_min"), e.get("clip_max"))

    def update_parm(self, **kwargs):
        self.parm_updates.append(kwargs)

    def apply_ensemble_updates(self, dict_parm_pert, list_parm2update, **kwargs):
        assert kwargs == {"cycle_nb": 0}  # the smoother always writes the "prior" key
        self.written.append({n: np.array(dict_parm_pert[n]["ini_perturbation"]) for n in list_parm2update})
        for j in range(self.NENS):
            self.member_params[j] = {n: dict_parm_pert[n]["ini_perturbation"][j] for n in list_parm2update}

    def _run_ensemble_hydrological_model(self, parallel, callexe):
        self.n_runs += 1

    def _check_before_analysis(self, ens_valid, threshold_rejected):
        bad = set(self.reject.get(self.n_runs, ()))
        return [(n not in ens_valid) or (n in bad) for n in range(self.NENS)]

    def read_outputs(self, filename, path=None):
        j = int(re.search(r"cathy_(\d+)", path).group(1)) - 1
        if j in self.unreadable:
            raise OSError(f"cannot read member {j + 1}")
        p = self.member_params[j]
        ks, phi = p["Ks"], p["porosity"]
        base = self.A * ks + self.B * phi  # (3 times, 3 nodes)
        keep = [i for i, t in enumerate([0.0, 3600.0, 7200.0]) if t in self.printed]
        df = pd.DataFrame(base[keep], index=self.printed)
        return df if filename == "psi" else (df, df.index)

    # ---- truth
    def observation(self, ks, phi):
        base = self.A * ks + self.B * phi
        return np.concatenate([base[1, NODES], base[2, NODES]])


def _problem(n_ens=40, seed=0, **stub_kw):
    sim = StubSim(n_ens=n_ens, seed=seed, **stub_kw)
    rng = np.random.default_rng(seed + 1)
    ks_true, phi_true = 2.0, 0.4
    truth = sim.observation(ks_true, phi_true)
    sigma = 0.01
    data = truth + sigma * rng.normal(size=truth.size)
    dict_obs = {
        3600.0: {"tensio1": {"data": data[:2], "mesh_nodes": NODES}},
        7200.0: {"tensio1": {"data": data[2:], "mesh_nodes": NODES}},
    }
    sim.stacked_data_cov = [np.full(2, sigma**2), np.full(2, sigma**2)]
    parm = {
        "Ks": {"ini_perturbation": np.exp(rng.normal(np.log(1.2), 0.5, n_ens)),
               "transf_type": "log", "clip_min": 0.05, "clip_max": 5.0},
        "porosity": {"ini_perturbation": rng.normal(0.3, 0.08, n_ens),
                     "transf_type": None, "clip_min": 0.05, "clip_max": 0.6},
    }
    return sim, dict_obs, parm, (ks_true, phi_true)


def _run(sim, dict_obs, parm, **kw):
    cfg = kw.pop("config", ESMDAConfig(n_iter=4, seed=0))
    return run_esmda_on_cathy(sim, dict_obs, ["St. var.", "Ks", "porosity"], parm, ["tensio"], config=cfg, **kw)


# --------------------------------------------------------------------------- orchestration
def test_stack_data_matches_plan_order():
    sim, dict_obs, parm, _ = _problem()
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["tensio"], True)
    plan = observation_plan(dict_obs, ["tensio"])
    obs = stack_data(sim, plan, ["tensio"])
    assert np.allclose(obs.data, np.concatenate([dict_obs[3600.0]["tensio1"]["data"],
                                                 dict_obs[7200.0]["tensio1"]["data"]]))
    assert obs.cov.shape == (4,) and np.allclose(obs.cov, 1e-4)
    assert obs.keep.all() and obs.use_log is False
    assert obs.layout == [(0, "tensio1", 2), (1, "tensio1", 2)]


def test_stack_data_bad_cov_shape():
    sim, dict_obs, parm, _ = _problem()
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["tensio"], True)
    sim.stacked_data_cov[0] = np.ones(3)
    with pytest.raises(ValueError, match="stacked_data_cov"):
        stack_data(sim, observation_plan(dict_obs, ["tensio"]), ["tensio"])


def test_predict_member_blocks_reads_each_time():
    sim, dict_obs, parm, _ = _problem()
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["tensio"], True)
    sim.member_params[0] = {"Ks": 2.0, "porosity": 0.4}
    plan = observation_plan(dict_obs, ["tensio"])
    blocks = predict_member_blocks(sim, 0, plan, [3600.0, 7200.0], time_offset=0.0)
    assert set(blocks) == {(0, "tensio1"), (1, "tensio1")}
    got = np.concatenate([blocks[(0, "tensio1")], blocks[(1, "tensio1")]])
    assert np.allclose(got, sim.observation(2.0, 0.4))
    # the requested state times (not the stamps) select the rows; time_offset shifts the output axis
    shifted = predict_member_blocks(sim, 0, plan, [4600.0, 8200.0], time_offset=1000.0)
    assert np.allclose(np.concatenate([shifted[(0, "tensio1")], shifted[(1, "tensio1")]]), got)


def test_predict_member_blocks_without_point_sensors_reads_nothing():
    sim, dict_obs, parm, _ = _problem()
    plan = [(3600.0, ["ERT"], [{"data": None, "data_format": "pygimli"}])]
    assert predict_member_blocks(sim, 0, plan, [3600.0]) == {}


def test_run_converges_and_bookkeeping():
    sim, dict_obs, parm, (ks_true, phi_true) = _problem()
    prior_ks = parm["Ks"]["ini_perturbation"].copy()
    prior_phi = parm["porosity"]["ini_perturbation"].copy()

    res = _run(sim, dict_obs, parm)

    # model run once per iteration + a final run
    assert sim.n_runs == 5 and len(sim.written) == 5
    # output times and end time written to every member before running
    assert len(sim.parm_updates) == sim.NENS
    assert all(u["TIMPRTi"] == [0.0, 3600.0, 7200.0] and u["TMAX"] == 10800.0 for u in sim.parm_updates)
    # posterior closer to the truth than the prior, misfit decreased
    assert res.misfit_history[-1] < res.misfit_history[0] / 5
    post = parm["Ks"]["update_nb1"], parm["porosity"]["update_nb1"]
    assert abs(post[0].mean() - ks_true) < abs(prior_ks.mean() - ks_true)
    assert abs(post[1].mean() - phi_true) < abs(prior_phi.mean() - phi_true)
    assert post[0].std() < prior_ks.std()
    # prior restored under its own key; posterior is physical (back-transformed) and bounded
    assert np.array_equal(parm["Ks"]["ini_perturbation"], prior_ks)
    assert np.array_equal(parm["porosity"]["ini_perturbation"], prior_phi)
    assert post[0].min() >= 0.05 and post[0].max() <= 5.0
    assert post[1].min() >= 0.05 and post[1].max() <= 0.6
    # the first write is the prior, later writes differ
    assert np.allclose(sim.written[0]["Ks"], prior_ks)
    assert not np.allclose(sim.written[1]["Ks"], prior_ks)
    # result kept on the simulation, ens_valid synchronised
    assert sim.esmda_result is res and sim.ens_valid == list(range(sim.NENS))
    assert res.param_names == ["Ks", "porosity"]


def test_run_is_reproducible_with_seed():
    a = _problem()
    b = _problem()
    ra, rb = _run(a[0], a[1], a[2]), _run(b[0], b[1], b[2])
    assert np.array_equal(ra.param, rb.param)


def test_rejected_member_is_dropped():
    sim, dict_obs, parm, _ = _problem(reject={2: [3, 9]})
    res = _run(sim, dict_obs, parm)
    assert not res.valid[3] and not res.valid[9] and res.valid.sum() == 38
    assert sim.ens_valid == [i for i in range(40) if i not in (3, 9)]
    assert np.all(np.isnan(res.prediction_history[-1][:, [3, 9]]))


def test_unreadable_member_output_is_dropped_not_fatal():
    sim, dict_obs, parm, _ = _problem(unreadable=[5])
    res = _run(sim, dict_obs, parm)
    assert not res.valid[5] and res.valid.sum() == 39


def test_missing_print_time_is_a_configuration_error():
    sim, dict_obs, parm, _ = _problem(drop_time=7200.0)
    with pytest.raises(ESMDATimeError, match="No output at t=7200"):
        _run(sim, dict_obs, parm)


def test_prior_restored_even_if_run_fails():
    sim, dict_obs, parm, _ = _problem(unreadable=range(40))
    prior_ks = parm["Ks"]["ini_perturbation"].copy()
    with pytest.raises(RuntimeError, match="no member"):
        _run(sim, dict_obs, parm)
    assert np.array_equal(parm["Ks"]["ini_perturbation"], prior_ks)


def test_unsupported_observation_fails_before_touching_the_simulation():
    sim, dict_obs, parm, _ = _problem()
    dict_obs[3600.0]["EM1"] = {"data": {}}
    with pytest.raises(NotImplementedError, match="EM1"):
        run_esmda_on_cathy(sim, dict_obs, ["Ks"], parm, ["all"], config=ESMDAConfig(n_iter=2))
    assert sim.prepare_called is False


@pytest.mark.parametrize("option", ["update_initial_conditions", "cleanup_between_iterations"])
def test_unsupported_options_raise(option):
    sim, dict_obs, parm, _ = _problem()
    with pytest.raises(NotImplementedError, match=option):
        _run(sim, dict_obs, parm, config=ESMDAConfig(n_iter=2, **{option: True}))
    assert sim.prepare_called is False


def test_no_parameter_selected_raises():
    sim, dict_obs, parm, _ = _problem()
    with pytest.raises(ValueError, match="No parameter"):
        run_esmda_on_cathy(sim, dict_obs, ["St. var."], parm, ["tensio"], config=ESMDAConfig(n_iter=2))


def test_adapter_module_does_not_import_cathy_DA():
    import ast
    import pathlib

    src = pathlib.Path(esmda_cathy.__file__).read_text()
    imported = {n.module for n in ast.walk(ast.parse(src)) if isinstance(n, ast.ImportFrom) and n.module}
    assert "pyCATHY.DA.cathy_DA" not in imported


def test_DA_exposes_run_DA_esmda():
    from pyCATHY.DA.cathy_DA import DA

    assert callable(DA.run_DA_esmda)


# --------------------------------------------------------------------------- hardening
def test_observation_plan_requires_data_and_mesh_nodes():
    obs = _dict_obs()
    del obs[7200.0]["tensio1"]["mesh_nodes"]
    with pytest.raises(ValueError, match=r"tensio1.*7200.*mesh_nodes"):
        observation_plan(obs, ["tensio"])
    obs = _dict_obs()
    del obs[3600.0]["tensio1"]["data"]
    with pytest.raises(ValueError, match="needs 'data'"):
        observation_plan(obs, ["tensio"])


def test_keyerror_in_prediction_is_not_swallowed_as_a_member_failure(monkeypatch):
    """A configuration/code error must stop the run, not silently drop every member."""
    sim, dict_obs, parm, _ = _problem()

    def boom(*args, **kwargs):
        raise KeyError("mesh3d_nodes")

    monkeypatch.setattr(esmda_cathy, "predict_member_blocks", boom)
    with pytest.raises(KeyError, match="mesh3d_nodes"):
        _run(sim, dict_obs, parm)


class GridStub(StubSim):
    """StubSim whose 3d mesh is not loaded yet, like a fresh DA object after prepare_DA."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid3d = {}
        self.events = []
        self.soil_reads = 0

    def run_processor(self, **kwargs):
        assert kwargs.get("IPRT1") == 3  # the dry run that writes output/grid3d
        self.events.append("grid_dry_run")

    def update_parm(self, **kwargs):
        self.events.append("parm")
        super().update_parm(**kwargs)


def _swc_problem(monkeypatch, sensor="swc1"):
    sim = GridStub(n_ens=40, seed=0)
    rng = np.random.default_rng(1)
    base_true = sim.A * 2.0 + sim.B * 0.4
    sigma = 0.01
    truth = 0.5 * np.concatenate([base_true[1, NODES], base_true[2, NODES]])
    data = truth + sigma * rng.normal(size=truth.size)
    dict_obs = {
        3600.0: {sensor: {"data": data[:2], "mesh_nodes": NODES}},
        7200.0: {sensor: {"data": data[2:], "mesh_nodes": NODES}},
    }
    sim.stacked_data_cov = [np.full(2, sigma**2), np.full(2, sigma**2)]
    parm = {
        "Ks": {"ini_perturbation": np.exp(rng.normal(np.log(1.2), 0.5, 40)),
               "transf_type": "log", "clip_min": 0.05, "clip_max": 5.0},
        "porosity": {"ini_perturbation": rng.normal(0.3, 0.08, 40),
                     "transf_type": None, "clip_min": 0.05, "clip_max": 0.6},
    }

    def fake_read_soil(*args, **kwargs):
        sim.soil_reads += 1
        return "SPP", None

    def fake_swc_mapper(state, nodes, grid3d, dem_parameters, spp, console):
        # the adapter passes grid3d["mesh3d_nodes"] (as DA._predict_member does): grid3d was loaded
        assert spp == "SPP" and len(grid3d) > 0
        return 0.5 * np.asarray(state[1])[nodes]

    monkeypatch.setattr(esmda_cathy.in_CT, "read_soil", fake_read_soil)
    monkeypatch.setattr(esmda_cathy.mapper, "swc_mapper", fake_swc_mapper)
    monkeypatch.setattr(
        esmda_cathy.out_CT, "read_grid3d",
        lambda path: {"mesh3d_nodes": np.zeros((5, 3)), "path": path},
    )
    return sim, dict_obs, parm


def test_swc_run_loads_grid3d_first_and_converges(monkeypatch):
    sim, dict_obs, parm = _swc_problem(monkeypatch)
    res = run_esmda_on_cathy(
        sim, dict_obs, ["Ks", "porosity"], parm, ["swc"], config=ESMDAConfig(n_iter=4, seed=0)
    )
    # one dry run, and before the member parm files are written
    assert sim.events.count("grid_dry_run") == 1
    assert sim.events.index("grid_dry_run") < sim.events.index("parm")
    assert sim.grid3d["path"].endswith("output/grid3d")
    # the soil of every member is read at every run (5 runs x 40 members)
    assert sim.soil_reads == 5 * 40
    assert res.misfit_history[-1] < res.misfit_history[0] / 5


def test_tensio_only_run_does_not_trigger_the_grid_dry_run(monkeypatch):
    sim, dict_obs, parm = _swc_problem(monkeypatch, sensor="tensio1")
    # the swc mapper must not be needed for tensio observations
    monkeypatch.setattr(esmda_cathy.mapper, "swc_mapper", lambda *a, **k: pytest.fail("swc mapper called"))
    # data above were built for swc (x0.5): the run only has to work, not fit
    run_esmda_on_cathy(
        sim, dict_obs, ["Ks", "porosity"], parm, ["tensio"], config=ESMDAConfig(n_iter=2, seed=0)
    )
    assert "grid_dry_run" not in sim.events
    assert sim.soil_reads == 0
    assert sim.grid3d == {}


# =========================================================================== ERT + conventions
from pyCATHY.DA import mapper  # noqa: E402
from pyCATHY.DA.esmda_cathy import (  # noqa: E402
    ObsStack,
    log10_cov,
    log10_obs,
    predict_ert_blocks,
    uses_ert,
)

N_Q = 3  # quadrupoles of the stub ERT


class ErtStub(StubSim):
    """StubSim + an ERT-like forward: rhoa_q(t) = 10 exp(0.3 (A Ks + B phi)[t, q]) > 0."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ert_inputs_valid = []
        self.ert_check_reports = []
        self.Archie_parms = {"porosity": 0.4}  # as after DA.set_Archie_parm

    def ert_model(self, ks, phi, time_index):
        base = self.A * ks + self.B * phi
        return 10.0 * np.exp(0.3 * base[time_index, :N_Q])

    def _ert_ensemble_inputs(self):
        self.ert_inputs_valid.append(list(self.ens_valid))
        paths = [f"/nowhere/stub/DA_Ensemble/cathy_{j + 1}" for j in self.ens_valid]
        return paths, np.ones(self.NENS, dtype=object)


def _ert_entry(rhoa, bad_q=None):
    df = pd.DataFrame({"rhoa": np.array(rhoa, dtype=float), "k": [1.0, 1.1, 0.9][: len(rhoa)], "valid": 1})
    if bad_q is not None:
        df.loc[bad_q, "rhoa"] = -5.0  # non-positive rhoa: unusable, flagged by the quality mask
    return {"data": df, "data_format": "pygimli", "forward_mesh_vtk_file": "mesh.vtk",
            "filename": "/x/survey.data", "data_err": 2, "dataErr": 2}


def _ert_problem(monkeypatch, n_ens=40, bad_q=None, stamps=(3600.0, 7200.0), truth_rows=(1, 2),
                 mixed=False, use_log=True, extra_rows=0, seed=0, **stub_kw):
    sim = ErtStub(n_ens=n_ens, seed=seed, **stub_kw)
    sim.use_log_transformed_obs = use_log
    rng = np.random.default_rng(seed + 1)
    ks_true, phi_true = 2.0, 0.4
    dict_obs, covs = {}, []
    for n, (stamp, row) in enumerate(zip(stamps, truth_rows)):
        rhoa = sim.ert_model(ks_true, phi_true, row) * (1 + 0.02 * rng.normal(size=N_Q))
        entry = {"ERT": _ert_entry(rhoa, bad_q if n == 0 else None)}
        var = (0.02 * np.abs(entry["ERT"]["data"]["rhoa"].to_numpy())) ** 2
        if mixed:
            psi = sim.observation(ks_true, phi_true)[2 * n:2 * n + 2] + 0.01 * rng.normal(size=2)
            entry = {"tensio1": {"data": psi, "mesh_nodes": NODES}, **entry}
            var = np.concatenate([np.full(2, 0.01**2), var])
        dict_obs[stamp] = entry
        covs.append(var)
    sim.stacked_data_cov = covs
    parm = {
        "Ks": {"ini_perturbation": np.exp(rng.normal(np.log(1.2), 0.5, n_ens)),
               "transf_type": "log", "clip_min": 0.05, "clip_max": 5.0},
        "porosity": {"ini_perturbation": rng.normal(0.3, 0.08, n_ens),
                     "transf_type": None, "clip_min": 0.05, "clip_max": 0.6},
    }
    calls = []

    def fake_map(dict_obs_, count_DA_cycle, project_name, Archie_parms, porosity, paths, **kw):
        row = kw["time_ass"]
        calls.append({"k": count_DA_cycle, "row": row, "n": len(paths), "kw": kw})
        ti = [0.0, 3600.0, 7200.0].index(sim.printed[row])
        cols = []
        for p in paths:
            j = int(re.search(r"cathy_(\d+)", p).group(1)) - 1
            m = sim.member_params[j]
            cols.append(sim.ert_model(m["Ks"], m["porosity"], ti))
        hx = np.array(cols).T
        if extra_rows:
            hx = np.vstack([hx, np.ones((extra_rows, hx.shape[1]))])
        return hx, None, None

    monkeypatch.setattr(esmda_cathy.mapper, "_map_ERT_parallel", fake_map)
    return sim, dict_obs, parm, (ks_true, phi_true), calls


def _run_ert(sim, dict_obs, parm, sensors=("ERT",), **kw):
    cfg = kw.pop("config", ESMDAConfig(n_iter=4, seed=0))
    return run_esmda_on_cathy(sim, dict_obs, ["Ks", "porosity"], parm, list(sensors), config=cfg, **kw)


# ---- conventions and helpers
def test_state_times_conventions():
    plan = [(0.0, ["ERT"], [{}]), (3600.0, ["ERT"], [{}])]
    atm = [0.0, 3600.0, 7200.0, 10800.0]
    assert state_times(plan, atm, "state_time") == [0.0, 3600.0]
    assert state_times(plan, atm, "window_start") == [3600.0, 7200.0]
    with pytest.raises(ValueError, match="not an atmbc time"):
        state_times([(1000.0, ["ERT"], [{}])], atm, "window_start")
    with pytest.raises(ValueError, match="last atmbc time"):
        state_times([(10800.0, ["ERT"], [{}])], atm, "window_start")
    with pytest.raises(ValueError, match="unknown"):
        state_times(plan, atm, "nope")


def test_log_helpers_match_the_library_formulas():
    from pyCATHY.DA.cathy_DA import _LOG_EPS, _log10_cov

    assert esmda_cathy.LOG_EPS == _LOG_EPS
    d = np.array([5.0, 50.0, 500.0, -3.0, 0.0])  # includes non-positive data (clipped like the library)
    var = (0.1 * np.abs(d) + 0.5) ** 2
    # the library takes a dense R; a 1-D vector of variances must give its diagonal
    assert np.allclose(log10_cov(d, var), np.diag(_log10_cov(d, np.diag(var))))
    full = np.diag(var) + 0.01
    assert np.allclose(log10_cov(d, full), _log10_cov(d, full))
    assert np.allclose(log10_obs(d), np.log10(np.clip(d, 0, None) + _LOG_EPS))
    assert np.isnan(log10_obs([np.nan])[0])  # NaN predictions stay NaN (row dropped by the driver)


def test_uses_ert():
    assert uses_ert(observation_plan(_ert_dict_obs_min(), ["ERT"]))
    assert not uses_ert(observation_plan(_dict_obs(), ["tensio"]))


def _ert_dict_obs_min():
    return {3600.0: {"ERT": _ert_entry([10.0, 20.0, 30.0])}}


def test_observation_plan_validates_ert_entries():
    obs = _ert_dict_obs_min()
    del obs[3600.0]["ERT"]["data_format"]
    with pytest.raises(ValueError, match="data_format"):
        observation_plan(obs, ["ERT"])
    # an ERT entry does not need mesh_nodes (unlike point sensors)
    assert observation_plan(_ert_dict_obs_min(), ["all"])


# ---- stacking: quality mask, log10, covariance
def test_stack_data_masks_bad_quadrupoles_and_logs(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, bad_q=1)
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["ERT"], True)
    plan = observation_plan(dict_obs, ["ERT"])
    obs = stack_data(sim, plan, ["ERT"])

    assert isinstance(obs, ObsStack) and obs.use_log
    assert obs.layout == [(0, "ERT", 3), (1, "ERT", 3)]
    assert obs.keep.tolist() == [True, False, True, True, True, True]  # bad quadrupole of survey 1
    raw = np.concatenate([dict_obs[3600.0]["ERT"]["data"]["rhoa"], dict_obs[7200.0]["ERT"]["data"]["rhoa"]])
    assert np.allclose(obs.data, log10_obs(raw[obs.keep]))
    assert obs.data.size == 5 and obs.cov.shape == (5,)
    # R propagated to log10 on the kept data (diagonal stays a vector)
    var = np.concatenate(sim.stacked_data_cov)[obs.keep]
    assert np.allclose(obs.cov, log10_cov(raw[obs.keep], var))
    # the same rows/transform for the predictions
    pred = np.arange(1.0, 7.0)[:, None] * np.ones((1, 4))
    assert np.allclose(obs.transform_prediction(pred), log10_obs(pred[obs.keep]))


def test_stack_data_without_log_or_mask(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, bad_q=1, use_log=False)
    sim.mask_bad_ert = False
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["ERT"], True)
    obs = stack_data(sim, observation_plan(dict_obs, ["ERT"]), ["ERT"])
    assert obs.keep.all() and not obs.use_log and obs.data.size == 6
    assert np.allclose(obs.cov, np.concatenate(sim.stacked_data_cov))


def test_stack_data_log_cov_can_be_switched_off(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    sim.log_transform_cov = False
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["ERT"], True)
    obs = stack_data(sim, observation_plan(dict_obs, ["ERT"]), ["ERT"])
    assert np.allclose(obs.cov, np.concatenate(sim.stacked_data_cov))  # R untouched, data still in log10


# ---- end-to-end with the stub ERT forward
def test_ert_run_maps_the_state_of_every_survey_time(monkeypatch):
    sim, dict_obs, parm, (ks_true, _), calls = _ert_problem(monkeypatch)
    prior_ks = parm["Ks"]["ini_perturbation"].copy()
    res = _run_ert(sim, dict_obs, parm)

    # the mapper is called once per survey per model run, with the sw-file row of the state time
    assert [c["row"] for c in calls] == [1, 2] * 5
    assert [c["k"] for c in calls] == [0, 1] * 5          # position of the survey in dict_obs
    assert all(c["kw"]["DA_cnb"] == c["k"] and c["kw"]["savefig"] is False for c in calls)
    # run on the members that passed the checks (all 40 here), read from sim.ens_valid
    assert all(c["n"] == 40 for c in calls) and sim.ert_inputs_valid == [list(range(40))] * 5
    # the filter works in log10: predictions are log10(rhoa), not rhoa (~10-50)
    assert np.nanmax(res.prediction_history[0]) < 3
    # convergence and bookkeeping as for the point sensors
    assert res.misfit_history[-1] < res.misfit_history[0] / 5
    post = parm["Ks"]["update_nb1"]
    assert abs(post.mean() - ks_true) < abs(prior_ks.mean() - ks_true)
    assert np.array_equal(parm["Ks"]["ini_perturbation"], prior_ks)
    assert sim.parm_updates[0]["TIMPRTi"] == [0.0, 3600.0, 7200.0]


def test_ert_run_drops_masked_quadrupoles_from_prediction_too(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, bad_q=1)
    res = _run_ert(sim, dict_obs, parm)
    assert res.prediction_history[0].shape == (5, 40)  # 6 raw quadrupoles - 1 masked
    assert res.misfit_history[-1] < res.misfit_history[0]


def test_ert_rejected_member_is_excluded_from_the_forward(monkeypatch):
    sim, dict_obs, parm, _, calls = _ert_problem(monkeypatch, reject={2: [3, 9]})
    res = _run_ert(sim, dict_obs, parm)
    assert not res.valid[3] and not res.valid[9]
    # run 1 (the second forward): the ERT forward only sees the 38 surviving members
    assert [c["n"] for c in calls] == [40, 40, 38, 38, 38, 38, 38, 38, 38, 38]
    assert sim.ert_inputs_valid[1] == [i for i in range(40) if i not in (3, 9)]


def test_mixed_ert_and_tensio_keep_the_data_order(monkeypatch):
    sim, dict_obs, parm, _, calls = _ert_problem(monkeypatch, mixed=True, use_log=False)
    sim.prepare_DA(dict_obs, parm, ["Ks"], ["all"], True)
    obs = stack_data(sim, observation_plan(dict_obs, ["all"]), ["all"])
    assert obs.layout == [(0, "tensio1", 2), (0, "ERT", 3), (1, "tensio1", 2), (1, "ERT", 3)]
    assert np.allclose(obs.data[:2], dict_obs[3600.0]["tensio1"]["data"])
    assert np.allclose(obs.data[2:5], dict_obs[3600.0]["ERT"]["data"]["rhoa"])

    sim2, dict_obs2, parm2, _, _ = _ert_problem(monkeypatch, mixed=True, use_log=False)
    res = _run_ert(sim2, dict_obs2, parm2, sensors=("all",))
    assert res.prediction_history[0].shape == (10, 40)
    assert res.misfit_history[-1] < res.misfit_history[0] / 3


def test_window_start_convention_reads_the_next_state(monkeypatch):
    # surveys ACQUIRED at 3600 s and 7200 s, stamped with the window start (0 s and 3600 s)
    sim, dict_obs, parm, _, calls = _ert_problem(
        monkeypatch, stamps=(0.0, 3600.0), truth_rows=(1, 2))
    res = _run_ert(sim, dict_obs, parm, config=ESMDAConfig(n_iter=3, seed=0, obs_time_convention="window_start"))
    assert [c["row"] for c in calls] == [1, 2] * 4            # states at 3600 s and 7200 s
    assert sim.parm_updates[0]["TIMPRTi"] == [0.0, 3600.0, 7200.0]
    assert res.misfit_history[-1] < res.misfit_history[0] / 3


def test_state_time_convention_with_the_same_stamps_reads_earlier_states(monkeypatch):
    sim, dict_obs, parm, _, calls = _ert_problem(monkeypatch, stamps=(0.0, 3600.0), truth_rows=(1, 2))
    _run_ert(sim, dict_obs, parm, config=ESMDAConfig(n_iter=2, seed=0))  # default: state_time
    assert [c["row"] for c in calls] == [0, 1] * 3            # one step earlier than the acquisition
    assert sim.parm_updates[0]["TIMPRTi"] == [0.0, 3600.0]


def test_ert_requires_parallel_and_fails_before_touching_the_simulation(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    with pytest.raises(NotImplementedError, match="parallel=True"):
        _run_ert(sim, dict_obs, parm, parallel=False)
    assert sim.prepare_called is False


def test_ert_forward_with_other_quadrupole_count_raises(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, extra_rows=1)
    with pytest.raises(ValueError, match="quadrupole sequence"):
        _run_ert(sim, dict_obs, parm)


def test_ert_missing_print_time_is_a_configuration_error(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, drop_time=7200.0)
    with pytest.raises(ESMDATimeError, match="No output at t=7200"):
        _run_ert(sim, dict_obs, parm)


def test_ert_check_reports_are_collected(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    monkeypatch.setattr(esmda_cathy.mapper, "pop_ert_reports", lambda: [{"cycle": 0, "members": {}}])
    _run_ert(sim, dict_obs, parm, config=ESMDAConfig(n_iter=2, seed=0))
    assert len(sim.ert_check_reports) == 3 * 2               # 3 model runs x 2 surveys
    assert sim.ert_check_reports[0]["stage"] == "forward"


def test_predict_ert_blocks_without_ert_does_nothing():
    sim = ErtStub(n_ens=4)
    plan = observation_plan(_dict_obs(), ["tensio"])
    assert predict_ert_blocks(sim, plan, [3600.0, 7200.0], [0, 1]) == {}


# ---- the library change: time_ass reaches the member forward, default unchanged
class _SerialPool:
    def __init__(self, *a, **k):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def map(self, func, iterable):
        return [func(x) for x in iterable]


class _Report:
    def warn(self):
        return self

    def summary(self):
        return "ok"

    def to_dict(self):
        return {}


def _library_call(monkeypatch, **kwargs):
    seen = {}

    def fake_forward(project_name, archie, poro, meta, path, **kw):
        seen[path] = kw
        pred = pd.DataFrame({"rhoa": [1.0, 2.0]})
        return pred, pd.DataFrame({"sw": [0.1, 0.2]}), None

    monkeypatch.setattr(mapper.multiprocessing, "Pool", _SerialPool)
    monkeypatch.setattr(mapper.Archie, "SW_2_ERT_ERa_DA", fake_forward)
    monkeypatch.setattr(mapper.ERT_checks, "check_ensemble_forward", lambda *a, **k: _Report())
    dict_obs = {3600.0: {"ERT": {**_ert_entry([1.0, 2.0, 3.0]), "fwdNoiseLevel": 0}}}
    mapper._map_ERT_parallel(
        dict_obs, 0, "prj", {}, np.ones(2, dtype=object),
        ["/x/cathy_1", "/x/cathy_2"], list_assimilated_obs="all", default_state="psi",
        DA_cnb=0, savefig=False, verbose=False, **kwargs,
    )
    mapper.pop_ert_reports()
    return seen


def test_mapper_forwards_time_ass_to_the_member_forward(monkeypatch):
    seen = _library_call(monkeypatch, time_ass=3)
    assert set(seen) == {"/x/cathy_1", "/x/cathy_2"}
    assert all(kw["time_ass"] == 3 for kw in seen.values())


def test_mapper_default_does_not_add_time_ass(monkeypatch):
    """Backwards compatible: without time_ass the forward gets exactly the keywords it got before."""
    seen = _library_call(monkeypatch)
    assert all("time_ass" not in kw for kw in seen.values())
    assert all({"DA_cnb", "savefig", "noise_level", "dict_ERT"} <= set(kw) for kw in seen.values())


# ---- against the REAL DA._get_data2assimilate (no stub of the data extraction)
def _real_da_shell(dict_obs, covs, **attrs):
    from pyCATHY.DA.cathy_DA import DA

    sim = DA.__new__(DA)  # no CATHY project: only what _get_data2assimilate / stack_data read
    sim.dict_obs = dict_obs
    sim.count_DA_cycle = 0
    sim.stacked_data_cov = covs
    for k, v in attrs.items():
        setattr(sim, k, v)
    return sim


def test_stack_data_agrees_with_the_real_data_extraction_for_ert():
    rhoa1, rhoa2 = [12.0, 30.0, 55.0], [14.0, 28.0, 60.0]
    dict_obs = {3600.0: {"ERT": _ert_entry(rhoa1)}, 7200.0: {"ERT": _ert_entry(rhoa2, bad_q=2)}}
    covs = [(0.05 * np.array(rhoa1)) ** 2, (0.05 * np.abs(_ert_entry(rhoa2, bad_q=2)["data"]["rhoa"].to_numpy())) ** 2]
    sim = _real_da_shell(dict_obs, covs, use_log_transformed_obs=True)
    # the values of every survey are the ones the library extracts
    for k in range(2):
        real, _ = sim._get_data2assimilate("all", time_ass=k)
        assert np.array_equal(real, esmda_cathy._entry_values("ERT", list(dict_obs.values())[k]["ERT"]))

    obs = stack_data(sim, observation_plan(dict_obs, "all"), "all")  # "all" as a string, as in the notebooks
    assert obs.layout == [(0, "ERT", 3), (1, "ERT", 3)]
    assert obs.keep.tolist() == [True, True, True, True, True, False]  # non-positive rhoa of survey 2
    assert np.allclose(obs.data, log10_obs(np.array(rhoa1 + rhoa2[:2])))
    assert obs.cov.shape == (5,)


def test_stack_data_agrees_with_the_real_data_extraction_for_point_sensors():
    dict_obs = {3600.0: {"tensio1": {"data": np.array([-1.2, -1.5]), "mesh_nodes": NODES}},
                7200.0: {"swc1": {"data": np.array([0.31]), "mesh_nodes": [4]}}}
    sim = _real_da_shell(dict_obs, [np.full(2, 1e-4), np.full(1, 1e-4)])
    for k, key in enumerate(("tensio1", "swc1")):
        real, _ = sim._get_data2assimilate("all", time_ass=k)
        assert np.array_equal(real, esmda_cathy._entry_values(key, list(dict_obs.values())[k][key]))


def test_stack_data_mixed_sensors_in_one_cycle_keep_the_key_order():
    """The library's own extraction cannot stack sensors of different lengths in one cycle
    (np.squeeze on a ragged list); the smoother extracts per sensor and keeps the key order."""
    psi = np.array([-1.2, -1.5])
    dict_obs = {3600.0: {"tensio1": {"data": psi, "mesh_nodes": NODES}, "ERT": _ert_entry([12.0, 30.0, 55.0])}}
    covs = [np.concatenate([np.full(2, 1e-4), np.full(3, 0.5)])]
    sim = _real_da_shell(dict_obs, covs)  # no log: pressure heads are negative
    obs = stack_data(sim, observation_plan(dict_obs, ["all"]), ["all"])
    assert obs.layout == [(0, "tensio1", 2), (0, "ERT", 3)]
    assert np.allclose(obs.data, [-1.2, -1.5, 12.0, 30.0, 55.0])
    assert obs.keep.all() and obs.cov.shape == (5,)


# ---- dry run
def test_dry_run_reports_the_plan_and_simulates_nothing(monkeypatch):
    sim, dict_obs, parm, _, calls = _ert_problem(monkeypatch, bad_q=1)
    prior = {n: parm[n]["ini_perturbation"].copy() for n in parm}
    cfg = ESMDAConfig(n_iter=3, alpha_schedule="geometric", seed=0)
    out = run_esmda_on_cathy(sim, dict_obs, ["Ks", "porosity"], parm, ["ERT"], config=cfg, dry_run=True)

    assert out["parameters"] == ["Ks", "porosity"] and out["n_members"] == 40 and out["n_iter"] == 3
    assert np.isclose(sum(1 / np.array(out["alphas"])), 1.0)
    assert out["surveys"] == [(3600.0, 3600.0), (7200.0, 7200.0)]       # (stamp, state time)
    assert out["print_times"] == [0.0, 3600.0, 7200.0]
    assert (out["n_data_raw"], out["n_data"], out["n_data_masked"]) == (6, 5, 1)
    assert out["log10"] is True and out["obs_time_convention"] == "state_time"
    # cost: n_iter updates + the final run
    assert out["n_model_runs"] == 4 and out["n_member_simulations"] == 4 * 40
    assert out["n_ert_forwards"] == 4 * 2 * 40                          # runs x surveys x members
    # nothing was simulated, read or modified
    assert sim.n_runs == 0 and calls == [] and sim.written == []
    assert all(np.array_equal(parm[n]["ini_perturbation"], prior[n]) for n in parm)
    assert not hasattr(sim, "esmda_result")
    assert len(sim.parm_updates) == sim.NENS                            # printing times were prepared


def test_dry_run_shows_the_state_times_of_the_chosen_convention(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, stamps=(0.0, 3600.0), truth_rows=(1, 2))
    cfg = ESMDAConfig(n_iter=2, seed=0, obs_time_convention="window_start")
    out = run_esmda_on_cathy(sim, dict_obs, ["Ks"], parm, ["ERT"], config=cfg, dry_run=True)
    assert out["surveys"] == [(0.0, 3600.0), (3600.0, 7200.0)]
    assert out["print_times"] == [0.0, 3600.0, 7200.0]


def test_dry_run_counts_point_sensors_without_ert_forwards():
    sim, dict_obs, parm, _ = _problem()
    out = run_esmda_on_cathy(
        sim, dict_obs, ["Ks", "porosity"], parm, ["tensio"], config=ESMDAConfig(n_iter=2), dry_run=True
    )
    assert out["n_ert_forwards"] == 0 and out["n_data"] == 4 and out["n_model_runs"] == 3
    assert sim.n_runs == 0


def test_dry_run_still_raises_setup_mistakes(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    with pytest.raises(ValueError, match="not perturbed parameters"):
        run_esmda_on_cathy(sim, dict_obs, ["Ks", "nope"], parm, ["ERT"], dry_run=True)
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch, drop_time=None)
    sim.stacked_data_cov = [np.ones(2), np.ones(3)]  # wrong size for the first survey
    with pytest.raises(ValueError, match="stacked_data_cov"):
        run_esmda_on_cathy(sim, dict_obs, ["Ks"], parm, ["ERT"], dry_run=True)


@pytest.mark.parametrize("dry_run", [True, False])
def test_ert_without_archie_parameters_is_refused(monkeypatch, dry_run):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    sim.Archie_parms = None
    with pytest.raises(ValueError, match="Archie parameters"):
        run_esmda_on_cathy(sim, dict_obs, ["Ks"], parm, ["ERT"], dry_run=dry_run)
    assert sim.n_runs == 0


def test_DA_run_DA_esmda_accepts_dry_run():
    import inspect

    from pyCATHY.DA.cathy_DA import DA

    assert inspect.signature(DA.run_DA_esmda).parameters["dry_run"].default is False


# ---- pre-run check of the state times (regression of the float32 atmbc read-back)
from pyCATHY.DA.esmda_cathy import check_state_times  # noqa: E402

# survey times as nc_model.build_times writes them (rounded to 3 decimals), seconds
_SURVEYS = [9240.001, 15419.999, 22020.001, 1222200.001]


def _plan_at(times):
    return [(t, ["ERT"], [{}]) for t in times]


def test_float32_read_back_of_the_atmbc_vector_is_not_a_mismatch():
    """DA.read_inputs('atmbc') replaces sim.atmbc by the float32 read-back (~1e-3 s at 1e4 s,
    ~0.06 s at 1e6 s): the check must accept it. An absolute 1e-6 s tolerance rejected the first survey."""
    exact = np.unique(np.concatenate([np.arange(0, 1230000, 3600.0), _SURVEYS]))
    read_back = exact.astype(np.float32).astype(float)
    assert np.abs(read_back - exact).max() > 1e-4                       # the read-back really differs
    assert not np.isclose(read_back, _SURVEYS[0], rtol=0, atol=1e-6).any()  # what the old check tested
    assert check_state_times(_plan_at(_SURVEYS), read_back, "state_time", _SURVEYS) == []
    assert check_state_times(_plan_at(_SURVEYS), exact, "state_time", _SURVEYS) == []


def test_a_survey_that_is_really_off_the_grid_is_reported_with_the_likely_cause():
    atm = np.arange(0, 20000, 3600.0)
    problems = check_state_times(_plan_at([5000.0]), atm, "state_time")
    assert len(problems) == 1
    assert "not an atmbc time" in problems[0] and "read_inputs('atmbc')" in problems[0]


def test_window_start_convention_check_and_acquisition_pairing():
    atm = np.arange(0, 20000, 3600.0)                                   # 0, 3600, 7200, ...
    # stamped with the window start, acquired at the next atmbc time: consistent
    assert check_state_times(_plan_at([0.0, 3600.0]), atm, "window_start", [3600.0, 7200.0]) == []
    # same stamps read with state_time: every survey is one step off its acquisition time
    wrong = check_state_times(_plan_at([0.0, 3600.0]), atm, "state_time", [3600.0, 7200.0])
    assert len(wrong) == 2 and all("obs_time_convention='state_time'" in w for w in wrong)
    # a stamp that is not an atmbc time cannot be a window start
    assert "not an atmbc time" in check_state_times(_plan_at([1000.0]), atm, "window_start")[0]


def test_state_times_tolerates_a_float32_atmbc_vector():
    exact = np.unique(np.concatenate([np.arange(0, 100000, 3600.0), [9240.001]]))
    read_back = exact.astype(np.float32).astype(float)
    # window_start: the stamp is the atmbc time before the survey, looked up in the read-back vector
    assert state_times(_plan_at([7200.0]), read_back, "window_start") == [float(read_back[np.searchsorted(read_back, 7200.0) + 1])]


# ---- console log of the run, as run_DA_sequential
@pytest.fixture
def script_like_logging():
    """Logging as in a plain script: no handler on the root nor on the pyCATHY logger
    (pytest installs its own root handlers, which would hide the stdout handler pyCATHY adds)."""
    names = ("", "pyCATHY")
    saved = {n: (logging.getLogger(n).handlers[:], logging.getLogger(n).level, logging.getLogger(n).propagate)
             for n in names}

    def clear():  # to call INSIDE the test: pytest re-adds its handlers when the test body starts
        for n in names:
            logging.getLogger(n).handlers = []

    clear()
    yield clear
    for n, (handlers, level, propagate) in saved.items():
        lg = logging.getLogger(n)
        lg.handlers, lg.level, lg.propagate = handlers, level, propagate


def _da_shell(tmp_path):
    from pyCATHY.DA.cathy_DA import DA

    sim = DA.__new__(DA)  # only what run_DA_esmda reads: no CATHY project
    sim.workdir, sim.project_name = str(tmp_path), "prj"
    return sim


def test_run_DA_esmda_writes_the_console_log_like_the_sequential_run(tmp_path, monkeypatch, script_like_logging, capsys):
    from pyCATHY.DA.runtime import _configure_logging

    def fake_run(sim, *args, **kwargs):
        script_like_logging()      # no application logging configuration, as in a plain script
        _configure_logging(False)  # what prepare_DA does at the start of a run
        print("a printed line")
        lg = logging.getLogger("pyCATHY.DA.esmda_cathy")
        lg.info("an info line")
        lg.warning("a warning line")
        return "done"

    monkeypatch.setattr(esmda_cathy, "run_esmda_on_cathy", fake_run)
    sim = _da_shell(tmp_path)
    assert sim.run_DA_esmda({}, ["Ks"], {}, ["swc"]) == "done"

    log = (tmp_path / "prj" / "console_log.txt").read_text()   # same file as run_DA_sequential
    assert "===== run started" in log and "===== run ended =====" in log
    for line in ("a printed line", "an info line", "a warning line"):
        assert line in log, f"{line!r} missing from the console log"
    assert "an info line" in capsys.readouterr().out            # and still shown on screen

    # a second run appends (as the sequential one), the first one is kept
    sim.run_DA_esmda({}, ["Ks"], {}, ["swc"])
    assert (tmp_path / "prj" / "console_log.txt").read_text().count("===== run started") == 2


def test_run_DA_esmda_logs_the_traceback_when_the_run_fails(tmp_path, monkeypatch, script_like_logging):
    def failing_run(sim, *args, **kwargs):
        print("before the failure")
        raise RuntimeError("boom")

    monkeypatch.setattr(esmda_cathy, "run_esmda_on_cathy", failing_run)
    sim = _da_shell(tmp_path)
    with pytest.raises(RuntimeError, match="boom"):
        sim.run_DA_esmda({}, ["Ks"], {}, ["swc"])
    log = (tmp_path / "prj" / "console_log.txt").read_text()
    assert "before the failure" in log and "===== run ended =====" in log  # the log is closed properly


# ---- the initial condition as a parameter
def test_global_initial_condition_is_estimated_like_any_parameter():
    sim, dict_obs, parm, _ = _problem()
    parm["ic"] = {"ini_perturbation": np.random.default_rng(0).normal(-1.5, 0.5, 40),
                  "transf_type": None, "clip_min": -5.0, "clip_max": -0.05, "pert_control_name": "global"}
    out = run_esmda_on_cathy(sim, dict_obs, ["St. var.", "ic"], parm, ["tensio"], config=ESMDAConfig(n_iter=2), dry_run=True)
    assert out["parameters"] == ["ic"]  # 'St. var.' alone would be refused, 'ic' makes it an estimation


def test_only_the_state_in_the_update_list_is_refused_with_a_clear_message():
    sim, dict_obs, parm, _ = _problem()
    with pytest.raises(ValueError, match=r"No parameter to update.*'St. var.' alone is not enough"):
        run_esmda_on_cathy(sim, dict_obs, ["St. var."], parm, ["tensio"], dry_run=True)


@pytest.mark.parametrize("dry_run", [True, False])
def test_layered_initial_condition_is_refused_instead_of_silently_ignored(dry_run):
    sim, dict_obs, parm, _ = _problem()
    parm["ic"] = {"ini_perturbation": np.zeros(40), "transf_type": None, "pert_control_name": "layers"}
    with pytest.raises(NotImplementedError, match="layered initial condition"):
        run_esmda_on_cathy(sim, dict_obs, ["ic"], parm, ["tensio"], dry_run=dry_run)
    assert sim.n_runs == 0


# ---- the recap of the run: atmbc period, assimilation windows, strategy
from pyCATHY.DA.esmda_cathy import build_run_summary, format_run_plan  # noqa: E402


def _summary_inputs(atm, stamps, n_params=1, **config_kw):
    names = [f"p{j}" for j in range(n_params)]
    return dict(
        names=names, prior_physical={n: np.linspace(1.0, 2.0, 20) for n in names}, n_members=20,
        config=ESMDAConfig(n_iter=3, seed=7, **config_kw), plan=_plan_at(stamps), times=list(stamps),
        print_times=[0.0, *stamps], atm_times=atm, n_data=9, n_data_raw=10, use_log=True,
        threshold_rejected=30, bounds_mode="reflect", parallel=True,
    )


def test_run_summary_counts_the_atmbc_times_of_every_window_even_for_a_float32_vector():
    atm = np.array([0, 1800, 3600, 5400, 7200, 9000], dtype=np.float32).astype(float) + 1e-4  # float32-like noise
    sm = build_run_summary(**_summary_inputs(atm, [3600.0, 9000.0]))
    assert sm["period"] == (float(atm[0]), float(atm[-1])) and sm["n_atmbc_times"] == 6
    # window 1 = (start, 3600]: 1800 and 3600 ; window 2 = (3600, 9000]: 5400, 7200, 9000
    assert [w[3] for w in sm["windows"]] == [2, 3]
    assert [(w[1], w[2]) for w in sm["windows"]] == [(3600.0, float(atm[0])), (9000.0, 3600.0)]
    assert sm["n_windows"] == 2 and sm["sensors"] == ["ERT"]
    assert sm["n_data_masked"] == 1 and sm["n_model_runs"] == 4 and sm["n_member_simulations"] == 80
    assert sm["prior"]["p0"] == pytest.approx((1.5, np.std(np.linspace(1, 2, 20))))
    assert (sm["threshold_rejected"], sm["param_bounds_mode"], sm["seed"], sm["localisation"]) == (30, "reflect", 7, False)


def test_format_run_plan_states_the_strategy_and_lists_every_window():
    atm = np.arange(0.0, 36001.0, 3600.0)
    stamps = [3600.0, 10800.0, 36000.0]
    sm = build_run_summary(**_summary_inputs(atm, stamps, alpha_schedule="geometric"))
    text = "\n".join(format_run_plan(sm))

    for expected in ("strategy", "3 update(s) + 1 final run", "alpha (geometric)", "sum(1/alpha) = 1.000",
                     "3 assimilation windows", "9 data (1 masked out of 10)", "log10: True",
                     "time convention  : state_time", "t = 0 -> 36000 s (0.00 -> 10.00 h, 0.42 d)",
                     "t = 0 -> 36000 s",
                     "11 atmbc times (10 atmbc windows)", "from t = start at every iteration",
                     "4 model runs x 20 members = 80 simulations", "seed 7", "bounds mode: reflect",
                     "exceeds 30", "p0 1.5"):
        assert expected in text, f"{expected!r} missing from the recap:\n{text}"
    rows = [l for l in format_run_plan(sm) if re.match(r"\s+\d+\s+\d", l)]
    assert len(rows) == 3
    # window 2: from 1 h to 3 h, 2 h long, atmbc times 7200 and 10800
    assert rows[1].split() == ["2", "3.00", "3.00", "1.00", "2.00", "2"]


def test_format_run_plan_truncates_a_long_parameter_list():
    sm = build_run_summary(**_summary_inputs(np.arange(0.0, 7201.0, 3600.0), [3600.0], n_params=15))
    assert "e+" not in format_run_plan(sm)[5]  # plain seconds, not 1.2e+06
    text = "\n".join(format_run_plan(sm))
    assert "15 estimated" in text and "p11 " in text and "p12 " not in text and "(+3 more)" in text


def test_dry_run_summary_carries_the_period_and_the_windows(monkeypatch):
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    out = run_esmda_on_cathy(sim, dict_obs, ["Ks", "porosity"], parm, ["ERT"], config=ESMDAConfig(n_iter=3), dry_run=True)
    assert out["period"] == (0.0, 10800.0) and out["n_atmbc_times"] == 4
    assert out["windows"] == [(3600.0, 3600.0, 0.0, 1), (7200.0, 7200.0, 3600.0, 1)]
    assert out["n_windows"] == 2 and out["sensors"] == ["ERT"]
    assert out["threshold_rejected"] == 100 and out["param_bounds_mode"] == "clip"
    assert set(out["prior"]) == {"Ks", "porosity"}


def test_the_recap_is_logged_before_the_first_model_run(monkeypatch, caplog):
    caplog.set_level(logging.INFO, logger="pyCATHY.DA")
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    seen = {}
    original = sim._run_ensemble_hydrological_model

    def run_and_check(parallel, callexe):
        seen.setdefault("recap_before_first_run", any("strategy" in r.getMessage() for r in caplog.records))
        return original(parallel, callexe)

    sim._run_ensemble_hydrological_model = run_and_check
    _run_ert(sim, dict_obs, parm, config=ESMDAConfig(n_iter=2, seed=0))
    assert seen["recap_before_first_run"] is True
    msgs = [r.getMessage() for r in caplog.records]
    recap = next(m for m in msgs if "strategy" in m)
    assert "assimilation windows" in recap and "simulated period" in recap and "2 update(s)" in recap
    first_iteration = next(i for i, m in enumerate(msgs) if m.startswith("ES-MDA iteration 1/2"))
    assert msgs.index(recap) < first_iteration


def test_dry_run_also_logs_the_recap(monkeypatch, caplog):
    caplog.set_level(logging.INFO, logger="pyCATHY.DA")
    sim, dict_obs, parm, _, _ = _ert_problem(monkeypatch)
    run_esmda_on_cathy(sim, dict_obs, ["Ks"], parm, ["ERT"], dry_run=True)
    assert any("strategy" in r.getMessage() and "assimilation windows" in r.getMessage() for r in caplog.records)
