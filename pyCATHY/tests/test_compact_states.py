"""Compact storage of the DA states (StateStore) behind ``DA.df_DA``: same content as the legacy long table,
a fraction of the memory, per-cycle chunks, node subsets, and a backup that no longer loses attributes."""
import os
import types

import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA.analysis import DAConfig
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.state_store import StateStore

NN, NE, NC = 30, 5, 6


def make_da(compact, tmp_path=None, **cfg):
    da = DA.__new__(DA)                      # no CATHY project needed: only the bookkeeping is exercised
    da.NENS = NE
    da.count_atmbc_cycle = 0
    da.workdir = str(tmp_path or ".")
    da.project_name = "prj"
    da._state_store = None
    da._df_DA_legacy = pd.DataFrame()
    da.config = DAConfig(compact_states=compact, **cfg)
    return da


def run_cycles(da, seed=0, analysis_at=(2, 5)):
    rng = np.random.RandomState(seed)
    truth = []
    for c in range(1, NC + 1):
        psi = -1.0 - rng.rand(NN, NE)
        sw = 0.4 + 0.1 * rng.rand(NN, NE)
        ana = psi + (0.05 if c in analysis_at else 0.0)
        rej = [0, 0, 1, 0, 0] if c == 3 else [0] * NE
        da.count_atmbc_cycle = c
        da._DA_df(state=[psi, sw], state_analysis=ana, rejected_ens=rej)
        truth.append((psi, sw, ana, rej))
    return truth


def test_same_content_as_legacy():
    legacy, compact = make_da(False), make_da(True)
    run_cycles(legacy)
    run_cycles(compact)
    a, b = legacy.df_DA, compact.df_DA
    assert len(a) == len(b) == NC * NE * NN
    for col in ("time", "Ensemble_nb", "psi_bef_update", "sw_bef_update_", "analysis", "OL", "rejected"):
        assert np.allclose(a[col].to_numpy(float), b[col].to_numpy(float), atol=1e-6), col
    assert (b["idnode"].to_numpy()[:NN] == np.arange(NN)).all()      # node id is explicit now


def test_memory_and_dtypes():
    legacy, compact = make_da(False), make_da(True)
    run_cycles(legacy)
    run_cycles(compact)
    assert compact._state_store.nbytes < 0.3 * legacy.df_DA.memory_usage(deep=True).sum()
    assert compact.df_DA["psi_bef_update"].dtype == np.float32 and compact.df_DA["time"].dtype == np.int32


def test_arrays_nan_for_rejected():
    da = make_da(True)
    truth = run_cycles(da)
    cyc, mem, nodes, sw = da._state_store.arrays("sw")
    assert sw.shape == (NC, NE, NN) and list(cyc) == list(range(1, NC + 1))
    assert np.isnan(sw[2, 2]).all() and not np.isnan(sw[2, 0]).any()          # member 2 rejected at cycle 3
    assert np.allclose(sw[0], truth[0][1].T, atol=1e-6)


def test_node_subset_and_analysis_only():
    da = make_da(True, store_full="analysis", store_nodes=[1, 7, 20])
    run_cycles(da, analysis_at=(2, 5))
    st = da._state_store
    assert [len(i) for i in st.node_ids] == [3, NN, 3, 3, NN, 3]              # full field at the analysis cycles only
    assert st.nbytes < NC * NE * NN * 12                                      # less than "everything at every cycle"
    _, _, nodes, a = st.arrays("psi")
    assert list(nodes) == [1, 7, 20] and a.shape == (NC, NE, 3)
    with pytest.raises(ValueError):
        DAConfig(store_full="analysis")                                       # needs store_nodes
    with pytest.raises(ValueError):
        DAConfig(store_full="sometimes")


def test_chunks_survive_and_reload(tmp_path):
    da = make_da(True, tmp_path, save_results=True)
    os.makedirs(os.path.join(tmp_path, "prj"))
    truth = run_cycles(da)
    chunk_dir = os.path.join(tmp_path, "prj", "prj_states")
    assert len(os.listdir(chunk_dir)) == NC                                    # one file per cycle, written as it goes
    back = StateStore.from_chunks(chunk_dir)
    assert np.allclose(back.arrays("sw")[3], da._state_store.arrays("sw")[3], equal_nan=True)
    da._finalize_state_backup()
    pkl = os.path.join(tmp_path, "prj", "prj_df_DA.pkl")
    nc = os.path.join(tmp_path, "prj", "prj_states.nc")
    assert os.path.isfile(nc) and not os.path.exists(pkl)                       # NetCDF by default, no legacy pickle
    import xarray as xr
    assert xr.load_dataset(nc)["sw"].shape == (NC, NE, NN)
    da.config.legacy_df_DA = True                                               # opt-in (deprecated)
    da._finalize_state_backup()
    assert os.path.isfile(pkl) and len(pd.read_pickle(pkl)) == NC * NE * NN
    # a new run of the same project removes stale chunks
    da2 = make_da(True, tmp_path, save_results=True)
    da2.count_atmbc_cycle = 1
    da2._get_state_store()
    assert not os.listdir(chunk_dir)


def test_setter_keeps_legacy_assignment():
    da = make_da(True)
    run_cycles(da)
    da.df_DA = pd.DataFrame({"time": [1.0]})                                   # load_backup style assignment
    assert len(da.df_DA) == 1 and da._state_store is None


def test_backup_saves_every_attribute_even_if_one_fails(tmp_path, monkeypatch):
    from pyCATHY import cathy_tools
    da = make_da(True, tmp_path)
    os.makedirs(os.path.join(tmp_path, "prj"))
    da.dict_parm_pert = {"a": 1}
    da.df_Archie = pd.DataFrame({"x": [1.0]})
    da.dict_obs = {}
    monkeypatch.setattr(cathy_tools, "dictObs_2pd", lambda d: (_ for _ in ()).throw(KeyError(0)))
    monkeypatch.setattr(cathy_tools, "print", lambda *a, **k: None, raising=False)   # rich is stubbed by another test
    cathy_tools.CATHY.backup_results_DA(da, {})          # the class method itself (other tests may stub the DA one)
    files = os.listdir(os.path.join(tmp_path, "prj"))
    assert "prj_dict_parm_pert.pkl" in files and any(f.startswith("prj_df_Archie") for f in files)


# ------------------------------------------------------------------ xarray / NetCDF boundary
def _small_store(nens=3, nnod=5, ncyc=4, **kw):
    st = StateStore(nens, **kw)
    rng = np.random.default_rng(0)
    for c in range(ncyc):
        psi = rng.normal(size=(nnod, nens)); sw = rng.uniform(size=(nnod, nens))
        st.add(c, psi=psi, sw=sw, analysis=psi + 1.0, rejected=[0, 1, 0][:nens])
    return st


def test_to_xarray_matches_arrays_and_masks_rejected():
    st = _small_store()
    ds = st.to_xarray()
    assert ds["sw"].dims == ("cycle", "member", "node")
    assert ds["sw"].shape == (4, 3, 5) and ds["rejected"].shape == (4, 3)
    _, _, _, a = st.arrays("sw")
    np.testing.assert_array_equal(ds["sw"].values, a)
    assert np.isnan(ds["sw"].values[:, 1, :]).all() and not np.isnan(ds["sw"].values[:, 0, :]).any()
    assert list(ds["cycle"].values) == [0, 1, 2, 3]


def test_netcdf_roundtrip_and_load_states_from_chunks(tmp_path):
    from pyCATHY.DA.state_store import load_states
    st = _small_store(chunk_dir=str(tmp_path / "chunks"))
    nc = st.save_netcdf(str(tmp_path / "x_states.nc"))
    a = load_states(nc)
    b = load_states(str(tmp_path / "chunks"))
    for f in ("psi", "sw", "analysis"):
        np.testing.assert_allclose(a[f].values, st.to_xarray()[f].values, rtol=1e-6)
        np.testing.assert_allclose(b[f].values, a[f].values, rtol=1e-6)


def test_df_DA_is_deprecated_but_still_rebuilt(recwarn):
    import warnings
    from pyCATHY.DA.cathy_DA import DA
    da = DA.__new__(DA)
    da._state_store = _small_store()
    da._df_DA_legacy = __import__("pandas").DataFrame()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        frame = da.df_DA
    assert any(issubclass(x.category, DeprecationWarning) for x in w)
    assert len(frame) == 4 * 3 * 5
    assert da.states_dataset()["psi"].shape == (4, 3, 5)
