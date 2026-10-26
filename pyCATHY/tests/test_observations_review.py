"""Review of `read_observations` and `make_data_cov` (round S19).

Every test here first FAILED on the code before the fix (the reason is in the docstring); they pin the behaviour the
data assimilation loop relies on: the covariance of a cycle is built with the SAME sensors, in the SAME order, as the
data vector of `DA._obs_key_select` / `DA._get_data2assimilate`.
"""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA import observations as OBS
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.observations import make_data_cov, read_observations
from pyCATHY.importers import sensors_measures as in_meas


# --------------------------------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------------------------------
def _entry(data, err, kind, **extra):
    e = {"data": data, "data_err": err, "data_type": kind, "units": "-", "mesh_nodes": [1], "sensor_name": kind}
    e.update(extra)
    return e


def _diag_of(stacked, k):
    c = np.asarray(stacked[k])
    return np.diag(c) if c.ndim == 2 else c


def _da_selection(dict_obs, cycle, sel):
    """What the DA loop assimilates at `cycle`: sensor names, in order."""
    da = SimpleNamespace(dict_obs=dict_obs, count_DA_cycle=cycle)
    return DA._obs_key_select(da, sel)[0]


# --------------------------------------------------------------------------------------------------------------------
# read_observations
# --------------------------------------------------------------------------------------------------------------------
def test_discharge_is_read(tmp_path):
    """`discharge` was an `if` apart from the `elif` chain: it fell into the 'unknown data_type' branch and then died on
    an undefined `obs_cov_type` / `filename`."""
    f = tmp_path / "q.tsv"
    pd.DataFrame({"t": [0, 1], "q": [0.1, 0.2]}).to_csv(f, sep="\t", index=False)
    d = read_observations({}, str(f), "discharge", 0.01, tA=0.0)
    e = d[0.0]["discharge"]
    assert e["filename"] == str(f) and e["data_type"] == "discharge" and len(e["data"]) == 2


def test_unknown_data_type_is_an_explicit_error():
    """It logged a warning and then raised UnboundLocalError on `filename` (message about a variable, not the cause)."""
    with pytest.raises(ValueError, match="data_type"):
        read_observations({}, 0.3, "not_a_sensor", 0.01, tA=0.0)


def test_filename_is_recorded_when_a_file_is_given(monkeypatch):
    """`filename` was only defined when `obs_2_add` was NOT a str, so reading a file died on an undefined variable."""
    monkeypatch.setattr(in_meas, "read_swc", lambda f: pd.DataFrame({"swc": [0.3]}), raising=False)
    d = read_observations({}, "some_swc.csv", "swc", 0.01, tA=0.0)
    assert d[0.0]["swc"]["filename"] == "some_swc.csv"


def test_reader_that_does_not_exist_says_so():
    """`in_meas.read_swc` / `read_scale` do not exist: AttributeError from deep inside; now an explicit message."""
    with pytest.raises(NotImplementedError, match="swc"):
        read_observations({}, "some_swc.csv", "swc", 0.01, tA=0.0)


def test_assimilation_time_is_required():
    """tA defaulted to `[]` (unhashable): `TypeError: unhashable type: 'list'` far from the cause."""
    with pytest.raises(ValueError, match="tA"):
        read_observations({}, 0.3, "swc", 0.01)


def test_other_type_at_same_time_does_not_overwrite_a_sensor():
    """The name was bumped only when `data_type in <LAST key>`: with keys {swc, tensio}, a second swc found 'tensio' as the
    last key, kept the name 'swc' and silently REPLACED the first swc."""
    d = {}
    d = read_observations(d, 0.30, "swc", 0.01, tA=0.0)
    d = read_observations(d, 0.31, "tensiometer", 0.01, tA=0.0)
    d = read_observations(d, 0.32, "swc", 0.01, tA=0.0)
    d = read_observations(d, 0.33, "swc", 0.01, tA=0.0)
    d = read_observations(d, 0.34, "tensiometer", 0.01, tA=0.0)
    assert list(d[0.0]) == ["swc", "tensiometer", "swc1", "swc2", "tensiometer1"]
    assert [d[0.0][k]["data"] for k in d[0.0]] == [0.30, 0.31, 0.32, 0.33, 0.34]
    assert all(d[0.0][k]["sensor_name"] == k for k in d[0.0])


def test_sensor_numbering_survives_a_gap():
    """`swc` then `swc2` present (no `swc1`): the next one must not collide with either."""
    d = {0.0: {"swc": _entry(1, 1, "swc"), "swc2": _entry(2, 1, "swc")}}
    d = read_observations(d, 0.5, "swc", 0.01, tA=0.0)
    assert list(d[0.0]) == ["swc", "swc2", "swc3"]


def test_ERT_without_meta_does_not_raise_keyerror(monkeypatch):
    """`'data_format' in kwargs['meta']` raised KeyError when `meta` was not given (optional kwarg)."""
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4], "m": [5, 6], "n": [7, 8], "resist": [1.0, 2.0],
                       "rhoa": [10.0, 20.0], "rec_err": [0.1, 0.1]})
    monkeypatch.setattr(in_meas, "read_ERT", lambda f, fmt: (df, {}))
    d = read_observations({}, "x.csv", "ERT", 5, tA=0.0)
    assert d[0.0]["ERT"]["data_format"] == "resipy"


# --------------------------------------------------------------------------------------------------------------------
# make_data_cov
# --------------------------------------------------------------------------------------------------------------------
def _two_sensor_dict():
    return {0.0: {"swc": _entry(0.3, 0.01, "swc"), "tensio": _entry(-1.0, 0.5, "tensio")},
            3600.0: {"swc": _entry(0.31, 0.01, "swc"), "tensio": _entry(-1.1, 0.5, "tensio")}}


def test_a_string_selects_that_sensor_only():
    """`list_assimilated_obs='swc'` (the form used in the documentation) was iterated CHARACTER by character: 's' in
    'tensio' is True, so the tensiometer was assimilated too."""
    _, _, st = make_data_cov(SimpleNamespace(), _two_sensor_dict(), list_assimilated_obs="swc")
    assert [len(_diag_of(st, k)) for k in range(2)] == [1, 1]
    np.testing.assert_allclose(_diag_of(st, 0), [0.01 ** 2])


@pytest.mark.parametrize("sel", ["all", ["all"], ["swc"], ["tensio"], ["swc", "tensio"], "tensio"])
def test_selection_is_the_one_of_the_da_loop(sel):
    """`make_data_cov` used substring matching and `== 'all'`; the DA uses the root name of the sensor and `'all' in`.
    ['all'] gave NO sensor (crash in np.hstack)."""
    d = _two_sensor_dict()
    _, _, st = make_data_cov(SimpleNamespace(), d, list_assimilated_obs=sel)
    for k in range(2):
        assert len(_diag_of(st, k)) == len(_da_selection(d, k, sel))


def test_root_name_matching_not_substring():
    """'tensio' must not select a sensor whose root is 'tensiometer' (the DA does not), otherwise R and the data differ."""
    d = {0.0: {"tensiometer": _entry(-1.0, 0.5, "tensiometer"), "swc": _entry(0.3, 0.01, "swc")}}
    assert _da_selection(d, 0, ["tensio"]) == []
    with pytest.raises(ValueError, match="no sensor"):
        make_data_cov(SimpleNamespace(), d, list_assimilated_obs=["tensio"])


def test_sensor_order_follows_the_dict_of_each_time():
    """The data vector of a cycle is stacked in the insertion order of the sensors AT THAT TIME. R was built in the order of
    first appearance over ALL times, so a sensor added in another order at a later time got the wrong variance (same total
    size: nothing raised)."""
    d = {0.0: {"swc": _entry(0.3, 0.01, "swc"), "tensio": _entry(-1.0, 0.5, "tensio")},
         3600.0: {"tensio": _entry(-1.1, 0.5, "tensio"), "swc": _entry(0.31, 0.01, "swc")}}
    _, _, st = make_data_cov(SimpleNamespace(), d)
    np.testing.assert_allclose(_diag_of(st, 0), [0.01 ** 2, 0.5 ** 2])
    assert _da_selection(d, 1, "all") == ["tensio", "swc"]
    np.testing.assert_allclose(_diag_of(st, 1), [0.5 ** 2, 0.01 ** 2])


def test_times_are_chronological_like_the_da():
    """The DA refuses an unsorted dict_obs, so the covariances are stacked in chronological order whatever the input order."""
    d = _two_sensor_dict()
    rev = {t: d[t] for t in sorted(d, reverse=True)}
    d[0.0]["swc"]["data_err"] = 0.02
    _, _, a = make_data_cov(SimpleNamespace(), d)
    _, _, b = make_data_cov(SimpleNamespace(), rev)
    for k in range(2):
        np.testing.assert_array_equal(_diag_of(a, k), _diag_of(b, k))
    np.testing.assert_allclose(_diag_of(a, 0)[0], 0.02 ** 2)


def test_vector_valued_sensor_gets_one_variance_per_datum():
    """A sensor holding n values got ONE variance (`len([array])` == 1): R was (1, 1) for n data."""
    d = {0.0: {"swc": _entry(np.array([0.3, 0.31, 0.32]), 0.01, "swc")}}
    _, _, st = make_data_cov(SimpleNamespace(), d)
    np.testing.assert_allclose(_diag_of(st, 0), [0.01 ** 2] * 3)


def test_resipy_ERT_uses_resist_for_a_relative_error():
    """The DA assimilates `resist` for resipy data (`rhoa` for pygimli); R read `rhoa` always: KeyError when the column is
    absent, a variance proportional to the wrong quantity when it is present."""
    df = pd.DataFrame({"resist": [2.0, 4.0], "rhoa": [200.0, 400.0]})
    d = {0.0: {"ERT": _entry(df, 10.0, "ERT", data_format="resipy", data_err_unit="percent")}}
    _, _, st = make_data_cov(SimpleNamespace(), d)
    np.testing.assert_allclose(_diag_of(st, 0), [(0.1 * 2.0) ** 2, (0.1 * 4.0) ** 2])
    d = {0.0: {"ERT": _entry(df[["resist"]], 10.0, "ERT", data_format="resipy")}}
    _, _, st = make_data_cov(SimpleNamespace(), d)
    assert len(_diag_of(st, 0)) == 2


def test_pygimli_ERT_uses_rhoa():
    df = pd.DataFrame({"rhoa": [20.0, 40.0], "resist": [1.0, 2.0]})
    d = {0.0: {"ERT": _entry(df, 10.0, "ERT", data_format="pygimli", data_err_unit="percent")}}
    _, _, st = make_data_cov(SimpleNamespace(), d)
    np.testing.assert_allclose(_diag_of(st, 0), [2.0 ** 2, 4.0 ** 2])


def test_EM_coils_come_from_the_entry():
    """The coils are stored in the entry by `read_observations`; asking for them again as a kwarg was mandatory."""
    df = pd.DataFrame({"HCP0.5": [1.0, 2.0], "HCP1.0": [3.0, 4.0]})
    d = {0.0: {"EM": _entry(df, 1.0, "EM", coils=["HCP0.5", "HCP1.0"])}}
    _, _, st = make_data_cov(SimpleNamespace(), d)
    assert len(_diag_of(st, 0)) == 4
    _, _, st2 = make_data_cov(SimpleNamespace(), d, coils=["HCP0.5"])      # an explicit kwarg still wins
    assert len(_diag_of(st2, 0)) == 2


def test_simu_receives_the_stacked_covariance():
    s = SimpleNamespace()
    _, _, st = make_data_cov(s, _two_sensor_dict())
    assert s.stacked_data_cov is st


# --------------------------------------------------------------------------------------------------------------------
# DA._obs_key_select now delegates to select_sensors: same result as the code it replaced (kept here as the reference)
# --------------------------------------------------------------------------------------------------------------------
def _old_obs_key_select(sensors, list_assimilated_obs):
    import re
    out = []
    for sensor in sensors:
        if "all" in list_assimilated_obs:
            out.append(sensor)
        elif re.sub(r"\d", "", sensor) in list_assimilated_obs:
            out.append(sensor)
    return out


@pytest.mark.parametrize("sel", [["all"], ["swc"], ["tensiometer"], ["ERT"], ["swc", "ERT"], ["EM"], [], ("swc",),
                                 ["swc", "all"], np.array(["swc", "ERT"])])
def test_obs_key_select_unchanged_for_lists(sel):
    d = {0.0: {"swc": _entry(1, 1, "swc"), "tensiometer": _entry(1, 1, "tensiometer"), "swc1": _entry(1, 1, "swc"),
               "ERT": _entry(1, 1, "ERT"), "swc12": _entry(1, 1, "swc")}}
    assert _da_selection(d, 0, sel) == _old_obs_key_select(list(d[0.0]), list(sel))


# --------------------------------------------------------------------------------------------------------------------
# round S19d: neighbours of the same bug classes
# --------------------------------------------------------------------------------------------------------------------
def test_unknown_obs_cov_type_is_an_explicit_error():
    """Any other string than the two known ones raised UnboundLocalError on `data_cov`."""
    from pyCATHY.DA.observations import init_obs_cov_matrice
    d = {"data": pd.DataFrame({"rec_err": [1.0, 2.0]}), "data_err": 3.0}
    with pytest.raises(ValueError, match="obs_cov_type"):
        init_obs_cov_matrice(d, "recip_err")
    assert init_obs_cov_matrice(d, "reciprocal_err").shape == (2, 2)
    assert init_obs_cov_matrice(d, "data_err").shape == (2, 2)


def _sync_inputs():
    dm = {100.0: {"swc": {"assimilation_times": 100.0, "datetime": 1}, "ERT": {"assimilation_times": 100.0}},
          250.0: {"swc": {"assimilation_times": 250.0}}}
    return dm, pd.DataFrame({"datetime": [pd.Timestamp("2023-01-01")]})


def test_resynchronise_shifts_keys_and_entries_without_touching_the_input():
    """`dict(data_measure)` was a shallow copy: the caller's entries were rewritten (key 100.0 but assimilation_times 50.0)."""
    from pyCATHY.cathy_utils import resynchronise_times
    dm, atm = _sync_inputs()
    out = resynchronise_times(dm, atm, 50.0)
    assert list(out) == [50.0, 200.0]
    assert out[50.0]["swc"]["assimilation_times"] == 50.0 and out[200.0]["swc"]["assimilation_times"] == 200.0
    assert dm[100.0]["swc"]["assimilation_times"] == 100.0 and list(dm) == [100.0, 250.0]      # input intact


@pytest.mark.parametrize("offset", [None, "x"])
def test_resynchronise_failure_names_the_real_cause(offset):
    """The bare `except:` handler indexed `data_measure[0]` (KeyError: 0), hiding the TypeError, or silently went on with
    the unsynchronised times."""
    from pyCATHY.cathy_utils import resynchronise_times
    dm, atm = _sync_inputs()
    with pytest.raises(ValueError, match="cannot synchronise") as exc:
        resynchronise_times(dm, atm, offset)
    assert isinstance(exc.value.__cause__, TypeError)
    assert dm[100.0]["swc"]["assimilation_times"] == 100.0


def test_legacy_performance_module_imports_os():
    import pyCATHY.DA.performance as perf
    assert hasattr(perf, "os")
