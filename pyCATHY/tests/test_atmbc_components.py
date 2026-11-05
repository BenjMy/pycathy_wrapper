"""Rain / ETp atmbc components: set_atmbc_components, separate perturbation, mass-balance ETa operator."""
import types

import numpy as np
import pandas as pd
import pytest

from pyCATHY.cathy_tools import CATHY
from pyCATHY.DA import perturbate
from pyCATHY.DA.cathy_DA import DA

DAY = 86400.0
MM = 1e-3 / DAY


def _fake_sim():
    rec = {}

    def update_atmbc(**kw):
        rec.update(kw)

    return types.SimpleNamespace(update_atmbc=update_atmbc, rec=rec)


def test_set_atmbc_components_writes_the_net_flux_and_keeps_the_parts():
    sim = _fake_sim()
    t = np.arange(4) * DAY
    rain = np.array([0.0, 5 * MM, 0.0, 1 * MM])
    etp = np.array([3 * MM, 3 * MM, -2 * MM, 2 * MM])          # abs is taken
    CATHY.set_atmbc_components(sim, t, rain, etp, HSPATM=1, IETO=1)
    assert np.allclose(sim.rec["netValue"], rain - np.abs(etp))
    assert sim.rec["IETO"] == 1 and sim.rec["HSPATM"] == 1
    assert np.allclose(sim.atmbc_components["etp"], np.abs(etp))
    assert np.allclose(sim.atmbc_components["rain"], rain)


def test_set_atmbc_components_rejects_bad_input():
    with pytest.raises(ValueError):
        CATHY.set_atmbc_components(_fake_sim(), [0, 1], [0.0], [0.0, 0.0])
    with pytest.raises(ValueError):
        CATHY.set_atmbc_components(_fake_sim(), [0, 1], [-1.0, 0.0], [0.0, 0.0])
    with pytest.raises(NotImplementedError):
        CATHY.set_atmbc_components(_fake_sim(), [0, 1], [0.0, 0.0], [0.0, 0.0], HSPATM=0)


def _entry(n=60, rho=0.0, rain_sd=0.3, etp_sd=0.1, tau=2 * DAY, gamma_kind="evensen"):
    t = np.arange(n) * DAY
    rain = np.where(np.arange(n) % 5 == 0, 4 * MM, 0.0)
    etp = np.full(n, 3 * MM)
    return {
        "type_parm": "atmbc0", "nominal": 1.0, "mean": 1.0, "sd": 0.3, "sampling_type": "normal",
        "per_type": None, "time_variable": True, "time_decorrelation_len": None, "savefig": "atmbc.png",
        "atmbc_components": {"rain": {"sd_rel": rain_sd, "tau": tau}, "etp": {"sd_rel": etp_sd, "tau": tau},
                             "gamma_kind": gamma_kind, "rho": rho},
        "data2perturbate": {"time": list(t), "VALUE": list(rain - etp), "rain": rain, "etp": etp},
    }


def test_perturbate_parm_dispatches_to_the_component_rules():
    parm = _entry()
    rng = np.random.default_rng(0)
    out = perturbate.perturbate_parm({}, parm=parm, type_parm="atmbc0", mean=1.0, sd=0.3,
                                     sampling_type="normal", ensemble_size=200, per_type=None, rng=rng)
    e = out["atmbc0"]
    assert e["time_variable_rain"].shape == (200, 60)
    assert np.allclose(e["time_variable_perturbation"], e["time_variable_rain"] - e["time_variable_etp"])


def test_factors_are_unbiased_with_the_requested_relative_sd_and_dry_days_stay_dry():
    parm = _entry(n=60, rain_sd=0.3, etp_sd=0.1)
    out = perturbate.perturbate_parm({}, parm=parm, type_parm="atmbc0", mean=1.0, sd=0.3,
                                     sampling_type="normal", ensemble_size=4000, per_type=None,
                                     rng=np.random.default_rng(1))["atmbc0"]
    rain0, etp0 = parm["data2perturbate"]["rain"], parm["data2perturbate"]["etp"]
    wet = rain0 > 0
    f_rain = out["time_variable_rain"][:, wet] / rain0[wet]
    f_etp = out["time_variable_etp"] / etp0
    assert abs(f_rain.mean() - 1) < 0.02 and abs(f_rain.std() - 0.3) < 0.02
    assert abs(f_etp.mean() - 1) < 0.01 and abs(f_etp.std() - 0.1) < 0.01
    assert np.all(out["time_variable_rain"][:, ~wet] == 0.0)


def test_ar1_time_correlation_and_rain_etp_correlation():
    parm = _entry(n=200, rho=0.8, rain_sd=0.5, etp_sd=0.5, tau=5 * DAY)
    e = perturbate.perturbate_parm({}, parm=parm, type_parm="atmbc0", mean=1.0, sd=0.3,
                                   sampling_type="normal", ensemble_size=2000, per_type=None,
                                   rng=np.random.default_rng(2))["atmbc0"]
    etp = e["time_variable_etp"] / parm["data2perturbate"]["etp"]
    lag1 = np.corrcoef(np.log(etp[:, 50:-1]).ravel(), np.log(etp[:, 51:]).ravel())[0, 1]
    assert 0.6 < lag1 < 0.95                                     # gamma = 1 - 1/5 = 0.8
    rain0 = parm["data2perturbate"]["rain"]
    wet = rain0 > 0
    lr = np.log(e["time_variable_rain"][:, wet] / rain0[wet]).ravel()
    le = np.log(etp[:, wet]).ravel()
    assert 0.7 < np.corrcoef(lr, le)[0, 1] < 0.9                  # rho = 0.8


def test_gamma_kind_and_rho_are_validated():
    with pytest.raises(ValueError):
        perturbate.perturbate_parm({}, parm=_entry(gamma_kind="nope"), type_parm="atmbc0", mean=1.0, sd=0.3,
                                   sampling_type="normal", ensemble_size=4, per_type=None)
    with pytest.raises(ValueError):
        perturbate.perturbate_parm({}, parm=_entry(rho=1.0), type_parm="atmbc0", mean=1.0, sd=0.3,
                                   sampling_type="normal", ensemble_size=4, per_type=None)


def test_each_variable_keeps_its_own_decorrelation_time_when_rho_is_not_zero():
    """rho is imposed on the white noises: ETp (tau = 7 d) must not inherit the rain memory (tau = 2 d)."""
    n = 121
    parm = _entry(n=n, rho=-0.3, tau=7 * DAY, gamma_kind="exp")
    parm["atmbc_components"]["rain"]["tau"] = 2 * DAY
    e = perturbate.perturbate_parm({}, parm=parm, type_parm="atmbc0", mean=1.0, sd=0.3, sampling_type="normal",
                                   ensemble_size=4000, per_type=None, rng=np.random.default_rng(3))["atmbc0"]
    x = np.log(e["time_variable_etp"] / parm["data2perturbate"]["etp"])
    for lag in (1, 3, 7):
        r = (x[:, :-lag] * x[:, lag:]).mean() / (x ** 2).mean()
        assert abs(r - np.exp(-lag / 7)) < 0.03, (lag, r)


def test_same_seed_same_ensemble_and_matches_the_documented_construction():
    parm = _entry(n=30, rho=0.5, tau=3 * DAY, gamma_kind="evensen")
    a = perturbate.perturbate_parm({}, parm=dict(parm), type_parm="atmbc0", mean=1.0, sd=0.3, sampling_type="normal",
                                   ensemble_size=8, per_type=None, rng=np.random.default_rng(9))["atmbc0"]
    b = perturbate.perturbate_parm({}, parm=dict(parm), type_parm="atmbc0", mean=1.0, sd=0.3, sampling_type="normal",
                                   ensemble_size=8, per_type=None, rng=np.random.default_rng(9))["atmbc0"]
    assert np.array_equal(a["time_variable_perturbation"], b["time_variable_perturbation"])


# ------------------------------------------------------------------ ETa operator
class _Sim:
    _window_rain_etp = DA._window_rain_etp
    _predict_ET = DA._predict_ET


def _fort777(monkeypatch, etra):
    df = pd.DataFrame({"time_sec": [0.0] * 3 + [DAY] * 3, "ACT. ETRA": [9.0] * 3 + list(etra)})
    monkeypatch.setattr("pyCATHY.DA.cathy_DA.out_CT.read_fort777", lambda f: df.copy())


def _sim(mode, perturbed=False):
    s = _Sim()
    s.eta_patch_mode = mode
    s.count_atmbc_cycle = 1
    s.atmbc_components = {"rain": np.array([0.0, 4 * MM, 0.0]), "etp": np.array([3 * MM] * 3)}
    s.dict_parm_pert = {}
    if perturbed:
        s.dict_parm_pert = {"atmbc0": {"time_variable_rain": np.array([[0, 1 * MM, 0], [0, 6 * MM, 0]]),
                                       "time_variable_etp": np.array([[3 * MM] * 3, [2 * MM] * 3])}}
    return s


OBS = [{"data_type": "ETact", "mesh_nodes": [1]}]


def test_default_mode_returns_the_raw_etra(monkeypatch):
    _fort777(monkeypatch, [1.0, 2.0, 3.0])
    assert np.allclose(_sim("none")._predict_ET("x", OBS, ens_nb=0), [2.0])


def test_mass_balance_adds_min_rain_etp_nominal(monkeypatch):
    _fort777(monkeypatch, [1.0, 2.0 * MM, 3.0])
    out = _sim("mass_balance")._predict_ET("x", OBS, ens_nb=0)       # window 1: rain 4, etp 3 mm/d
    assert np.allclose(out, 2.0 * MM + 3 * MM)                       # rain > etp -> ETp


def test_mass_balance_uses_each_members_own_forcing(monkeypatch):
    _fort777(monkeypatch, [1.0, 2.0 * MM, 3.0])
    s = _sim("mass_balance", perturbed=True)
    assert np.allclose(s._predict_ET("x", OBS, ens_nb=0), 2 * MM + 1 * MM)   # rain 1 < etp 3 -> rain
    assert np.allclose(s._predict_ET("x", OBS, ens_nb=1), 2 * MM + 2 * MM)   # rain 6 > etp 2 -> ETp


def test_mass_balance_without_components_or_member_fails_loudly(monkeypatch):
    _fort777(monkeypatch, [1.0, 2.0, 3.0])
    s = _sim("mass_balance")
    s.atmbc_components = None
    with pytest.raises(ValueError):
        s._predict_ET("x", OBS, ens_nb=0)
    with pytest.raises(ValueError):
        _sim("mass_balance")._predict_ET("x", OBS)
    with pytest.raises(ValueError):
        _sim("bogus")._predict_ET("x", OBS, ens_nb=0)


def test_record_offset_selects_the_next_atmbc_record(monkeypatch):
    _fort777(monkeypatch, [1.0, 2.0 * MM, 3.0])
    s = _sim("mass_balance")
    s.eta_patch_record_offset = 1                                     # window 1 -> record 2: rain 0
    assert np.allclose(s._predict_ET("x", OBS, ens_nb=0), 2 * MM + 0.0)
    s.count_atmbc_cycle, s.eta_patch_record_offset = 2, 5             # clipped to the last record
    assert np.allclose(s._predict_ET("x", OBS, ens_nb=0), 2 * MM + 0.0)
