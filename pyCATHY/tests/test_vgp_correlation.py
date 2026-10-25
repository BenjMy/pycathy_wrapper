"""Correlated van Genuchten parameters (Botto et al. 2018): vgp_correlation.py."""
import numpy as np
import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA import vgp_correlation as vc

R = np.array([[1.0, 0.6, -0.3], [0.6, 1.0, 0.2], [-0.3, 0.2, 1.0]])


def entry(name, mean, sd, lo, hi, transf, n):
    return dict(type_parm=name, mean=mean, nominal=mean, sd=sd, clip_min=lo, clip_max=hi, transf_type=[transf] * 3,
                sampling_type="normal", ensemble_size=n, ini_perturbation=np.random.normal(mean, sd, n),
                sampling=np.zeros(n))


def make(n, groups=(0, 1)):
    d = {}
    for g in groups:
        d[f"SATCELL_VG{g}"] = entry(f"SATCELL_VG{g}", 0.5, 0.1, 0.05, 2.0, "log", n)
        d[f"n_VG{g}"] = entry(f"n_VG{g}", 1.6, 0.1, 1.05, 2.5, "log-ratio", n)
        d[f"thetar_VG{g}"] = entry(f"thetar_VG{g}", 0.05, 0.02, 0.0, 0.15, "None", n)
    return d


def test_sample_correlation_matches_target_in_transformed_space():
    np.random.seed(0)
    d = make(20000)
    d, rep = vc.apply_vgp_correlation(d, R)
    assert np.allclose(rep[0]["sample"], R, atol=0.03) and np.allclose(rep[1]["sample"], R, atol=0.03)
    for k in ("SATCELL_VG0", "n_VG0", "thetar_VG0"):                 # bounds respected, ensemble stored in both keys
        a = d[k]["ini_perturbation"]
        assert a.min() >= d[k]["clip_min"] and a.max() <= d[k]["clip_max"]
        assert np.array_equal(a, d[k]["sampling"])
    assert d["n_VG0"]["ini_perturbation"].mean() == pytest.approx(1.6, abs=0.03)       # centred on the nominal
    assert d["thetar_VG0"]["ini_perturbation"].std() == pytest.approx(0.02, rel=0.1)    # marginal sd of the scenario


def test_identity_gives_uncorrelated_and_per_group_matrices():
    np.random.seed(1)
    d = make(20000)
    d, rep = vc.apply_vgp_correlation(d, {0: np.eye(3), "default": R})
    off = rep[0]["sample"] - np.eye(3)
    assert np.abs(off).max() < 0.03
    assert np.allclose(rep[1]["sample"], R, atol=0.03)


def test_sigma_vg_scales_the_spread_and_seed_reproduces():
    np.random.seed(2); a = vc.apply_vgp_correlation(make(500), R, sigma_vg=0.5)[0]
    np.random.seed(2); b = vc.apply_vgp_correlation(make(500), R, sigma_vg=0.5)[0]
    np.random.seed(2); c = vc.apply_vgp_correlation(make(500), R, sigma_vg=1.0)[0]
    assert np.array_equal(a["n_VG0"]["ini_perturbation"], b["n_VG0"]["ini_perturbation"])
    assert a["thetar_VG0"]["ini_perturbation"].std() < 0.75 * c["thetar_VG0"]["ini_perturbation"].std()


def test_invalid_matrices_are_rejected_and_semidefinite_repaired():
    with pytest.raises(ValueError):
        vc.check_correlation(np.array([[1, 0.5], [0.4, 1]]))                 # not symmetric
    with pytest.raises(ValueError):
        vc.check_correlation(np.array([[2.0, 0], [0, 1]]))                   # diagonal != 1
    with pytest.raises(ValueError):
        vc.apply_vgp_correlation(make(10), np.eye(2))                          # wrong size
    bad = np.array([[1, 0.9, -0.9], [0.9, 1, 0.9], [-0.9, 0.9, 1]])             # indefinite
    with pytest.warns(RuntimeWarning):
        fixed = vc.check_correlation(bad)
    assert np.linalg.eigvalsh(fixed).min() > 0 and np.allclose(np.diag(fixed), 1)


def test_incomplete_group_is_left_independent():
    np.random.seed(3)
    d = make(50, groups=(0,))
    del d["n_VG0"]
    before = d["SATCELL_VG0"]["ini_perturbation"].copy()
    d, rep = vc.apply_vgp_correlation(d, R)
    assert rep == {} and np.array_equal(d["SATCELL_VG0"]["ini_perturbation"], before)


def test_johnson_helpers_roundtrip_and_accept_both_names():
    x = np.array([0.2, 0.5, 1.5])
    for kind in ("LN", "SB", "SU", "Johnson_LN", "Johnson_SB", "Johnson_SU"):
        y = perturbate.Johnson1970_transform(x, kind, A=0.0, B=2.0)
        back = perturbate.inverse_Johnson1970_transform(y, kind, A=0.0, B=2.0)
        assert np.allclose(back, x)


def make_botto(n, groups=(0, 1)):
    d = {}
    for g in groups:
        d[f"SATCELL_VG{g}"] = entry(f"SATCELL_VG{g}", 0.44, 0.1, 0.10, 2.0, "log-ratio", n)
        d[f"n_VG{g}"] = entry(f"n_VG{g}", 1.4, 0.1, 1.05, 2.5, "log-ratio", n)
        d[f"thetar_VG{g}"] = entry(f"thetar_VG{g}", 0.04, 0.02, 0.0, 0.15, "log-ratio", n)
    return d


@pytest.mark.parametrize("cls", ["sandy_loam", "silt_loam"])
def test_botto_table2_correlation_reproduced(cls):
    np.random.seed(3)
    d, rep = vc.apply_botto_vgp(make_botto(40000), cls, sigma_vg=0.25)
    cov, corr = vc.botto_covariance(cls)
    assert np.allclose(np.diag(corr), 1) and np.allclose(cov, cov.T)
    assert np.allclose(rep[0]["sample"], corr, atol=0.05)
    for k in ("SATCELL_VG0", "n_VG0", "thetar_VG0"):
        a = d[k]["ini_perturbation"]
        assert a.min() >= d[k]["clip_min"] and a.max() <= d[k]["clip_max"]
        assert np.array_equal(a, d[k]["sampling"])


def test_botto_sigma_scales_and_per_group_class():
    np.random.seed(4)
    d1 = vc.apply_botto_vgp(make_botto(5000), {0: "sandy_loam", 1: "silt_loam"}, sigma_vg=0.1)[0]
    d2 = vc.apply_botto_vgp(make_botto(5000), {0: "sandy_loam", 1: "silt_loam"}, sigma_vg=0.25)[0]
    assert d2["n_VG0"]["ini_perturbation"].std() > 1.5 * d1["n_VG0"]["ini_perturbation"].std()
    assert d1["SATCELL_VG1"]["correlated_with"]
    assert not hasattr(vc, "PLACEHOLDER_CORR")


def test_botto_links_to_petro_carsel_parrish():
    cp = vc.carsel_parrish_bounds("silt_loam")
    assert cp["thetar_VG"] == (0.0, 0.11, "log-ratio") and cp["n_VG"] == (1.35, 2.0, "log-ratio") and cp["alpha"][2] == "log"
    assert vc.carsel_parrish_bounds("sandy_loam") is None                 # not tabulated in petro: scenario bounds
    from pyCATHY import petro
    assert petro.Carsel_Parrish_1988("SIL")[2]["VGPSATCELL"] == pytest.approx(0.5)      # 1/(0.02 cm-1) = 0.5 m
    np.random.seed(5)
    d, rep = vc.apply_botto_vgp(make_botto(20000), "silt_loam", sigma_vg=0.25)
    assert d["n_VG0"]["ini_perturbation"].min() >= 1.35 and d["thetar_VG0"]["ini_perturbation"].max() <= 0.11
    assert np.allclose(rep[0]["sample"], rep[0]["target"], atol=0.06)


def _stack(n, spec, sigma=None, seed=7):
    """Build the VGP entries through perturbate_parm (as nb7 does) with the correlation set in perturbate."""
    from pyCATHY.DA import perturbate as P
    np.random.seed(seed)
    P.set_vgp_correlation(spec, sigma)
    try:
        out = {}
        for nm, mean, sd, lo, hi in (("SATCELL_VG0", 0.44, 0.1, 0.10, 2.0), ("n_VG0", 1.4, 0.1, 1.05, 2.5),
                                     ("thetar_VG0", 0.04, 0.02, 0.0, 0.15)):
            parm = dict(type_parm=nm, nominal=mean, mean=mean, sd=sd, clip_min=lo, clip_max=hi, sampling_type="normal",
                        per_type=None, ensemble_size=n, transf_type=["log-ratio"] * 3, savefig=False)
            out = P.perturbate_parm(out, parm=parm, type_parm=nm, mean=mean, sd=sd, sampling_type="normal",
                                    ensemble_size=n, per_type=None)
    finally:
        P.set_vgp_correlation(None)
    return out


def test_hook_in_perturbate_parm_default_vs_botto(capsys):
    d0 = _stack(3000, None)
    assert "DEFAULT" in capsys.readouterr().out
    keys = [k for k in d0 if k.startswith(("SATCELL_VG", "n_VG", "thetar_VG"))]
    assert keys and all(d0[k].get("vgp_corr_mode") == "independent" for k in keys)
    d1 = _stack(3000, "botto:silt_loam")
    out = capsys.readouterr().out
    assert "NOT default" in out and "Botto 2018 silt_loam" in out
    assert all(d1[k]["vgp_corr_mode"].startswith("botto:silt_loam") for k in keys)
    a, b = d1["SATCELL_VG0"]["ini_perturbation"], d1["thetar_VG0"]["ini_perturbation"]
    assert abs(np.corrcoef(1 / a, b)[0, 1]) > 0.3        # alpha-theta_r correlation of Table 2 (0.017/0.538 scale) is present
