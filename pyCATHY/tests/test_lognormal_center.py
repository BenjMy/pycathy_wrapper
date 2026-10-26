"""Lognormal prior: nominal value = median (default, historical) or = ensemble mean (``center="mean"``)."""
import numpy as np
import pytest

from pyCATHY.DA import perturbate as P


def test_default_is_bit_identical_to_the_historical_draw():
    np.random.seed(0); a = P.sampling_dist("lognormal", 1.3, 0.7, 50)
    np.random.seed(0); b = np.random.lognormal(np.log(1.3), 0.7, 50)
    np.testing.assert_array_equal(a, b)
    assert P.lognormal_log_location(1.3, 0.7) == np.log(1.3)


@pytest.mark.parametrize("sd", [0.3, 0.7, 1.0])
def test_median_centred_mean_is_biased_up_and_mean_centred_is_not(sd):
    np.random.seed(1)
    med = P.sampling_dist("lognormal", 2.0, sd, 400_000)
    mean = P.sampling_dist("lognormal", 2.0, sd, 400_000, center="mean")
    assert med.mean() / 2.0 == pytest.approx(np.exp(sd ** 2 / 2), rel=0.02)       # +5 %, +28 %, +65 %
    assert np.median(med) / 2.0 == pytest.approx(1.0, rel=0.01)
    assert mean.mean() / 2.0 == pytest.approx(1.0, rel=0.02)
    assert np.median(mean) / 2.0 == pytest.approx(np.exp(-sd ** 2 / 2), rel=0.01)


def test_location_formula():
    assert P.lognormal_log_location(2.0, 0.7, "mean") == pytest.approx(np.log(2.0) - 0.245)
    assert P.lognormal_log_location(2.0, 0.7, "median") == np.log(2.0)


def test_unknown_center_is_an_error():
    with pytest.raises(ValueError, match="lognormal_center"):
        P.lognormal_log_location(1.0, 0.3, "average")
    with pytest.raises(ValueError):
        P.set_lognormal_center("mode")


def test_context_manager_switches_and_restores():
    with P.using_lognormal_center("mean"):
        assert P.lognormal_log_location(1.0, 1.0) == pytest.approx(-0.5)
    assert P.lognormal_log_location(1.0, 1.0) == 0.0
    with pytest.raises(RuntimeError):
        with P.using_lognormal_center("mean"):
            raise RuntimeError
    assert P.lognormal_log_location(1.0, 1.0) == 0.0                              # restored after an error


def test_entry_key_beats_the_process_default():
    with P.using_lognormal_center("mean"):
        np.random.seed(2); a = P.sampling_dist("lognormal", 1.0, 0.7, 10, center="median")
    np.random.seed(2); b = P.sampling_dist("lognormal", 1.0, 0.7, 10)
    np.testing.assert_array_equal(a, b)


def test_layered_sampler_follows_the_option():
    kw = dict(per_type=None, nominal=2.0, sampling_type="lognormal")
    np.random.seed(0); med = P.perturbate_parm_by_layers("Ks", 200_000, 3, mean=2.0, sd=0.7, **kw)
    np.random.seed(0); mean = P.perturbate_parm_by_layers("Ks", 200_000, 3, mean=2.0, sd=0.7, lognormal_center="mean", **kw)
    assert med.mean() / 2.0 == pytest.approx(np.exp(0.245), rel=0.02) and mean.mean() / 2.0 == pytest.approx(1.0, rel=0.02)
    per_layer = P.perturbate_parm_by_layers("Ks", 100_000, 3, mean=np.array([[2.], [4.], [8.]]),
                                            sd=np.array([[.7], [.3], [1.]]), lognormal_center="mean", **kw)
    np.testing.assert_allclose(per_layer.mean(axis=1) / np.array([2, 4, 8]), 1.0, rtol=0.03)


def test_truncated_vg_path_reads_the_entry_key():
    base = dict(sampling_type="lognormal", mean=2.0, sd=0.5, per_type=None, nominal=2.0, clip_min=None, clip_max=None)
    out = {}
    for key in (None, "median", "mean"):
        parm = dict(base, **({} if key is None else {"lognormal_center": key}))
        np.random.seed(3)
        out[key], _ = P.VG_pert_rules({}, parm, "n_VG", 200_000, 2.0, 0.5, None, "lognormal")
    np.testing.assert_array_equal(out[None], out["median"])                       # default unchanged
    assert out["mean"].mean() / 2.0 == pytest.approx(1.0, rel=0.02)
    assert out["median"].mean() / 2.0 == pytest.approx(np.exp(0.125), rel=0.02)
