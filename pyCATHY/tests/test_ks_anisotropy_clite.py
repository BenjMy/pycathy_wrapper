"""Ks "c-lite": PERMX = PERMY = Ks, PERMZ = Ks * base kz/kx of each row (anisotropy kept)."""
import numpy as np
import pandas as pd

from pyCATHY.DA.param_tables import _assign_soil_value

IDX = pd.MultiIndex.from_product([[1, 2], range(3)], names=["zone", "layer"])
COLS = ["PERMX", "PERMY", "PERMZ"]


def table(with_ratio=True):
    kx = np.array([3e-7] * 3 + [0.5e-8] * 3)
    kz = np.array([3e-7] * 3 + [0.5e-8 * 0.017] * 3)
    df = pd.DataFrame({"PERMX": kx, "PERMY": kx, "PERMZ": kz}, index=IDX)
    ratio = kz / kx
    df[COLS] = np.nan
    if with_ratio:
        df.attrs["kz_ratio"] = ratio
    return df, ratio


def test_zone_material_keeps_ratio():
    df, ratio = table()
    _assign_soil_value(df, "zone_material", 0, COLS, 1e-6, layer_groups=[[2, [0, 1]]])
    sel = [(2, 0), (2, 1)]
    assert np.allclose(df.loc[sel, "PERMX"], 1e-6) and np.allclose(df.loc[sel, "PERMY"], 1e-6)
    assert np.allclose(df.loc[sel, "PERMZ"], 1e-6 * 0.017)
    assert df.loc[(1, 0), COLS].isna().all()          # other rows untouched


def test_all_rows_and_layers():
    df, ratio = table()
    _assign_soil_value(df, "layers", 1, COLS, 2e-6)
    for z in (1, 2):
        assert np.isclose(df.loc[(z, 1), "PERMZ"], 2e-6 * ratio[(z - 1) * 3 + 1])
    df2, _ = table()
    _assign_soil_value(df2, "all", None, COLS, 5e-7)
    assert np.allclose(df2["PERMZ"], 5e-7 * ratio)


def test_without_ratio_isotropic_as_before():
    df, _ = table(with_ratio=False)
    _assign_soil_value(df, "zone", 1, COLS, 1e-6)
    assert np.allclose(df.loc[(2, slice(None)), "PERMZ"], 1e-6)


def test_other_parameters_unaffected():
    df, _ = table()
    df["POROS"] = np.nan
    _assign_soil_value(df, "layers", 0, ["POROS"], 0.4)
    assert np.allclose(df.loc[(slice(None), 0), "POROS"], 0.4)
    assert df[COLS].isna().all().all()


def test_ks_scenario_sampling_type_is_honoured_only_for_ks():
    from pyCATHY.DA import perturbate as P
    sc = {"per_name": ["ic", "Ks"], "sampling_type": ["lognormal", "lognormal"]}
    kw = dict(scenario_nom=1.0, scenario_mean=1.0, scenario_sd=0.5, clip_min=0.0, clip_max=9.0, NENS=4, per_type=None)
    ks = P.create_dict_entry("Ks", sc, 1, 0, **kw)
    ic = P.create_dict_entry("ic", sc, 0, 0, **kw)
    assert ks["sampling_type"] == "lognormal" and ic["sampling_type"] == "normal"      # ic: historical behaviour
    assert P.create_dict_entry("Ks", {"per_name": ["Ks"]}, 0, 0, **kw)["sampling_type"] == "normal"


def test_lognormal_ks_prior_is_centred_in_log_space_and_bounded():
    import numpy as np
    import pytest
    from pyCATHY.DA import perturbate as P
    np.random.seed(0)
    parm = dict(mean=2.1e-6, nominal=2.1e-6, sd=0.5, clip_min=2.1e-8, clip_max=2.1e-4, sampling_type="lognormal", per_type=None,
                ensemble_size=20000, transf_type=["log"], pert_control_name="material")
    out = P.VG_pert_rules({}, parm, "Ks0", 20000, 2.1e-6, 0.5, None, "lognormal")
    x = np.asarray(out[0] if isinstance(out, tuple) else out)
    assert x.min() >= 2.1e-8 and x.max() <= 2.1e-4 and (x > 0).all()
    assert np.exp(np.log(x).mean()) == pytest.approx(2.1e-6, rel=0.02)
    assert np.log(x).std() == pytest.approx(0.5, rel=0.05)


def test_layered_lognormal_ks_prior_and_normal_path_unchanged():
    import numpy as np
    import pytest
    from pyCATHY.DA import perturbate as P
    np.random.seed(0)
    x = P.perturbate_parm_by_layers("Ks0", 20000, 3, mean=2.1e-6, sd=0.5, nominal=2.1e-6, clip_min=2.1e-8,
                                    clip_max=2.1e-4, sampling_type="lognormal")
    assert x.shape == (3, 20000) and (x > 0).all() and x.min() >= 2.1e-8 and x.max() <= 2.1e-4
    assert np.exp(np.log(x).mean()) == pytest.approx(2.1e-6, rel=0.02) and np.log(x).std() == pytest.approx(0.5, rel=0.03)
    np.random.seed(0)
    a = P.perturbate_parm_by_layers("Ks0", 50, 3, mean=2.1e-6, sd=2.1e-6, nominal=2.1e-6, clip_min=2.1e-8, clip_max=2.1e-4)
    np.random.seed(0)
    b = P.perturbate_parm_by_layers("Ks0", 50, 3, mean=2.1e-6, sd=2.1e-6, nominal=2.1e-6, clip_min=2.1e-8, clip_max=2.1e-4,
                                    sampling_type="normal")
    assert np.array_equal(a, b)
