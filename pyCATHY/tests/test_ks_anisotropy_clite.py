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
