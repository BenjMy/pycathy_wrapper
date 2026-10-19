"""make_data_cov(diagonal=True): same analysis input, a fraction of the memory."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.observations import make_data_cov


def _dict_obs(n=40, cycles=(0.0, 3600.0, 7200.0), err=5.0, unit=None):
    rng = np.random.default_rng(0)

    def entry():
        e = {"data": pd.DataFrame({"rhoa": rng.uniform(10, 40, n), "k": np.full(n, 30.0),
                                   "valid": np.ones(n)}),
             "data_err": err, "data_type": "ERT", "units": "Ohm", "mesh_nodes": None,
             "sensor_name": "ERT"}
        if unit:
            e["data_err_unit"] = unit
        return e

    return {t: {"ERT": entry()} for t in cycles}


@pytest.mark.parametrize("unit", [None, "percent"])
def test_diagonal_storage_rebuilds_exactly_the_dense_matrix(unit):
    d = _dict_obs(unit=unit)
    _, _, dense = make_data_cov(SimpleNamespace(), d)
    _, _, diag = make_data_cov(SimpleNamespace(), d, diagonal=True)
    assert [np.ndim(x) for x in dense] == [2, 2, 2] and [np.ndim(x) for x in diag] == [1, 1, 1]
    for cycle, (a, b) in enumerate(zip(dense, diag)):
        da = SimpleNamespace(stacked_data_cov=[b] * 3, count_DA_cycle=cycle)
        np.testing.assert_array_equal(DA._data_cov_of_cycle(da, 40), a)  # identical, bit for bit


def test_dense_is_still_the_default():
    _, _, stacked = make_data_cov(SimpleNamespace(), _dict_obs())
    assert all(np.shape(x) == (40, 40) for x in stacked)
    da = SimpleNamespace(stacked_data_cov=stacked, count_DA_cycle=1)
    assert DA._data_cov_of_cycle(da, 40) is stacked[1]  # untouched: no copy, no conversion


def test_memory_is_linear_in_the_number_of_data():
    d = _dict_obs(n=500, cycles=tuple(float(i) for i in range(10)))
    _, _, dense = make_data_cov(SimpleNamespace(), d)
    _, _, diag = make_data_cov(SimpleNamespace(), d, diagonal=True)
    dense_bytes, diag_bytes = sum(x.nbytes for x in dense), sum(x.nbytes for x in diag)
    assert dense_bytes == 10 * 500 * 500 * 8 and diag_bytes == 10 * 500 * 8
    assert dense_bytes / diag_bytes == 500


def test_size_mismatch_is_still_reported_for_a_diagonal():
    da = SimpleNamespace(stacked_data_cov=[np.ones(10)], count_DA_cycle=0)
    with pytest.raises(ValueError, match="Wrong stacked_data_cov shape"):
        DA._data_cov_of_cycle(da, 12)
