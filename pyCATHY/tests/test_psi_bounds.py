import numpy as np
import pytest
from pyCATHY.DA.analysis import DAConfig
from pyCATHY.DA.param_tables import _clip_pressure_heads


def test_clip_two_sided_counts_and_keeps_nan():
    psi = np.array([-2000.0, -5.0, 0.3, np.nan])
    out, n = _clip_pressure_heads(psi, (-1e3, 0.0))
    assert n == 2
    np.testing.assert_array_equal(out[:3], [-1000.0, -5.0, 0.0])
    assert np.isnan(out[3])


def test_open_side_and_none():
    psi = np.array([-5.0, 3.0])
    assert _clip_pressure_heads(psi, (None, 0.0))[0].tolist() == [-5.0, 0.0]
    out, n = _clip_pressure_heads(psi, None)
    assert n == 0 and out.tolist() == [-5.0, 3.0]


def test_config_validation_still_raises():
    with pytest.raises(ValueError):
        DAConfig(psi_bounds=(0.0, -1.0))
    DAConfig(psi_bounds=(-1e3, 0.0))
