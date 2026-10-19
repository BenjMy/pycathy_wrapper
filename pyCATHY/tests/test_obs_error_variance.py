import numpy as np
import pytest

from pyCATHY.DA.observations import obs_error_variance


def test_absolute_is_default():
    d = np.array([10.0, 100.0])
    np.testing.assert_allclose(obs_error_variance(d, 5), [25.0, 25.0])


def test_percent_scales_with_data():
    d = np.array([15.0, 100.0, -20.0])
    np.testing.assert_allclose(obs_error_variance(d, 5, "percent"),
                               (0.05 * np.abs(d)) ** 2)


def test_unknown_unit_raises():
    with pytest.raises(ValueError):
        obs_error_variance(np.ones(2), 5, "ohm")
