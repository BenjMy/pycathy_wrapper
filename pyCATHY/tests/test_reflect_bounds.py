import numpy as np
from pyCATHY.DA.param_tables import _reflect_into_bounds


def test_reflect_two_sided():
    out = _reflect_into_bounds(np.array([-0.2, 0.5, 1.3, 2.5, np.nan]), 0.0, 1.0)
    assert np.allclose(out[:4], [0.2, 0.5, 0.7, 0.5]) and np.isnan(out[4])


def test_reflect_one_sided_and_none():
    assert np.allclose(_reflect_into_bounds(np.array([-1.0, 2.0]), 0.0, None), [1.0, 2.0])
    assert np.allclose(_reflect_into_bounds(np.array([3.0, -1.0]), None, 0.0), [-3.0, -1.0])
    assert np.allclose(_reflect_into_bounds(np.array([5.0]), None, None), [5.0])


def test_reflect_keeps_spread_unlike_clip():
    v = np.array([-0.1, -0.2, -0.3])
    assert len(set(np.clip(v, 0, 1))) == 1
    assert len(set(np.round(_reflect_into_bounds(v, 0.0, 1.0), 9))) == 3
