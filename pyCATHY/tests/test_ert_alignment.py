import numpy as np
import pytest
from pyCATHY.ERT.alignment import (ERTAlignmentError, bounds_from_nodes, check_ert_cathy_alignment,
                                   check_k_sign, suggest_frame, transform_points)

CB = np.array([0, 5, 0, 5, -2, 0.001])


def _grid(zsign=1, offset=(0.0, 0.0, 0.0)):
    """Electrode line on a 5x5x2 mesh in a raw frame; mesh nodes include the electrodes."""
    xs, ys = np.meshgrid(np.linspace(1, 4, 7), np.linspace(1, 4, 4))
    e = np.c_[xs.ravel(), ys.ravel(), np.zeros(xs.size)]
    zs = np.linspace(0, -2, 5)
    mesh = np.array([[x, y, z] for x in np.linspace(0, 5, 11) for y in np.linspace(0, 5, 11) for z in zs])
    mesh = np.vstack([mesh, e])
    f = lambda p: p * np.array([1, 1, zsign]) + np.array(offset)
    return f(e), f(mesh)


def test_aligned_passes():
    e, m = _grid()
    e, m = e + [0, 0, 0.001], m + [0, 0, 0.001]
    assert check_ert_cathy_alignment(e, CB, m, on_fail='raise')['ok']


def test_z_down_is_detected_and_suggested():
    e, m = _grid(zsign=-1)
    assert not check_ert_cathy_alignment(e, CB, m, on_fail='ignore')['ok']
    signs, shift, ov = suggest_frame(e, m, CB)
    assert signs[2] == -1 and ov > 0.95
    e2, m2 = transform_points(e, signs, shift), transform_points(m, signs, shift)
    assert check_ert_cathy_alignment(e2, CB, m2, on_fail='raise')['ok']


def test_electrodes_not_moved_with_mesh_raise():
    e, m = _grid()
    with pytest.raises(ERTAlignmentError):
        check_ert_cathy_alignment(e + [0.3, 0, 0.001], CB, m + [0, 0, 0.001], on_fail='raise')


def test_k_sign():
    assert check_k_sign(np.ones(50)) == 0
    with pytest.raises(ERTAlignmentError):
        check_k_sign(np.r_[np.ones(50), -np.ones(50)])
