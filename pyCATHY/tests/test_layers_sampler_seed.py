import warnings
import numpy as np
from pyCATHY.DA.perturbate import perturbate_parm_by_layers


def draw(**kw):
    return perturbate_parm_by_layers("n_VG", 64, 3, mean=1.5, sd=0.2, per_type=None, nominal=1.5, **kw)


def test_reproducible_with_global_seed_and_no_sobol_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        np.random.seed(7); a = draw()
        np.random.seed(7); b = draw()
    assert a.shape == (3, 64) and np.array_equal(a, b)


def test_statistics_and_bounds():
    np.random.seed(1)
    x = draw()
    assert abs(x.mean() - 1.5) < 0.05 and abs(x.std() - 0.2) < 0.05
    y = draw(clip_min=1.4, clip_max=1.6)
    assert y.min() >= 1.4 and y.max() <= 1.6


def test_sobol_still_available():
    assert draw(sampler="sobol").shape == (3, 64)
