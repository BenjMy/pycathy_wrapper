"""ES-MDA diagnostic plots: structure of the figures, no display needed."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pyCATHY.DA.esmda import ESMDAConfig, run_esmda
from pyCATHY.plotters.esmda_plots import (
    _active_masks,
    plot_esmda_diagnostics,
    plot_misfit,
    plot_parameter_evolution,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _result(n_iter=3, fail_at=None, n_ens=30):
    rng = np.random.default_rng(0)
    G = rng.normal(size=(5, 2))
    d = G @ np.array([1.0, -0.5]) + 0.05 * rng.normal(size=5)
    calls = {"n": 0}

    def forward(p, valid):
        ok = np.ones(p.shape[1], bool)
        if fail_at is not None and calls["n"] == fail_at:
            ok[[2, 4]] = False
        calls["n"] += 1
        return G @ p, ok

    prior = rng.normal(0.0, 1.0, size=(2, n_ens))
    return run_esmda(prior, forward, d, np.full(5, 0.05**2), ESMDAConfig(n_iter=n_iter, seed=0))


def test_active_masks_follow_valid_history():
    res = _result(n_iter=3, fail_at=1)
    masks = _active_masks(res)
    assert len(masks) == 4
    assert masks[0].all()
    assert masks[1].all()  # update 0 had every member
    assert masks[2].sum() == 28 and masks[3].sum() == 28
    assert np.array_equal(masks[3], res.valid)


def test_plot_misfit_data_and_labels():
    res = _result()
    ax = plot_misfit(res)
    line = ax.get_lines()[0]
    assert np.allclose(line.get_ydata(), res.misfit_history)
    assert np.allclose(line.get_xdata(), np.arange(len(res.misfit_history)))
    assert ax.get_yscale() == "log"
    assert "misfit" in ax.get_ylabel()
    assert plot_misfit(res, logy=False).get_yscale() == "linear"


def test_plot_parameter_evolution_structure_and_truth():
    res = _result()
    truth = {"Ks": 1.0, "porosity": -0.5}
    axes = plot_parameter_evolution(res, names=["Ks", "porosity"], truth=truth)
    assert len(axes) == 2
    assert [a.get_ylabel() for a in axes] == ["Ks", "porosity"]
    # x axes are shared: tick labels are shown on the bottom axis
    labels = [t.get_text() for t in axes[-1].get_xticklabels()]
    assert labels == ["prior", "1", "2", "3"]
    for ax, name in zip(axes, ["Ks", "porosity"]):
        flat = [l for l in ax.get_lines() if len(l.get_xdata()) == 2 and l.get_ydata()[0] == l.get_ydata()[1]]
        assert any(np.isclose(l.get_ydata()[0], truth[name]) for l in flat)  # the truth line
        assert ax.get_legend() is not None


def test_truth_as_sequence_and_default_names():
    res = _result()
    axes = plot_parameter_evolution(res, truth=[1.0, -0.5])
    assert [a.get_ylabel() for a in axes] == ["p0", "p1"]


def test_validation_errors():
    res = _result()
    with pytest.raises(ValueError, match="names has"):
        plot_parameter_evolution(res, names=["only_one"])
    with pytest.raises(ValueError, match="truth has"):
        plot_parameter_evolution(res, truth=[1.0])
    fig, axes = plt.subplots(3, 1)
    with pytest.raises(ValueError, match="axes has"):
        plot_parameter_evolution(res, axes=axes)


def test_full_diagnostics_figure():
    res = _result(fail_at=1)
    fig = plot_esmda_diagnostics(res, names=["Ks", "porosity"], truth=[1.0, -0.5])
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) == 3  # misfit + 2 parameters
    assert fig.axes[0].get_title() == "Data misfit"


def test_plots_do_not_modify_result():
    res = _result()
    before = [p.copy() for p in res.param_history]
    plot_esmda_diagnostics(res)
    assert all(np.array_equal(a, b) for a, b in zip(before, res.param_history))


def test_plots_use_the_parameter_names_of_the_result():
    res = _result()
    assert res.param_names is None  # the generic driver does not know them
    res.param_names = ["ic", "Ks"]
    axes = plot_parameter_evolution(res)
    assert [a.get_ylabel() for a in axes] == ["ic", "Ks"]
    # explicit names still win
    axes = plot_parameter_evolution(res, names=["a", "b"])
    assert [a.get_ylabel() for a in axes] == ["a", "b"]
