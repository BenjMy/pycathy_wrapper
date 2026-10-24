"""Diagnostic plots of an ES-MDA run (``esmda.ESMDAResult``).

Only matplotlib and numpy are needed (no CATHY output is read), so the plots also work on the
result of ``esmda.run_esmda`` with any forward model.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def _active_masks(result):
    """Members shown at each of the ``n_iter + 1`` parameter snapshots (prior ... posterior).

    The prior shows every member; snapshot ``i`` shows the members that were still active in
    update ``i - 1``; the posterior shows the members valid at the end.
    """
    n_snap = len(result.param_history)
    n_ens = result.param_history[0].shape[1]
    masks = [np.ones(n_ens, dtype=bool)]
    for i in range(1, n_snap):
        masks.append(result.valid.copy() if i == n_snap - 1 else result.valid_history[i - 1].copy())
    return masks


def _check_names(result, names):
    n_param = result.param.shape[0]
    if names is None:
        names = getattr(result, "param_names", None)
    if names is None:
        return [f"p{j}" for j in range(n_param)]
    names = list(names)
    if len(names) != n_param:
        raise ValueError(f"names has {len(names)} entries for {n_param} parameter(s).")
    return names


def plot_misfit(result, ax=None, logy=True):
    """Reduced chi-square of every forward run (about 1 when the data are fitted).

    The last point is the run with the final parameters when ``final_forward=True``.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3.2))
    runs = np.arange(len(result.misfit_history))
    ax.plot(runs, result.misfit_history, "o-", color="tab:blue")
    ax.axhline(1.0, color="grey", ls="--", lw=1, label="chi$^2$/n = 1")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("forward run (0 = prior)")
    ax.set_ylabel("misfit  $r^T R^{-1} r / n_{obs}$")
    ax.set_xticks(runs)
    ax.legend(frameon=False)
    ax.set_title("Data misfit")
    return ax


def plot_parameter_evolution(result, names=None, truth=None, axes=None):
    """Ensemble of every parameter at the prior and after each update (box plots).

    Parameters
    ----------
    names : list of str, optional
        Parameter names (one per row of ``result.param``).
    truth : dict or sequence, optional
        True values drawn as horizontal lines (``{name: value}`` or one value per parameter).
    axes : sequence of Axes, optional
        One axis per parameter; created when omitted.
    """
    names = _check_names(result, names)
    n_param = len(names)
    if truth is not None and not isinstance(truth, dict):
        if len(truth) != n_param:
            raise ValueError(f"truth has {len(truth)} entries for {n_param} parameter(s).")
        truth = dict(zip(names, truth))
    if axes is None:
        _, axes = plt.subplots(n_param, 1, figsize=(5, 2.4 * n_param), squeeze=False, sharex=True)
        axes = axes[:, 0]
    axes = np.atleast_1d(axes)
    if len(axes) != n_param:
        raise ValueError(f"axes has {len(axes)} entries for {n_param} parameter(s).")

    masks = _active_masks(result)
    n_snap = len(result.param_history)
    labels = ["prior"] + [str(i) for i in range(1, n_snap)]
    for j, (name, ax) in enumerate(zip(names, axes)):
        values = [result.param_history[i][j, masks[i]] for i in range(n_snap)]
        ax.boxplot(values, positions=range(n_snap), widths=0.5, showfliers=False)
        if truth is not None and name in truth:
            ax.axhline(truth[name], color="tab:red", ls="--", lw=1, label="truth")
            ax.legend(frameon=False)
        ax.set_ylabel(name)
        ax.set_xticks(range(n_snap))
        ax.set_xticklabels(labels)
    axes[-1].set_xlabel("update")
    return axes


def plot_esmda_diagnostics(result, names=None, truth=None, logy=True):
    """Misfit and parameter evolution in one figure. Returns the ``Figure``."""
    names = _check_names(result, names)
    n_param = len(names)
    fig = plt.figure(figsize=(6, 3.0 + 2.2 * n_param), constrained_layout=True)
    grid = fig.add_gridspec(1 + n_param, 1, height_ratios=[1.3] + [1.0] * n_param)
    plot_misfit(result, ax=fig.add_subplot(grid[0]), logy=logy)
    axes = [fig.add_subplot(grid[1 + j]) for j in range(n_param)]
    plot_parameter_evolution(result, names=names, truth=truth, axes=axes)
    return fig
