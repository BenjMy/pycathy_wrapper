"""
ES-MDA on a synthetic problem
=============================

Ensemble Smoother with Multiple Data Assimilation (ES-MDA, Emerick & Reynolds, 2013) on a
toy forward model, **without running CATHY**. It shows the pieces you also use on a CATHY
ensemble (``DA.run_DA_esmda``): the alpha schedule, the forward-model callable, the result
object and the diagnostic plots.

The same data are assimilated ``n_iter`` times with an observation error inflated by
``alpha_i`` (``sum(1/alpha_i) = 1``); the forward model is re-run from t = 0 between
iterations.

*Estimated time to run the notebook = a few seconds*
"""
import matplotlib.pyplot as plt
import numpy as np

from pyCATHY.DA.esmda import ESMDAConfig, compute_alphas, run_esmda
from pyCATHY.DA.esmda_plots import plot_esmda_diagnostics

#%% Toy forward model
# -------------------
# A decaying signal ``amplitude * exp(-rate * t)`` observed at 12 times: non-linear in
# ``rate``, so several iterations help.

t_obs = np.linspace(0.2, 3.0, 12)


def g(param):
    amplitude, rate = param
    return amplitude[None, :] * np.exp(-rate[None, :] * t_obs[:, None])


truth = {"amplitude": 2.0, "rate": 0.8}
sigma = 0.03  # observation error (standard deviation)
rng = np.random.default_rng(42)
data = g(np.array([[truth["amplitude"]], [truth["rate"]]]))[:, 0] + sigma * rng.normal(size=t_obs.size)

#%% Prior ensemble
# ----------------
n_ens = 100
prior = np.vstack([rng.normal(1.4, 0.6, n_ens), rng.normal(1.1, 0.4, n_ens)])

#%% Alpha schedule
# ----------------
# ``constant``: alpha_i = n_iter. ``geometric``: large alphas first (small, safe updates),
# then smaller ones. Both satisfy sum(1/alpha) = 1.
for schedule in ("constant", "geometric"):
    alphas = compute_alphas(4, schedule)
    print(f"{schedule:10s} alphas = {np.round(alphas, 2)}   sum(1/alpha) = {np.sum(1 / alphas):.3f}")

#%% Run ES-MDA
# ------------
# ``forward(param, valid)`` runs the model for the members flagged in ``valid`` and
# returns the predictions and a flag per member telling which ones ran correctly.
# With CATHY this is a full simulation of every member; here it is a formula.


def forward(param, valid):
    return g(param), np.ones(param.shape[1], dtype=bool)


def project(param):
    """Physical bounds, applied after every update (CATHY: the scenario bounds)."""
    return np.clip(param, 0.05, 5.0)


config = ESMDAConfig(n_iter=4, alpha_schedule="geometric", seed=1)
result = run_esmda(
    prior,
    forward,
    data,
    np.full(t_obs.size, sigma**2),  # observation error variances
    config,
    project=project,
)

print("misfit (reduced chi-square) per forward run:", np.round(result.misfit_history, 2))
for name, row in zip(truth, result.param):
    print(f"{name:10s} posterior {row.mean():.3f} +/- {row.std():.3f}   (truth {truth[name]})")

#%% Diagnostics
# -------------
# Misfit per forward run (about 1 when the data are fitted within their error) and the
# ensemble of every parameter after each update.
plot_esmda_diagnostics(result, names=list(truth), truth=truth)
plt.show()

#%% Predictions before / after
# ----------------------------
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(t_obs, g(prior), color="tab:blue", alpha=0.08)
ax.plot(t_obs, g(result.param), color="tab:orange", alpha=0.15)
ax.plot([], [], color="tab:blue", label="prior ensemble")
ax.plot([], [], color="tab:orange", label="posterior ensemble")
ax.errorbar(t_obs, data, yerr=sigma, fmt="k.", label="data")
ax.set_xlabel("time")
ax.set_ylabel("observation")
ax.legend(frameon=False)
plt.show()
