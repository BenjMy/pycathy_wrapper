"""
Sensitivity analysis
=====================

Before run a Data Assimilation it is often necessary to evaluate the sensitivity of the model parameters with respect to a given scenario.
In this example, we use the Weil et al dataset and generate 24 possible trajectories varying PERMX and POROS parameters respectively the hydraulic
conductivity and the porosity of the soil. 
 
*Estimated time to run the notebook = 5min*

"""
import multiprocessing
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SALib.analyze import morris as ma
from SALib.plotting import morris as mp
from SALib.sample import morris as ms

from pyCATHY import cathy_tools
from pyCATHY.cathy_tools import subprocess_run_multi

#%% Create an observation scenario  and run the hydrological modelling
prj_name = "test0"
path2prj = "weil_exemple_sensitivityAnalysis"  # add your local path here
simu = cathy_tools.CATHY(
    dirName=path2prj, prj_name=prj_name, clear_src=False, clear_outputs=True
)
simu.run_preprocessor(verbose=False)
simu.run_processor(verbose=True)

SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)
dpsi = simu.read_outputs("psi")
dsw, _ = simu.read_outputs("sw")

obs_data = np.hstack(dsw)
len(obs_data)

#%% The Morris Problem
# number of variables, their names, plausible range and if you want to group or not the variable

morris_problem = {
    # There are six variables
    "num_vars": 2,
    # These are their names
    "names": ["PERMX", "POROS"],
    # Plausible ranges over which we'll move the variables
    "bounds": [
        [1e-5, 1e-3],  # Ks
        [0.4, 0.6],  # porosity
    ],
    # I don't want to group any of these variables together
    "groups": None,
}

#%% Sampling and plot

number_of_trajectories = 2
sample = ms.sample(morris_problem, number_of_trajectories, num_levels=4)
# sample = saltelli.sample(morris_problem, number_of_trajectories, num_levels=4)
# sample = saltelli.sample(problem, 1024)
df_sample = pd.DataFrame(sample, columns=morris_problem["names"])
df_sample.index.name = "sample"
for p in morris_problem["names"]:
    df_sample["dev_" + p] = 1e2 * (df_sample[p] - SPP_map[p]) / SPP_map[p]
fig, ax = plt.subplots()
mp.sample_histograms(fig, sample, morris_problem)

#%% creating subfolder for all the trajectories

pathexe_list = []
simu_ensemble = np.zeros((len(obs_data), len(sample)))
for ii in range(0, len(sample)):
    path_exe = os.path.join(
        simu.workdir, prj_name + "_sensitivity", "sample" + str(ii + 1)
    )
    pathexe_list.append(path_exe)
    if os.path.exists(
        os.path.join(simu.workdir, prj_name + "_sensitivity", "sample" + str(ii + 1))
    ):
        continue
    else:
        shutil.copytree(
            prj_name,
            os.path.join(
                simu.workdir, prj_name + "_sensitivity", "sample" + str(ii + 1)
            ),
        )
#%% mapping soil physical properties versus trajectory

for ii in range(0, len(sample)):
    PERMX = PERMY = PERMZ = sample[ii, 0]
    POROS = sample[ii, 1]
    SoilPhysProp = {
        "PERMX": PERMX,
        "PERMY": PERMY,
        "PERMZ": PERMZ,
        "ELSTOR": SPP_map["ELSTOR"],
        "POROS": POROS,
        "VGNCELL": SPP_map["VGNCELL"],
        "VGRMCCELL": SPP_map["VGRMCCELL"],
        "VGPSATCELL": SPP_map["VGPSATCELL"],
    }

    simu.update_soil(SPP_map=SoilPhysProp, path=pathexe_list[ii] + "/input/")

#%% running all the trajectories

with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
    result = pool.map(subprocess_run_multi, pathexe_list)

#%% read results and stack into a matrice

simu_ensemble = np.zeros((len(obs_data), len(sample)))
for ii in range(0, len(sample)):
    dpsi = simu.read_outputs("psi", path=pathexe_list[ii] + "/output/")
    dsw, _ = simu.read_outputs("sw", path=pathexe_list[ii] + "/output/")
    simu_ensemble[:, ii] = np.hstack(dsw)

# Perform the sensitivity analysis using the model output
# Specify which column of the output file to analyze (zero-indexed)
# Si = ma.analyze(
#     morris_problem,
#     sample,
#     rmse,
#     conf_level=0.95,
#     print_to_console=True,
# )
# Returns a dictionary with keys 'mu', 'mu_star', 'sigma', and 'mu_star_conf'
# e.g. Si['mu_star'] contains the mu* value for each parameter, in the
# same order as the parameter file

from SALib.plotting.morris import (
    covariance_plot,
    horizontal_bar_plot,
    sample_histograms,
)

# fig, (ax1, ax2) = plt.subplots(1, 2)
# horizontal_bar_plot(ax1, Si, {}, sortby="mu_star", unit=r"tCO$_2$/year")
# covariance_plot(ax2, Si, {}, unit=r"tCO$_2$/year")

# fig2 = plt.figure()
# sample_histograms(fig2, sample, morris_problem, {"color": "y"})
# plt.show()


#%% Define an objective function: here I use the error weighted rmse


def err_weighted_rmse(sim, obs, noise):
    y = np.divide(sim - obs, noise)  # weighted data misfit
    y = np.sqrt(np.inner(y, y))
    return y


#%% Compute the weighted data misfit y
rmse = np.zeros((1, len(sample)))
for ii in range(0, len(sample)):
    rmse[0, ii] = err_weighted_rmse(
        simu_ensemble[:, ii], obs_data, 0.025 * obs_data
    )  # assume 2.5% noise in the data

#%%
Si = ma.analyze(morris_problem, sample, rmse, print_to_console=True)

print("{:20s} {:>7s} {:>7s} {:>7s}".format("Name", "mean(EE)", "mean(|EE|)", "std(EE)"))
for name, s1, st, mean in zip(
    morris_problem["names"], Si["mu"], Si["mu_star"], Si["sigma"]
):
    print("{:20s} {:=7.3f} {:=7.3f} {:=7.3f}".format(name, s1, st, mean))

# total_Si, first_Si, second_Si = Si.to_df()
#%% Plot covariance

fig, ax = plt.subplots()
mp.covariance_plot(ax, Si)

#%% Plot Distribution of Elementary effects
# The higher mean |EE|, the more important factor
# line within the dashed envelope means nonlinear or interaction effects dominant

fig, ax = plt.subplots()
ax.scatter(Si["mu_star"], Si["sigma"])
# ax.plot(Si['mu_star'],2*Si['sigma']/np.sqrt(number_of_trajectories),'--',alpha=0.5)
# ax.plot(np.array([0,Si['mu_star'][0]]),2*np.array([0,Si['sigma'][0]/np.sqrt(number_of_trajectories)]),'--',alpha=0.5)

plt.title("Distribution of Elementary effects")
plt.xlabel("mean(|EE|)")
plt.ylabel("std($EE$)")
for i, txt in enumerate(Si["names"]):
    ax.annotate(txt, (Si["mu_star"][i], Si["sigma"][i]))
