"""
Sensitivity analysis
=====================

Before running a Data Assimilation it is often necessary to evaluate the
sensitivity of the model parameters with respect to a given scenario.
Here we use the Weil et al dataset and generate 24 trajectories varying
PERMX and POROS (hydraulic conductivity and porosity of the soil).

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

#%% Create an observation scenario and run the hydrological modelling
prj_name = "test0"
path2prj = "weil_exemple_sensitivityAnalysis"  # add your local path here

simu = cathy_tools.CATHY(dirName=path2prj, prj_name=prj_name)

simu.run_preprocessor(verbose=False)
simu.run_processor(
    IPRT1=2,
    DTMIN=1e-2,
    DTMAX=1e2,
    DELTAT=5,
    TRAFLAG=0,
    verbose=False,
)

SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)

dsw, _ = simu.read_outputs("sw")
obs_data = dsw.to_numpy()  # (n_times, n_nodes)
obs = obs_data[-1, :]      # last time step, (n_nodes,)
n_nodes = obs.size

#%% The Morris problem
morris_problem = {
    "num_vars": 2,
    "names": ["PERMX", "POROS"],
    "bounds": [
        [1e-5, 1e-3],  # Ks
        [0.4, 0.6],    # porosity
    ],
    "groups": None,
}

#%% Sampling and plot
# Morris gives N * (num_vars + 1) samples: 8 * (2 + 1) = 24
number_of_trajectories = 8
sample = ms.sample(morris_problem, number_of_trajectories, num_levels=4)

df_sample = pd.DataFrame(sample, columns=morris_problem["names"])
df_sample.index.name = "sample"

# SPP_map is (zone, layer)-indexed, dtype object: cast to float and
# collapse to one reference value (homogeneous soil)
for p in morris_problem["names"]:
    ref = SPP_map[p].astype(float).mean()
    df_sample["dev_" + p] = 1e2 * (df_sample[p] - ref) / ref

fig = plt.figure()
mp.sample_histograms(fig, sample, morris_problem)

#%% Create one sub-folder per trajectory
pathexe_list = []
for ii in range(len(sample)):
    path_exe = os.path.join(
        simu.workdir, prj_name + "_sensitivity", "sample" + str(ii + 1)
    )
    pathexe_list.append(path_exe)
    if os.path.exists(path_exe):
        continue
    shutil.copytree(os.path.join(simu.workdir, prj_name), path_exe)

#%% Map soil physical properties onto each trajectory
simu.update_veg_map()
for ii in range(len(sample)):
    # fresh defaults each iteration (avoids carry-over between samples)
    SoilPhysProp = simu.set_SOIL_defaults(SPP_map_default=True)
    SoilPhysProp[["PERMX", "PERMY", "PERMZ"]] = sample[ii, 0]
    SoilPhysProp["POROS"] = sample[ii, 1]

    simu.update_soil(SPP_map=SoilPhysProp, path=pathexe_list[ii] + "/input/")

#%% Run all the trajectories
with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
    result = pool.map(subprocess_run_multi, pathexe_list)

#%% Read results (last time step) -> (n_nodes, n_samples)
simu_ensemble = np.zeros((n_nodes, len(sample)))
for ii in range(len(sample)):
    dsw_ii, _ = simu.read_outputs("sw", path=pathexe_list[ii] + "/output/")
    simu_ensemble[:, ii] = dsw_ii.to_numpy()[-1, :]

#%% Objective function: error-weighted misfit (weighted L2 norm)
def err_weighted_rmse(sim, obs, noise):
    y = np.divide(sim - obs, noise)  # weighted data misfit
    return np.sqrt(np.inner(y, y))


#%% One scalar per trajectory -> 1D array (n_samples,), as required by SALib
noise = 0.025 * obs  # assume 2.5% noise in the data

rmse = np.array([
    err_weighted_rmse(simu_ensemble[:, ii], obs, noise)
    for ii in range(len(sample))
])
assert rmse.shape == (len(sample),)

print("rmse:", rmse)
print("max spread across samples:", np.ptp(simu_ensemble, axis=1).max())
print("return codes:", result)

# 1) did the soil files actually differ between samples?
for p in pathexe_list[:3]:
    print(p, os.listdir(os.path.join(p, "input")))
    print(open(os.path.join(p, "input", "soil")).read()[:300])  # adjust filename if needed

# 2) are the outputs different files / freshly written?
for p in pathexe_list[:3]:
    f = os.path.join(p, "output", "sw")
    print(f, os.path.exists(f), os.path.getmtime(f) if os.path.exists(f) else None)
    
#%% Morris analysis
Si = ma.analyze(morris_problem, sample, rmse, print_to_console=True)

print("{:20s} {:>10s} {:>10s} {:>10s}".format("Name", "mean(EE)", "mean(|EE|)", "std(EE)"))
for name, mu, mu_star, sigma in zip(
    morris_problem["names"], Si["mu"], Si["mu_star"], Si["sigma"]
):
    print("{:20s} {:10.3f} {:10.3f} {:10.3f}".format(name, mu, mu_star, sigma))

#%% Covariance plot
fig, ax = plt.subplots()
mp.covariance_plot(ax, Si)

#%% Distribution of elementary effects
# The higher mean|EE|, the more important the factor;
# a high std(EE) means nonlinear or interaction effects dominate
fig, ax = plt.subplots()
ax.scatter(Si["mu_star"], Si["sigma"])
plt.title("Distribution of Elementary effects")
plt.xlabel("mean(|EE|)")
plt.ylabel("std($EE$)")
for i, txt in enumerate(Si["names"]):
    ax.annotate(txt, (Si["mu_star"][i], Si["sigma"][i]))

plt.show()