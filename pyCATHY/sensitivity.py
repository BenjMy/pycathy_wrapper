#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dependency to the SAlib python package
Iwanaga, T., Usher, W., & Herman, J. (2022). Toward SALib 2.0: Advancing the accessibility and interpretability of global sensitivity analyses. Socio-Environmental Systems Modelling, 4, 18155. doi:10.18174/sesmo.18155
Herman, J. and Usher, W. (2017) SALib: An open-source Python library for sensitivity analysis. Journal of Open Source Software, 2(9). doi:10.21105/joss.00097

"""

import matplotlib.pyplot as plt
import numpy as np
# from SALib.analyze import morris as ma
from SALib.plotting import morris as mp
from SALib.sample import morris as ms


def define_Morris():
    """
    morris_problem = {
    # There are six variables
    'num_vars': 2,
    # These are their names
    'names': ['PERMX', 'POROS'
              ],
    # Plausible ranges over which we'll move the variables
    'bounds': [[1e-5, 1e-3], # Ks
                [0.4, 0.6], # porosity
              ],
    # I don't want to group any of these variables together
    'groups': None
    }
    """
    morris_problem = {}
    return morris_problem


def samples_generation(number_of_trajectories=24, num_levels=4):
    samples = ms.sample(morris_problem, number_of_trajectories, num_levels)
    return samples


def df_Morris(morris_problem, ref_parms):
    df_samples = pd.DataFrame(sample, columns=morris_problem["names"])
    df_samples.index.name = "sample"

    for p in morris_problem["names"]:
        df_samples["dev_" + p] = 1e2 * (df_samples[p] - ref_parms[p]) / SPP[p]
    return


def _prepare_CATHY_folders(prj_name):
    pathexe_list = []
    for ii in range(0, len(sample)):
        path_exe = os.path.join(prj_name + "_sensitivity", "sample" + str(ii + 1))
        pathexe_list.append(path_exe)
        # print(os.path.join(prj_name,prj_name + '_sensitivity/sample' + str(ii+1)))
        if os.path.exists(
            os.path.join(prj_name + "_sensitivity", "sample" + str(ii + 1))
        ):
            continue
        else:
            shutil.copytree(
                prj_name,
                os.path.join(prj_name + "_sensitivity", "sample" + str(ii + 1)),
            )
    return pathexe_list


def _update_CATHY_inputs(samples):

    for ii in range(0, len(samples)):
        print(str(ii + 1) + "/" + str(len(samples)))
        print(os.getcwd())

    PERMX = PERMY = PERMZ = sample[ii, 0]
    POROS = sample[ii, 1]
    SoilPhysProp = {
        "PERMX": PERMX,
        "PERMY": PERMY,
        "PERMZ": PERMZ,
        "ELSTOR": SPP["ELSTOR"],
        "POROS": POROS,
        "VGNCELL": SPP["VGNCELL"],
        "VGRMCCELL": SPP["VGRMCCELL"],
        "VGPSATCELL": SPP["VGPSATCELL"],
    }

    simu.update_soil(SPP=SoilPhysProp, path=pathexe_list[ii] + "/input/")
    pass


def run_sensitivity_analysis(prj_name):

    pathexe_list = _prepare_CATHY_folders(prj_name)
    _update_CATHY_inputs(samples)

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        result = pool.map(subprocess_run_multi, pathexe_list)

    simu_ensemble = np.zeros((len(obs_data), len(sample)))
    for ii in range(0, len(sample)):
        dpsi = simu.read_outputs("psi", path=pathexe_list[ii] + "/output/")
        dsw = simu.read_outputs("sw", path=pathexe_list[ii] + "/output/")
        simu_ensemble[:, ii] = np.hstack(dsw)

    return simu_ensemble


def analysis_Morris(morris_problem, obs_data, simu_ensemble, samples, err=0.025):

    rmse = np.zeros((1, len(samples)))
    for ii in range(0, len(samples)):
        rmse[0, ii] = err_weighted_rmse(
            simu_ensemble[:, ii], obs_data, err * obs_data
        )  # assume 2.5% noise in the data

    Si = ma.analyze(morris_problem, samples, rmse, print_to_console=True)

    return Si


def print_Morris_out(morris_problem):

    print(
        "{:20s} {:>7s} {:>7s} {:>7s}".format(
            "Name", "mean(EE)", "mean(|EE|)", "std(EE)"
        )
    )
    for name, s1, st, mean in zip(
        morris_problem["names"], Si["mu"], Si["mu_star"], Si["sigma"]
    ):
        print("{:20s} {:=7.3f} {:=7.3f} {:=7.3f}".format(name, s1, st, mean))

    pass


def err_weighted_rmse(sim, obs, noise):
    y = np.divide(sim - obs, noise)  # weighted data misfit
    y = np.sqrt(np.inner(y, y))
    return y


def plot_Morris():

    fig, ax = plt.subplots()
    ax.scatter(Si["mu_star"], Si["sigma"])
    # ax.plot(Si['mu_star'],2*Si['sigma']/np.sqrt(number_of_trajectories),'--',alpha=0.5)
    # ax.plot(np.array([0,Si['mu_star'][0]]),2*np.array([0,Si['sigma'][0]/np.sqrt(number_of_trajectories)]),'--',alpha=0.5)

    plt.title("Distribution of Elementary effects")
    plt.xlabel("mean(|EE|)")
    plt.ylabel("std($EE$)")
    for i, txt in enumerate(Si["names"]):
        ax.annotate(txt, (Si["mu_star"][i], Si["sigma"][i]))

    pass
