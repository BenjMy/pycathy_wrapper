"""
Test Data Assimlitation

This example shows long-term assimilation of tensiometer data to invert for the hydraulic condutivity. 
The data are real and come from a study in Noble (USA). 
The hydraulic simulation runs a RWU for a single plant for a volume of 1m3.
Soil is homogeneous in properties. 
Plant root water uptake at the center of the plot. 
Input file were previously created, see CATHY inputs example for more details
"""
import os
from collections import OrderedDict

import numpy as np
import pandas as pd

from pyCATHY.DA import perturbate
from pyCATHY.DA.cathy_DA import DA

dirName = "simu_test_DA"
prj_name = "tensiometers_DA"


# Load existing project object and fetch observation data
# ---------------------------------------------------------
# input file were previously created, see CATHY inputs example for more details

simu_test_DA = DA(dirName=dirName, prj_name=prj_name + "_DA", notebook=False)  #

data = pd.read_csv("./doc_test_data/doc_test_point_data.csv")
data_ass_time_s = pd.read_csv(
    "./doc_test_data/doc_test_point_data_assimilation_times.csv"
)


# Define a scenario
# -----------------
# A dictionnary describing how model parameters are pertubarted
# and what observation to assimilate
# for more information on how the parameters are perturbated: DA.perturbate_parm.__doc__
scenario = {
    "per_type": [None, "multiplicative"],
    "per_name": ["ic", "ks"],
    "per_nom": [-1.5, 1.880e-04],
    "per_mean": [-1.5, 1.880e-04],
    "per_sigma": [0.75, 1.75],
    "transf_type": [None, None],
    "listUpdateParm": ["St. var.", "ks"],
    "listAssimilatedObs": ["tensio"],
}

# Define simulation parameters
# ----------------------------
# Some more mandatory simulation parameters
ensemble_size = 32
DA_type = "enkf_Evensen2009_Sakov"
open_loop_run = True

# simu.update_parm(TMAX=1,IPRT1=3)
# simu.run_processor(IPRT1=3,verbose=True)


# Import data to assimilate
# -------------------------

data_err = 5
typ = "tensiometer"
colname = "a"
teros_1_depth = 7 * 0.0254
teros_1_mesh_node_pos = simu_test_DA.find_nearest_node([14.25, 1, -teros_1_depth])

for i, tt in enumerate(data_ass_time_s):

    dict_obs = simu_test_DA.read_observations(
        data[colname].iloc[i],
        typ,
        data_err,
        show=False,
        tA=tt,
        mesh_nodes=teros_1_mesh_node_pos,
        datetime=data["datetime"].iloc[i],
    )

sorted_dict_obs = OrderedDict(sorted(dict_obs.items()))

obs_df = dict_obs.dictObs_2pd()
print(obs_df)  # take a look at the data imported

# Make covariance matrice
# -----------------------
# cov_tensio = [  err_tensio      0               0
#                 0               err_tensio.     0
#                 0               0               err_tensio
#             ]

data_cov_diag = np.zeros([len(sorted_dict_obs), len(sorted_dict_obs)])
np.fill_diagonal(data_cov_diag, data_err)

simu_test_DA.prepare_observations(stacked_data_cov=[data_cov_diag])


# Synchronise time observation with atmbc times
# ---------------------------------------------
# The observation initial time is shifted compared to the first atmbc time

simu_test_DA.resynchronise_times(
    sorted_dict_obs,
    simu_test_DA.atmbc["time"],
)

# pertubate parameters
# -------------------------

list_pert = perturbate.perturbate(simu_test_DA, scenario, ensemble_size)

for dp in list_pert:
    # need to call perturbate_var as many times as variable to perturbate
    # return a dict merging all variable perturbate to parse into prepare_DA
    parm_per = DA.perturbate_parm(
        parm=dp,
        type_parm=dp["type_parm"],  # can also be VAN GENUCHTEN PARAMETERS
        mean=dp["mean"],
        sd=dp["sd"],
        sampling_type=dp["sampling_type"],
        ensemble_size=dp["ensemble_size"],  # size of the ensemble
        per_type=dp["per_type"],
        savefig=os.path.join(
            simu_test_DA.project_name, simu_test_DA.project_name + dp["savefig"]
        ),
    )


# run simulation + DA
# -------------------
DA.run_DA_sequential(
    parallel=True,
    dict_obs=sorted_dict_obs,
    list_assimilated_obs=scenario["listAssimilatedObs"],  # default
    list_update_parm=scenario["listUpdateParm"],
    DA_type=DA_type,  # default
    dict_parm_pert=parm_per,
    open_loop_run=open_loop_run,
)
