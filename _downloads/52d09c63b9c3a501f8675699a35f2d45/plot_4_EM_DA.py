"""
Read EM dataset observations & assimilate
===========================================

The notebook illustrate how to read EM sensors dataset to be prepare for DA

*Estimated time to run the notebook = 2min*

"""
import numpy as np
from pyCATHY.DA.cathy_DA import DA
import pandas as pd
import matplotlib.pyplot as plt
from pyCATHY.DA.cathy_DA import DA, dictObs_2pd
from pyCATHY.DA.observations import read_observations, prepare_observations, make_data_cov
from pathlib import Path
import pyvista as pv

#%% Create a CATHY project
# -----------------------
simuWithDA = DA(
                dirName='./DA_with_swc',
                prj_name='import_SMC',
                notebook=True,
                )

#%% Archie transformation



ER_converted_ti, df_Archie, sw_nodes =   Archie.SW_2_ER0_DA(
                                        project_name,
                                        Archie_parms,
                                        POROS_mesh_nodes_ensi,
                                        EM_meta_dict,
                                        path_fwd_CATHY,
                                        DA_cnb=count_DA_cycle,  # kwargs
                                        Ens_nbi=ens_nb,  # kwargs
                                        # savefig=savefig,  # kwargs
                                        noise_level=EM_meta_dict["data_err"],  # kwargs
                                    )
df_Archie["OL"] = np.ones(len(df_Archie["time"])) * False
EC_converted_ti_mS_m = (1.0 / ER_converted_ti) * 1000

print("EC_converted_ti_mS_m stats:")
print(f"  min = {EC_converted_ti_mS_m.min():.2f} mS/m")
print(f"  max = {EC_converted_ti_mS_m.max():.2f} mS/m")

print('Build forward EM model')
depths, conds, xy_coords = build_forward_profiles(EC_converted_ti_mS_m,
                                                  grid3d['mesh3d_nodes'],
                                                  var="EC")

from emagpy import Problem
k = Problem()
k.setModels(depths,conds)

print('Forward EM model')
# dfsFSeq = k.forward(forwardModel='FSeq', coils=coils, noise=5)
EM_fwd_model_array = k.forward(forwardModel='FSlin',
                     coils=EM_meta_dict["coils"],
                     noise=EM_meta_dict["fwdEMnoise"]
                     )
# np.shape(dfsFSlin)
EM_fwd_model_2darray = np.vstack(EM_fwd_model_array)
EM_fwd_model_1darray = np.hstack(EM_fwd_model_2darray)
