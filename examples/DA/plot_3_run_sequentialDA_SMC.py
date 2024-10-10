"""
Read SMC sensors observations to assimilate
===========================================

The notebook illustrate how to read SMC sensors dataset to be prepare for DA

*Estimated time to run the notebook = 2min*

"""
import numpy as np
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.observations import make_data_cov
from pyCATHY.DA.cathy_DA import DA, dictObs_2pd
from pyCATHY.DA import perturbate
import pickle

#%% Create a CATHY project
# -----------------------
simuWithDA = DA(
                dirName='./DA_with_swc',
                prj_name='DA_SMC',
                notebook=True,
                )


#%% Set absolute error level abnd build dictionnary of observations
abs_data_err = 1e-1 # constant error does not vary with time
dict_obs = {} # initiate the dictionnary

with open('./DA_with_swc/obs_prepared_SMC.pkl', 'rb') as fp:
    dict_obs = pickle.load(fp)
data_measure_df = dictObs_2pd(dict_obs) 
    # data_measure_df.index
#%% Create observation covariance matrices for all assimilation time steps
# By default, there is no correlation between sensors
# Therefore, the covariance matrices are diagonal with the error values on the diagonals

_,_, stacked_data_cov = make_data_cov(
                                        simuWithDA,
                                        dict_obs,
                                        list_assimilated_obs = 'all',
                                        nb_assimilation_times=len(dict_obs)
                                        )
print(np.shape(stacked_data_cov))
simuWithDA.stacked_data_cov = stacked_data_cov

#%%
DEM, _ = simuWithDA.read_inputs('dem')
simuWithDA.DEM = DEM
simuWithDA.update_dem_parameters()
simuWithDA.update_veg_map()
simuWithDA.update_soil()

NENS = 5

# ZROOT
# -------------------
pert_nom_ZROOT = 1
pert_sigma_ZROOT = 0.35e-9
minZROOT = 0
maxZROOT = 2

scenario = {'per_type': [None],
             'per_name':['ZROOT'],
             'per_nom':[pert_nom_ZROOT],
             'per_mean':[pert_nom_ZROOT],    
             'per_sigma': [pert_sigma_ZROOT],
             'per_bounds': [
                            {'min':minZROOT,'max':maxZROOT}
                            ],
             'sampling_type': ['normal'],
             'transf_type':[None],
             'listUpdateParm': ['St. var.', 'ZROOT'],
             'listObAss': ['SMC'],
             }
            
scenario['per_name']  
                           
list_pert = perturbate.perturbate(simuWithDA, 
                                  scenario, 
                                  NENS
                                  )   

#%% Parameters perturbation
# stop
import os

var_per_dict_stacked = {}
for dp in list_pert:
    savefig = os.path.join(
                            simuWithDA.workdir,
                            simuWithDA.project_name,
                            simuWithDA.project_name + dp['savefig']
                            )
    np.random.seed(1)
    # need to call perturbate_var as many times as variable to perturbate
    # return a dict merging all variable perturbate to parse into prepare_DA
    var_per_dict_stacked = perturbate.perturbate_parm(
                                                    var_per_dict_stacked,
                                                    parm=dp, 
                                                    type_parm = dp['type_parm'], # can also be VAN GENUCHTEN PARAMETERS
                                                    mean =  dp['mean'],
                                                    sd =  dp['sd'],
                                                    sampling_type =  dp['sampling_type'],
                                                    ensemble_size =  dp['ensemble_size'], # size of the ensemble
                                                    per_type= dp['per_type'],
                                                    savefig=savefig
                                                    )
    
    
#%% Run assimilation
# f
# simuWithDA.parm
# simuWithDA.read_inputs('atmbc')
atmbc_times = data_measure_df.index.get_level_values(1).unique().to_list()
simuWithDA.update_atmbc(HSPATM=1,IETO=0,
                        time=atmbc_times, 
                        netValue=[0]*len(atmbc_times)
                        )

# simuWithDA.update_parm()
# simuWithDA.read_inputs('atmbc')


#%%
# simuWithDA.atmbc

# simuWithDA.run_DA_smooth(
#                           VTKF=2,
#                           TRAFLAG=0,
#                           dict_obs= dict_obs,
#                           list_assimilated_obs='all', # default
#                           list_parm2update= ['St. var.', 'ZROOT0'],
#                           DA_type='enkf_Evensen2009',
#                           dict_parm_pert=var_per_dict_stacked,
#                         )

#%%

# aa
simuWithDA.run_DA_sequential(
                              VTKF=2,
                              TRAFLAG=0,
                              dict_obs= dict_obs,
                              list_assimilated_obs='all', # default
                              list_parm2update= ['St. var.', 'ZROOT0'],
                              DA_type='enkf_Evensen2009',
                              dict_parm_pert=var_per_dict_stacked,
                            )
