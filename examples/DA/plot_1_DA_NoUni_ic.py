"""
DA with Random Initial Conditions on Soil Layers
============================================

This notebook demonstrates how to create an ensemble of models with varying initial conditions for different soil layers. By introducing randomness into the initial conditions, we can better simulate the natural variability in soil properties and improve the robustness of data assimilation (DA) processes.

*Estimated time to run the notebook = 2min*

"""

import os
import numpy as np
import pyvista as pv
from pyCATHY.DA.perturbate import perturbate_parm
from pyCATHY.DA import perturbate
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.plotters import cathy_plots as CTp


#%% Define scenario where initial conditions are perturbated by layers
# -----------------------
nlay = 6
scenario = {
            'per_name':['ic'],
            # 'per_type': [[None]*nlay],
            'per_type': [['multiplicative']*nlay],
            'per_nom':[[-15]*nlay],
            # 'per_mean':[[-15]*nlay],
            'per_mean':[[1]*nlay],
            'per_sigma': [[3.75]*nlay],
            'per_bounds': [[None]*nlay],
            'sampling_type': [['normal']*nlay],
            'transf_type':[[None]*nlay],
            'listUpdateParm': ['St. var.'],
            'pert_control_name': ['layers']*nlay
            }



#%% Create a CATHY project
# -----------------------
simuWithDA = DA(
                        dirName='./',
                        prj_name= 'DA_with_non_uniform_ic',
                        notebook=True,
                    )

# linear z depth
# ---------------
zb = np.linspace(0, 2, nlay)
nstr = len(zb)
zr = list(np.ones(len(zb))/nstr)

simuWithDA.update_prepo_inputs(
                                nstr=nstr,
                                zratio=zr,
                                base=max(zb),
                                )
#%% Add dem and soil/veg properties to the simuWithDA object
# ----------------------------------------------------------
simuWithDA.update_dem_parameters()
simuWithDA.update_veg_map()
simuWithDA.update_soil()

#%%  Perturbate ic according to the scenario
# ------------------------------------------

simuWithDA.NENS = 3
list_pert = perturbate.perturbate(
                                    simuWithDA,
                                    scenario,
                                    simuWithDA.NENS,
                                 )

var_per_dict_stacked = {}
for dp in list_pert:
    var_per_dict_stacked = perturbate_parm(
                                var_per_dict_stacked,
                                parm=dp,
                                type_parm = dp['type_parm'], # can also be VAN GENUCHTEN PARAMETERS
                                mean =  dp['mean'],
                                sd =  dp['sd'],
                                sampling_type =  dp['sampling_type'],
                                ensemble_size =  dp['ensemble_size'], # size of the ensemble
                                per_type= dp['per_type'],
                                nlayers = nlay,
                                savefig= os.path.join(simuWithDA.workdir,
                                                      simuWithDA.project_name,
                                                      simuWithDA.project_name + dp['savefig'])
                                )

#%% Update ensemble and plot/save ic vtk
# This in normally directly called when using run_DA_sequential()

simuWithDA.ens_valid =  list(np.arange(0, simuWithDA.NENS))
simuWithDA._create_subfolders_ensemble()
simuWithDA.apply_ensemble_updates(var_per_dict_stacked,
                                  var_per_dict_stacked.keys(),
                                  cycle_nb=0
                                  )
#%% Plot results
# -----------------------

pl = pv.Plotter(shape=(1,2))
for i, ensi in enumerate([1,3]):
    DApath = f'DA_Ensemble/cathy_{ensi}/vtk/'
    path = os.path.join(simuWithDA.workdir,
                        simuWithDA.project_name,
                        DApath,
                        simuWithDA.project_name + '.vtk'
                        )

    pl.subplot(0,i)
    CTp.show_vtk(path,
                 'ic_nodes',
                 ax=pl,
                 clim = [-25,-5],
                 #show_scalar_bar=True,
                 )
    _ = pl.add_legend('')
    pl.add_title(f'Ensemble nb:{ensi}')

pl.show()
