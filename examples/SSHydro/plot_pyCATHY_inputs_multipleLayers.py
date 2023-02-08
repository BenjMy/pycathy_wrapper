"""
Input plots
===========

This example shows how to use pyCATHY object to plot inputs of the hydrological model.

*Estimated time to run the notebook = 5min*

"""

# map_prop_veg ?
# map_prop2zone



#%% Import packages

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import numpy as np

#%% add CATHY object

path2prj ='weil_exemple_inputs_multipleLayers' # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj)


#%%
simu.update_parm()
simu.update_zone()
simu.update_soil()
df_soil, _ = simu.read_inputs('soil')
df = simu.read_inputs('soil')
simu.show_input(prop='soil', yprop='PERMX')

#%% Add layers

import matplotlib.pyplot as plt
import pandas as pd


#%%

get_layer_depths(simu)


layers_names = ['sub','bedrock']
layers_depths = [[0,1],[1,3]]
propertie_names = simu.soil_SPP['SPP_map'].keys()
p_subsurface_default = [1e-6,1e-6,1e-6,1.00E-05,0.48,1.914e+00,1.296e-01,1.240e+00]

SPP_perlayer = []
for li in range(len(layers_names)):
    SPP_p = []
    for i, p in enumerate(propertie_names):
        if li==1:
            SPP_p.append(p_subsurface_default[i]+1)
        else:
            SPP_p.append(p_subsurface_default[i])
    SPP_perlayer.append(SPP_p)


#%%
SoilPhysProp, SPP_map_dict = het_soil_layers_mapping_generic(simu,
                                                            propertie_names,
                                                            SPP_perlayer,
                                                            layers_names,
                                                            layers_depths,
                                                            )

#%%

simu.update_soil(IVGHU=0,
                  SPP_map=SPP_map_dict,
                  soil_heteregeneous_in_z=True,
                  # FP_map=FP_map,
                  )       
    


#%%
simu.show_input(prop='soil', yprop='PERMX', layer_nb=14)

        