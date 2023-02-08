#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 17:35:10 2023

@author: ben
"""

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
import pyCATHY.meshtools as mt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Init CATHY model
# ------------------------
path2prj ='../SSHydro/' # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj,prj_name='soil3d_het_slope')
    

#%% Plot again dem in 3d

simu.update_parm()
simu.update_veg_map()
simu.update_soil()    
    
fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop='dem', ax=ax)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()

# plt.savefig('../figs/weil_dem.png', 
#             dpi=400, bbox_inches='tight', pad_inches=0)

#%% This is creating grid3d file
# grid3d = utils_CATHY.create_mesh(simu)

#%% This creates a vtk file from grid3d
simu.create_mesh_vtk()
simu.update_zone()
zone3d = mt.zone3d(simu)

SPP_map =     {
                    'PERMX': [1e-4, 1e-4],
                    'PERMY': [1e-4, 1e-4],
                    'PERMZ': [1e-4, 1e-4],
                    'ELSTOR': [1e-05, 1e-05],
                    'POROS': [0.45, 0.55],
                    'VGNCELL': [1.46, 1.46],
                    'VGRMCCELL': [0.15, 0.15],
                    'VGPSATCELL': [0.03125, 0.03125]
            }


#%%

zone3d_top = []
zone3d_bot = []
for li in range(simu.dem_parameters["nstr"]):
    top, bot = mt.get_layer_depth(simu,li)
    print(top)
    zone3d_top.append(zone3d[li]*top)
    zone3d_bot.append(zone3d[li]*bot)

zone3d_top = np.array(zone3d_top)
zone3d_bot = np.array(zone3d_bot)

# Get coords dem
xcoord_dem, ycoord_dem, dem_mat = cplt.get_dem_coords(hapin=simu.hapin,
                                                        workdir=simu.workdir, 
                                                        project_name=simu.project_name,
                                                        )
dem_mat = np.flipud(dem_mat)

# Adjust altitude of dem
# dem_mat3d_layers = [dem_mat - zz for zz in zone3d_top] #-(zone3d_top-zone3d_bot)/2
dem_mat3d_layers = [dem_mat - zz for zz in zone3d_top] #-(zone3d_top-zone3d_bot)/2

np.shape(dem_mat3d_layers)
np.shape(zone3d_top)

zone3d_top


#%% Add flag zone
slope_bedrockvar = 0.2


# build a 2 grid
xgrid, ygrid = np.meshgrid(ycoord_dem,xcoord_dem)
grid_coords = [xgrid, ygrid]

zgrid_bedrock = (xgrid-min(xgrid[:,0]))*slope_bedrockvar
zgrid_bedrock = abs(zgrid_bedrock)

zone3d_topflag = []
for li in range(simu.dem_parameters["nstr"]):
    
    zone3d_topbot = zone3d_top #+ (zone3d_top-zone3d_bot)/2
    booltop = zone3d_topbot[li]<zgrid_bedrock
    zone3d_defaultflag = np.ones(np.shape(zone3d_top[0]))
    zone3d_defaultflag[booltop]=2
    zone3d_topflag.append(zone3d_defaultflag)

np.shape(booltop)
np.shape(zone3d_topflag)
zone3d_topflag = np.array(zone3d_topflag)

#%%


# zone3d_topflag_test[:,:,0:5]=2

simu.update_zone(zone3d_topflag[0])
simu.update_soil(SPP_map=SPP_map,
                 zone3d=zone3d_topflag,
                 )

#%%

simu.run_processor(IPRT1=2,verbose=True)


