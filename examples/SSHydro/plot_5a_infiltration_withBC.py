#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update and show Boundary conditions
===================================

This example shows how to use pyCATHY object to update a 3d BC properties from a DEM and run the hydrological model.

1st config: Dirichlet runs: 
    - The outlet nodes at the bottom layer form a constant head boundary of zero pressure head 
    and the nodes above along the outlet face have a no-flow condition imposed
2nd config: return flow runs:
    - the entire outlet face is a no-flow boundary and water is allowed 
    to leave the system only by exfiltration (return flow) through the surface.



*Estimated time to run the notebook = 5min*

"""


# !! run preprocessor change the DEM shape !
# dtm_13 does not have the same shape anymore!

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as cplt
import pyvista as pv

#%% Init CATHY model
# ------------------------
path2prj = "./SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name="bc_noflow"
                         )
#%% Fetch and show initial DEM
dem_mat, str_hd_dem = in_CT.read_dem(
    os.path.join(simu.workdir, simu.project_name, "prepro/dem"),
    os.path.join(simu.workdir, simu.project_name, "prepro/dtm_13.val"),
)

maxdepth = 10
zb = np.geomspace(1e-1, maxdepth, num=15)
nstr=len(zb)
zr = [abs(zb[0]/maxdepth)]
zr.extend(list(abs(np.diff(zb)/maxdepth)))

dem = np.ones(np.shape(dem_mat))
dem[-1,-1] = 1-1e-3

simu.update_prepo_inputs(
    DEM=np.ones(np.shape(dem_mat)),
    nstr=nstr,
    zratio=zr,
    base=max(zb),
)
simu.create_mesh_vtk(verbose=True)

#%% Show results

# add nodes of interest
node_bot, node_bot_pos = simu.find_nearest_node([5,5,-9])
node_top, node_top_pos = simu.find_nearest_node([5,5,1])
node_leftxmin, node_leftxmin_pos = simu.find_nearest_node([0,5,-9/2])
NOI = [node_bot,node_top,node_leftxmin]
NOI_pos = [node_bot_pos,node_top_pos,node_leftxmin_pos]
NOI_labels = ['bot','top','left_xmin']
NOI_colors = ['b','orange','green']

pl = pv.Plotter(notebook=False)
mesh = pv.read(os.path.join(simu.workdir,
                     simu.project_name,
                     'vtk',simu.project_name + '.vtk'
                  )
        )

pl.add_mesh(mesh,
            style='wireframe',
            opacity=0.1,
            color='k'
            )
for i, n in enumerate(NOI):
    pl.add_points(NOI_pos[i][0],
                  color=NOI_colors[i],
                  label=NOI_labels[i],
                  point_size=30
                  )
pl.show_bounds()
pl.view_zx()
pl.show()



#%% We simulate here a constant infiltration of 5.5e-6 m/s during [0,12e3] seconds

df_atmbc = simu.read_inputs('atmbc')
time = [0,12e3,18e3]
simu.update_atmbc(   
                    HSPATM=1,
                    IETO=1,
                    time=time,
                    netValue=[5.5e-06, 0, 0]
                    )
simu.show_input('atmbc')

#%% Init boudary conditions
# This mesh is automatically created when updating with update_nansfneubc() 
simu.create_mesh_bounds_df(
                            'nansfdirbc', 
                            simu.grid3d["mesh3d_nodes"], 
                            time, 
                            )
print(simu.mesh_bound_cond_df.head())
 
#%% No flow boundary conditions 
simu.update_nansfneubc(no_flow=True)
simu.update_nansfdirbc(no_flow=True)
simu.update_sfbc(no_flow=True)

#%% Show Boundary Conitions at time 0
simu.show_bc(time=0)
meshbc = simu.mesh_bound_cond_df
cplt.plot_mesh_bounds('nansfdirbc',meshbc, time=0)

#%% Update initial conditons
simu.update_ic(INDP=0,
               IPOND=0,
               pressure_head_ini=-15
               )

#%% Run  hydrological model
simu.update_parm(DELTAT=1e3)
simu.run_processor(IPRT1=2, 
                    TRAFLAG=0, 
                    verbose=True
                    )

#%% Get saturation water results
df_sw, _ = simu.read_outputs('sw')

#%%
import copy
# We create a second CATHY object to compare when changing Neumann BC
simu_with_neubc = copy.copy(simu)

#%%
simu_with_neubc.update_nansfneubc(
                                time = time,
                                bot_bound=1.2e-6
                                )


#%%
simu_with_neubc.show_bc(time=0)
meshbc = simu_with_neubc.mesh_bound_cond_df
cplt.plot_mesh_bounds('nansfneubc',
                      meshbc, 
                      time=0
                      )

#%%
simu_with_neubc.run_processor(IPRT1=2, 
                              TRAFLAG=0, 
                              verbose=True
                    )
#%% Get saturation data

df_sw_neumBC, _ = simu_with_neubc.read_outputs('sw')
df_sw.head()

fig, axs = plt.subplots(2,1,
                        sharey=True
                        )

for i, n in enumerate(NOI):
    df_sw[n].plot(ax=axs[0],
                  label=NOI_labels[i],
                  c=NOI_colors[i]
                     )
    df_sw_neumBC[n].plot(ax=axs[1],
                  label=NOI_labels[i],
                  c=NOI_colors[i],
                     )
axs[0].legend('')
axs[1].legend(NOI_labels)
axs[1].set_xlabel('time (s)')
axs[1].set_ylabel('saturation (-)')
