#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:39:59 2024

@author: z0272571a
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

#%% Init CATHY model
# ------------------------
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
			prj_name="meshing_from_weill_test2"
				)

rootpath = os.path.join(simu.workdir + simu.project_name)


#%% Fetch and show initial DEM
# the dimension of the mesh is squared (20,20)


dem_mat, str_hd_dem = in_CT.read_dem(
    os.path.join(simu.workdir, simu.project_name, "prepro/dem"),
    os.path.join(simu.workdir, simu.project_name, "prepro/dtm_13.val"),
)

fig, ax = plt.subplots(1)
img = ax.imshow(dem_mat)
plt.colorbar(img)


# simu.show_input(prop="dem")


# simu.update_prepo_inputs(
#     DEM=dem_mat,
#     # N=np.shape(dem_mat)[1],
#     # M=np.shape(dem_mat)[0],
# )

# fig = plt.figure()
# ax = plt.axes(projection="3d")
# simu.show_input(prop="dem", ax=ax)
# simu.create_mesh_vtk(verbose=True)

delta_x = 0.2
delta_y = 0.2
ncellx = int(np.round(14.2/delta_x))
ncelly = int(np.round(2/delta_y))
DEM = np.ones([ncelly,ncellx])*1e-1 # raster of 20x20
DEM[-1,-1]=0.1-1e-3
maxdepth = 2-1e-1
zb = np.linspace(0, maxdepth, 15)   
# zb = np.linspace(0, maxdepth, 12)
# zb = [np.round(zi,1) for zi in zb]


nstr = len(zb)
zr = list((np.ones(len(zb))) / (nstr))
len(zr)
    
# maxdepth = 10

# # the fraction of total grid height that each layer is to occupy
# # log z depth
# # -------------------------------------------------------------
# zb = np.geomspace(1e-1, maxdepth, num=20)
# nstr = len(zb)
# zr = [abs(zb[0] / maxdepth)]
# zr.extend(list(abs(np.diff(zb) / maxdepth)))


simu.update_prepo_inputs(
    DEM=dem_mat,
    xllcorner=1e4,
    yllcorner=4e3,
    # nstr=20,
    zratio=zr,
    base=max(zb),
)

simu.update_parm(TRAFLAG=0)
simu.create_mesh_vtk(verbose=True)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)