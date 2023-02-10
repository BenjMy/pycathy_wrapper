#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshing from a Digital Elevation Model (DEM)
============================================

This example shows how to use pyCATHY object to mesh from a DEM and run the hydrological model.

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

#%% Init CATHY model
# ------------------------
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, prj_name="meshing_from_weil", clear_src=True)

rootpath = os.path.join(simu.workdir + simu.project_name)
#%% Fetch and show initial DEM
# the dimension of the mesh is squared (20,20)


dem_mat, str_hd_dem = in_CT.read_dem(
    os.path.join(simu.workdir, simu.project_name, "prepro/dem"),
    os.path.join(simu.workdir, simu.project_name, "prepro/dtm_13.val"),
)

simu.show_input(prop="dem")

print(dem_mat)

simu.update_prepo_inputs(
    DEM=dem_mat,
    # N=np.shape(dem_mat)[1],
    # M=np.shape(dem_mat)[0],
)

fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=True)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)


#%% Crop the mesh
# the new dimension of the mesh is rectangle (10,20)

dem_crop = dem_mat[0:10, :]
print("DEM shape is {}".format(np.shape(dem_crop)))

simu.update_prepo_inputs(
    DEM=dem_crop,
)

simu.update_zone()
simu.update_veg_map()

fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=True)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)


#%% Change mesh resolution


#%% Change coordinates offset

simu.update_prepo_inputs(
    DEM=dem_crop,
    xllcorner=1e4,
    yllcorner=4e3,
)

simu.update_zone()
simu.update_veg_map()

fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=False)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)

#%% Flip y axis

dem_crop_flipy = np.flipud(dem_crop)

simu.update_prepo_inputs(
    DEM=dem_crop_flipy,
    xllcorner=1e4,
    yllcorner=4e3,
)

simu.update_zone()
simu.update_veg_map()
fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=False)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)


#%% Change number of layers and maximum depth

dem_crop_3layers = np.flipud(dem_crop)
maxdepth = 10

# linear z depth
# -------------------------------------------------------------
zb = np.linspace(0, maxdepth, 3)
nstr = len(zb) - 1
zr = list((np.ones(len(zb))) / (nstr))


simu.update_prepo_inputs(
    DEM=dem_crop,
    xllcorner=1e4,
    yllcorner=4e3,
    nstr=nstr,
    zratio=zr,
    base=max(zb),
)
fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=False)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)


#%% Change number of layers (log)

# the fraction of total grid height that each layer is to occupy
# log z depth
# -------------------------------------------------------------
zb = np.geomspace(1e-1, maxdepth, num=15)
nstr = len(zb)
zr = [abs(zb[0] / maxdepth)]
zr.extend(list(abs(np.diff(zb) / maxdepth)))


simu.update_prepo_inputs(
    DEM=dem_crop,
    xllcorner=1e4,
    yllcorner=4e3,
    nstr=nstr,
    zratio=zr,
    base=max(zb),
)
simu.create_mesh_vtk(verbose=False)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)


#%% Run  hydrological model

simu.run_processor(IPRT1=2, verbose=True)
