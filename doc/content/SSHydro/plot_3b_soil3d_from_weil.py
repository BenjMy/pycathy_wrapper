#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Soil 3d from a Digital Elevation Model (DEM)
============================================

This example shows how to use pyCATHY object to build a 3d soil properties from a DEM and run the hydrological model.

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
simu = cathy_tools.CATHY(dirName=path2prj, prj_name="soil3d_from_weil", clear_src=True)

rootpath = os.path.join(simu.workdir + simu.project_name)
#%% Fetch and show initial DEM
# the dimension of the mesh is squared (20,20)


dem_mat, str_hd_dem = in_CT.read_dem(
    os.path.join(simu.workdir, simu.project_name, "prepro/dem"),
    os.path.join(simu.workdir, simu.project_name, "prepro/dtm_13.val"),
)

simu.show_input(prop="dem")

maxdepth = 10
# # linear z depth
# # -------------------------------------------------------------
zb = np.linspace(0, maxdepth, 10)
nstr = len(zb)
zr = list((np.ones(len(zb))) / (nstr))

sum(zr)

# zb = np.geomspace(1e-1, maxdepth, num=15)
# nstr=len(zb)
# zr = [abs(zb[0]/maxdepth)]
# zr.extend(list(abs(np.diff(zb)/maxdepth)))

# np.shape(dem_mat)
simu.update_prepo_inputs(
    DEM=dem_mat,
    nstr=nstr,
    zratio=zr,
    base=max(zb),
)
fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)

# simu.update_soil()
simu.create_mesh_vtk(verbose=True)

#%% Define and map layers
simu.update_zone()

layers = {1: [0, 2], 2: [2, 6], 3: [6, 10]}

zone3d_flag = mt.map_layers_2_DEM(layers, simu.DEM, simu.zone, simu.dem_parameters)

#%% Define a dictionnary of Soil Physical Properties

SPP_map = {
    "PERMX": [0.000188] * 3,
    "PERMY": [0.000188] * 3,
    "PERMZ": [0.000188] * 3,
    "ELSTOR": [1e-05] * 3,
    "POROS": [0.55, 0.65, 0.5],
    "VGNCELL": [1.46, 1.46, 1.46],
    "VGRMCCELL": [0.15, 0.15, 0.15],
    "VGPSATCELL": [0.03125] * 3,
}


#%% Update the soil file

simu.update_soil(
    SPP_map=SPP_map,
    zone3d=zone3d_flag,
)


#%% Run  hydrological model

simu.run_processor(IPRT1=2, verbose=True)
