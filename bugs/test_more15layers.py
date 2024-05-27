#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Changing mesh properties (resolution, coordinates, ...)
=======================================================

Weill, S., et al. « Coupling Water Flow and Solute Transport into a Physically-Based Surface–Subsurface Hydrological Model ». 
Advances in Water Resources, vol. 34, no 1, janvier 2011, p. 128‑36. DOI.org (Crossref), 
https://doi.org/10.1016/j.advwatres.2010.10.001.

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
simu = cathy_tools.CATHY(dirName=path2prj, 
			prj_name="meshing_from_weill"
				)

rootpath = os.path.join(simu.workdir + simu.project_name)

#%%



#%% Fetch and show initial DEM
# the dimension of the mesh is squared (20,20)


dem_mat, str_hd_dem = in_CT.read_dem(
    os.path.join(simu.workdir, simu.project_name, "prepro/dem"),
    os.path.join(simu.workdir, simu.project_name, "prepro/dtm_13.val"),
)

#%% Change number of layers and maximum depth

maxdepth = 1

maproot, _ = simu.read_inputs('root_map')
dem_parameters = in_CT.read_dem_parameters(
                                                os.path.join(simu.workdir, 
                                                             simu.project_name,
                                                             "input",
                                                             "dem_parameters")
                                                )
nstr = dem_parameters['nstr']

simu.update_prepo_inputs(
    DEM=dem_mat,
    nstr=nstr,
)
SPP, FP = simu.read_inputs('soil', MAXVEG=len(np.unique(maproot)))
extended_layers = 20
extended_index = pd.MultiIndex.from_product([[1], range(extended_layers)], names=['zone', 'layer'])
soil_FP_newlayers = SPP.reindex(extended_index).fillna(method='ffill')

# linear z depth
# -------------------------------------------------------------
zb = np.linspace(0, maxdepth, extended_layers)
nstr = len(zb)
zr = list((np.ones(len(zb))) / (nstr))


simu.update_prepo_inputs(
    DEM=dem_mat,
    xllcorner=1e4,
    yllcorner=4e3,
    nstr=nstr,
    zratio=zr,
    base=max(zb),
)

simu.update_soil(SPP=soil_FP_newlayers,
                 FP=FP)

# fig = plt.figure()
# ax = plt.axes(projection="3d")
# simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=True)
#%% Plot mesh
meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"
import pyvista as pv

mesh2plot = pv.read(meshfile)
mesh2plot.plot(show_edges=True, show_axes=True, show_bounds=True)