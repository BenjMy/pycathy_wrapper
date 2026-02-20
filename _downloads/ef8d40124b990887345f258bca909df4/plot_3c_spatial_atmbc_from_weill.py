#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spatially heterogeneous atmbc (Weil et al)
==========================================

Weill, S., et al. « Coupling Water Flow and Solute Transport into a Physically-Based Surface–Subsurface Hydrological Model ».
Advances in Water Resources, vol. 34, no 1, janvier 2011, p. 128‑36. DOI.org (Crossref),
https://doi.org/10.1016/j.advwatres.2010.10.001.


This example shows how to use pyCATHY object to create spatially and temporally variable atmbc conditions

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
                         prj_name="atmbc_spatially_from_weill"
                         )

simu.create_mesh_vtk(verbose=True)

#%% spatially variable atmospheric boundary condition inputs

grid3d = simu.read_outputs('grid3d')

# np.shape(simu.DEM)

# DEM, dem_header = simu.read_inputs('dem')
t_atmbc = [0,86400]
v_atmbc = np.zeros(int(grid3d['nnod']))
v_atmbc[0:int(len(np.zeros(int(grid3d['nnod'])))/2)] = 1e-7

v_atmbc_mat = np.reshape(v_atmbc,[np.shape(simu.DEM)[0]+1,
                                  np.shape(simu.DEM)[1]+1
                                  ])
fig, ax = plt.subplots()
ax.imshow(v_atmbc_mat)

# np.shape([v_atmbc]*len(t_atmbc))

simu.update_atmbc(
                    HSPATM=0,
                    IETO=0,
                    time=t_atmbc,
                    netValue=[v_atmbc]*len(t_atmbc)
                  )

#%%
simu.run_processor(IPRT1=2,
                    DTMIN=1e-2,
                    DTMAX=1e2,
                    DELTAT=5,
                   TRAFLAG=0,
                   verbose=False
                   )

# cplt.show_spatial_atmbc()

#%%

cplt.show_vtk(
    unit="pressure",
    timeStep=1,
    notebook=False,
    path=simu.workdir + simu.project_name + "/vtk/",
    savefig=True,
)

#%%  3d visualiation of the water saturation for the  all time steps

cplt.show_vtk_TL(
                unit="pressure",
                notebook=False,
                path=simu.workdir + simu.project_name + "/vtk/",
                show=False,
                x_units='days',
                clim = [0.55,0.70],
                savefig=True,
            )
