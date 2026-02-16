#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshing from a Digital Elevation Model (DEM)
============================================

GRwater project 

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

import rioxarray
import pyvista as pv



#%% Init CATHY model
# ------------------------
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name="meshing_from_subcachment"
                         )

rootpath = os.path.join(simu.workdir + simu.project_name)


# Path to the directory containing the .adf file (not the file itself)
adf_folder = "../data/dtmplot1/"

# Open the raster (typically named 'hdr.adf', but you only need the folder)
raster_DEM = rioxarray.open_rasterio(adf_folder, masked=True).isel(band=0)

# Create a mask of valid (non-NaN) data
valid_mask = ~np.isnan(raster_DEM)

# Apply the mask to get valid coordinates
valid_x = raster_DEM['x'].where(valid_mask.any(dim='y'), drop=True)
valid_y = raster_DEM['y'].where(valid_mask.any(dim='x'), drop=True)

# Get min and max valid coordinates
min_lon, max_lon = float(valid_x.min()), float(valid_x.max())
min_lat, max_lat = float(valid_y.min()), float(valid_y.max())

raster_DEM_masked = raster_DEM.where(
    (raster_DEM['x'] >= min_lon) & (raster_DEM['x'] <= max_lon), drop=True
).where(
    (raster_DEM['y'] >= min_lat) & (raster_DEM['y'] <= max_lat), drop=True
)
       
raster_DEM_masked = np.where(np.isnan(raster_DEM_masked), -9999, raster_DEM_masked)
np.shape(raster_DEM_masked)
np.shape(raster_DEM)

#%% Fetch and show initial DEM

fig, ax = plt.subplots(1)
img = ax.imshow(raster_DEM_masked)
plt.colorbar(img)

simu.show_input(prop="dem")

simu.update_prepo_inputs(
    DEM=raster_DEM_masked,
    delta_x=5,
    delta_y=5,
    ivert=1,
)

fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=True)

#%% Plot mesh
# meshfile = rootpath + "/vtk/" + simu.project_name + ".vtk"

# mesh2plot = pv.read(meshfile)

# pl = pv.Plotter(off_screen=True)

# pl.add_mesh(mesh2plot)
# pl.show_bounds()

# pl.show()
# mesh2plot.plot(show_edges=True, 
#                show_axes=True, 
#                show_bounds=True
#                )

#%%

simu.run_preprocessor(verbose=False)
# simu.run_processor(IPRT1=3,verbose=True)

# # simu.read_inputs('atmbc')
simu.update_parm(TIMPRTi=[1800,7200])

grid3d = simu.read_outputs('grid3d')

t_atmbc = [0,86400,86400*4]
v_atmbc_t_rain = 1e-7
netValue = np.zeros(len(t_atmbc))
netValue[0] = v_atmbc_t_rain
# fig, ax = plt.subplots()
# ax.imshow(v_atmbc)


#% Create an empty dataframe of SPP and set default SPP properties 
df_SPP_map = simu.init_soil_SPP_map_df(nzones=1,nstr=15)
SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)


SPP_map['PERMX'] = 6.88e-4
SPP_map['PERMY'] = 6.88e-4
SPP_map['PERMZ'] = 6.88e-4

#% Update soil file
simu.update_soil(SPP_map=SPP_map)


simu.update_ic(INDP=0,
               IPOND=0,
               pressure_head_ini=-15
                  )


simu.update_atmbc(
                    HSPATM=1,
                    IETO=1,
                    time=t_atmbc,
                    netValue=netValue
                  )

simu.update_parm(
                        IPRT=4,
                        VTKF=2, # dont write vtk files
                        )

#%%

int(grid3d['nnod'])


# simu.run_processor(IPRT1=2,
#                     DTMIN=1e-2,
#                     DTMAX=1e3,
#                     DELTAT=1e2,
#                     TRAFLAG=0,
#                     verbose=True
#                     )

#%%


# # pl = pv.Plotter(notebook=True)
# pl = pv.Plotter(off_screen=True)



# cplt.show_vtk(unit="pressure",
#               timeStep=1,
#               path=simu.workdir + "/meshing_from_subcachment/vtk/",
#               ax=pl,
#               )
# # pl.show()
# image = pl.screenshot('test.png')


