#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapping atmbc from coarse EO
============================

- 1. create a CATHY object 
- 2. create a CATHY mesh from a catchment DEM from a file .adf containing no values -9999
- 3. wrap the mesh properties (from hapin file) to an xarray dataset

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
import xarray as xr



#%% Init CATHY model
# ------------------------
path2prj = "../EOdata/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name="mapping_atmbc_from_coarse_EO_dataset"
                         )

rootpath = os.path.join(simu.workdir + simu.project_name)


# Path to the directory containing the .adf file (not the file itself)
adf_folder = "../../data/dtmplot1/"

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

simu.update_prepo_inputs(
    DEM=raster_DEM_masked,
    delta_x=5,
    delta_y=5,
    ivert=1,
)



#%%
simu.create_mesh_vtk(verbose=True)
#%%
ds_mesh = mt.build_mesh_dataset(simu, 
                                raster_DEM_masked=raster_DEM_masked, 
                                plot_grid=True
                                )

#%% Define EO ETp (synthetic)


def create_et_dataset(simu, ETa_2D=None, ETp_2D=None, start_date="2026-02-01", n_days=30):
    """
    Create a 3D xarray.Dataset of ETa and ETp from a simu.hapin grid.

    Parameters
    ----------
    simu : object
        Simulation object with hapin attributes (delta_x, delta_y, N, M, xllcorner, yllcorner).
    ETa_2D : 2D numpy array, optional
        Daily ETa pattern (shape M x N). If None, a synthetic pattern will be generated.
    ETp_2D : 2D numpy array, optional
        Daily ETp pattern (shape M x N). If None, a synthetic pattern will be generated.
    start_date : str, optional
        Start date for the time axis. Default "2026-02-01".
    n_days : int, optional
        Number of days to generate. Default 30.

    Returns
    -------
    ds_et : xarray.Dataset
        Dataset with variables "ETa" and "ETp" (time, y, x).
    """

    # Use hapin grid
    hapin = simu.hapin
    N = hapin["N"]
    M = hapin["M"]
    dx = hapin["delta_x"]
    dy = hapin["delta_y"]
    x0 = hapin["xllcorner"]
    y0 = hapin["yllcorner"]

    # Cell centers
    x = x0 + (np.arange(N) + 0.5) * dx
    y = y0 + (np.arange(M) + 0.5) * dy

    # Create meshgrid
    X, Y = np.meshgrid(x, y)

    # Synthetic patterns if none provided
    if ETa_2D is None:
        ETa_2D = 3 + 2 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())
    if ETp_2D is None:
        ETp_2D = 4 + 1.5 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())

    # Time axis
    time = pd.date_range(start_date, periods=n_days, freq="D")

    # Expand 2D arrays to 3D (time, y, x)
    ETa_3D = np.tile(ETa_2D[None, :, :], (len(time), 1, 1))
    ETp_3D = np.tile(ETp_2D[None, :, :], (len(time), 1, 1))

    # Create xarray Dataset
    ds_et = xr.Dataset(
        data_vars={
            "ETa": (("time", "y", "x"), ETa_3D),
            "ETp": (("time", "y", "x"), ETp_3D),
        },
        coords={
            "time": time,
            "x": x,
            "y": y,
        },
        attrs={
            "resolution_x": dx,
            "resolution_y": dy,
            "xllcorner": x0,
            "yllcorner": y0,
        }
    )

    return ds_et

ds_et = create_et_dataset(simu, start_date="2026-02-01", n_days=30)
ds_et.rio.resolution()
# ds_et.rio.crs

#%% Coarsen EO (so its resolution is lower than CATHY mesh resolution)

# Define coarsening factor
# Example: make resolution 2x coarser in x and y
factor_x = 5  # combine 2 grid cells along x
factor_y = 5 # combine 2 grid cells along y

# Coarsen and take mean
ds_et_coarse = ds_et.coarsen(x=factor_x, y=factor_y, boundary='trim').mean()

# Update attributes to reflect new resolution
ds_et_coarse.attrs['resolution_x'] = ds_et.attrs['resolution_x'] * factor_x
ds_et_coarse.attrs['resolution_y'] = ds_et.attrs['resolution_y'] * factor_y

# print(ds_et_coarse)

print(ds_et_coarse.rio.resolution())

#%%

import xarray as xr
import numpy as np

# ds_mapped = mt.map_grid_to_mesh(ds_et_coarse, ds_mesh, variables=["ETp", "ETa"])
ds_mapped = mt.map_grid_to_mesh(ds_et, ds_mesh, variables=["ETp", "ETa"])
# ds_mapped = map_grid_to_mesh(ds_era5, ds_mesh, method="linear")
# ds_mapped = map_grid_to_mesh(ds_modis, ds_mesh, variables=["NDVI"],
#                              x_dim="lon", y_dim="lat")

ds_mesh_surfacenodes = ds_mapped.isel(node=slice(0,int(simu.grid3d['nnod'])))

t0 = ds_mesh_surfacenodes.isel(time=0)['ETp']

fig,ax= plt.subplots(1,2)
ds_et_coarse['ETp'].isel(time=0).plot.imshow(ax=ax[0])
sc = ax[1].scatter(t0.x, t0.y, c=t0.values, cmap='viridis', s=5)
plt.colorbar(sc, ax=ax[1], label='ETp')
ax[1].set_title(str(t0.time.values)[:10])
ax[1].set_xlabel('x'); ax[1].set_ylabel('y')
plt.show()


#%% Update atmbc

# Convert xarray DataArray to NumPy array
ETp_nodes = ds_mesh_surfacenodes["ETp"].values  # shape: (time, node)


len(ETp_nodes)
# Reference time: first time step
t0 = ds_et.time[0].values

# Convert all times to seconds since t0
time_sec = (ds_et.time.values - t0) / np.timedelta64(1, 's')  # seconds

# Add as a new coordinate if you want
ds_et = ds_et.assign_coords(time_sec=("time", time_sec))


simu.update_atmbc(HSPATM=0,
                  IETO=1,
                  time = ds_et.time_sec,
                  netValue = -ETp_nodes*0.001 / 86400
                  )


#%% Run processor


simu.update_ic(INDP=0,pressure_head_ini=-10)

simu.update_parm(TIMPRTi=ds_et.time_sec.values,
                 VTKF=2
                 )

simu.run_processor(IPRT1=2,
                    DTMIN=1e-2,
                    DTMAX=1e2,
                    DELTAT=5,
                    TRAFLAG=0,
                    verbose=True)
