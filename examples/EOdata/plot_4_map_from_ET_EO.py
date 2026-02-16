#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapping atmbc from Earth Observation xarray dataset
===================================================

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



#%% Init CATHY model
# ------------------------
path2prj = "../EOdata/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name="mapping_atmbc_from_EO_dataset"
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


#%%
simu.create_mesh_vtk(verbose=True)
#%%


import xarray as xr

simu.mesh_pv_attributes

# # --- Convert to xarray.Dataset without ETa ---
ds_mesh = xr.Dataset(
    coords={
        "node": np.arange(20336),
        "x": ("node", simu.mesh_pv_attributes.points[:, 0]),
        "y": ("node", simu.mesh_pv_attributes.points[:, 1]),
        "z": ("node", simu.mesh_pv_attributes.points[:, 2])
    },
    attrs={
        "N_cells": 105660,
        # "nodes_per_element": nodes_per_element
    }
)


ds_mesh["mask"] = (("y", "x"), raster_DEM_masked)

# ---- 2D structured grid from hapin ----
fig, ax = plt.subplots(figsize=(6, 5))
ax.set_title("Structured 2D grid (hapin)")
ax.set_aspect("equal")

# Plot grid lines
for xval in ds_mesh.x.values:
    ax.plot([xval, xval], [ds_mesh.y.values.min(), ds_mesh.y.values.max()], color="lightgrey", lw=0.5)
for yval in ds_mesh.y.values:
    ax.plot([ds_mesh.x.values.min(), ds_mesh.x.values.max()], [yval, yval], color="lightgrey", lw=0.5)

ax.set_xlabel("X")
ax.set_ylabel("Y")

plt.show()

ds_mesh["mask"].plot.imshow()

X_nodes = ds_mesh.x.values
Y_nodes = ds_mesh.y.values
Z_nodes = ds_mesh.z.values


# 2D mask from DEM
# mask2d = raster_DEM_masked != -9999  # shape (M, N)
# ds_mesh["mask2d"] = (("y", "x"), mask2d)

# Extract parameters from hapin
dx = simu.hapin["delta_x"]
dy = simu.hapin["delta_y"]
x0 = simu.hapin["xllcorner"]
y0 = simu.hapin["yllcorner"]


# Convert to indices
ix = ((X_nodes - x0) / dx).astype(int)
iy = ((Y_nodes - y0) / dy).astype(int)

# Clip indices inside raster bounds
ix = np.clip(ix, 0, simu.hapin["N"] - 1)
iy = np.clip(iy, 0, simu.hapin["M"] - 1)


# Build boolean mask for nodes
mask_array = ds_mesh["mask"].values  # convert to NumPy array
bool_mask_nodes = mask_array[iy, ix]      # shape: (node,)

# Add to xarray
ds_mesh["mask_node"] = (("node",), bool_mask_nodes)

#%%

import numpy as np
import xarray as xr

# Use simu.hapin
hapin = simu.hapin

# Grid info
N = hapin["N"]
M = hapin["M"]
dx = hapin["delta_x"]
dy = hapin["delta_y"]
x0 = hapin["xllcorner"]
y0 = hapin["yllcorner"]

# Cell centers
x = x0 + (np.arange(N) + 0.5) * dx
y = y0 + (np.arange(M) + 0.5) * dy

# # Create meshgrid for smooth patterns
X, Y = np.meshgrid(x, y)

# Example ETa / ETp arrays (replace with your actual values)
# Example ETa / ETp 2D arrays
ETa_2D = 3 + 2 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())
ETp_2D = 4 + 1.5 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())


# Create daily time axis for a full month (30 days)
time = pd.date_range("2026-02-01", periods=30, freq="D")

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

#%% Map ETa/ETp to nodes, set NaN where masked

ETp_nodes = np.where(ds_mesh['mask_node'], ETp_3D[:, iy, ix], np.nan)
ds_mesh["ETp"] =  (("time","node"), ETp_nodes)
ds_mesh["ETp_surfacenodes"] = ds_mesh["ETp"].isel(node=slice(0,int(simu.grid3d['nnod'])))


#%%

# a
import matplotlib.pyplot as plt

# Select time step 0
t = 0

# --- Raster ETp from ds_et ---
ETp_raster = ds_et["ETp"].isel(time=t).values

# --- Node ETp from ds_mesh ---
ETp_nodes = ds_mesh["ETp_surfacenodes"].isel(time=t).values
x_nodes = ds_mesh["x"].values[:len(ETp_nodes)]
y_nodes = ds_mesh["y"].values[:len(ETp_nodes)]

# --- Plot side by side ---
fig, axes = plt.subplots(1, 2, figsize=(16,6))

# Raster
im0 = axes[0].imshow(ETp_raster, origin='lower',
                     extent=[x.min(), x.max(), y.min(), y.max()],
                     aspect='auto', cmap='viridis')
axes[0].set_title("ETp (raster ds_et)")
axes[0].set_xlabel("x [m]")
axes[0].set_ylabel("y [m]")
fig.colorbar(im0, ax=axes[0], label="ETp [mm/day]")

# Nodes (unstructured mesh)
sc = axes[1].scatter(x_nodes, y_nodes, c=ETp_nodes, cmap='viridis', s=20)
axes[1].set_title("ETp (mesh nodes ds_mesh)")
axes[1].set_xlabel("X [m]")
axes[1].set_ylabel("Y [m]")
axes[1].axis('equal')
fig.colorbar(sc, ax=axes[1], label="ETp [mm/day]")

plt.tight_layout()
plt.show()


#%% Update atmbc

# Convert xarray DataArray to NumPy array
ETp_nodes = ds_mesh["ETp_surfacenodes"].values  # shape: (time, node)

# Identify columns (nodes) that are all NaN
valid_nodes = ~np.all(np.isnan(ETp_nodes), axis=0)  # True for nodes with any valid value

# Keep only valid nodes
ETp_nodes_clean = ETp_nodes[:, valid_nodes]

# Reference time: first time step
t0 = ds_et.time[0].values

# Convert all times to seconds since t0
time_sec = (ds_et.time.values - t0) / np.timedelta64(1, 's')  # seconds

# Add as a new coordinate if you want
ds_et = ds_et.assign_coords(time_sec=("time", time_sec))


simu.update_atmbc(HSPATM=0,
                  IETO=1,
                  time = ds_et.time_sec,
                  netValue = -ETp_nodes_clean*0.001 / 86400
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
