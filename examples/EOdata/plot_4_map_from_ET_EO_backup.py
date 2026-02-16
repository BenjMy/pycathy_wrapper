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

#%%
import numpy as np

# Map ETa/ETp to nodes, set NaN where masked
# ETa_nodes = np.where(mask_nodes, ETa_3D[:, iy, ix], np.nan)
ETp_nodes = np.where(ds_mesh['mask_node'], ETp_3D[:, iy, ix], np.nan)

# Add to 3D mesh dataset
# ds_mesh["ETa"] = (("time","node"), ETa_nodes)
ds_mesh["ETp"] =  (("time","node"), ETp_nodes)
# ds_mesh["mask3d"] =  (("time","node"), mask_nodes_time)

# ds_mesh["ETp"]
# ds_mesh["ETp"].isel(node=slice(0,int(simu.grid3d['nnod'])))

ds_mesh["ETp_surfacenodes"] = ds_mesh["ETp"].isel(node=slice(0,int(simu.grid3d['nnod'])))

#%%
import numpy as np

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
                  netValue = ETp_nodes_clean
                  )


#%%

simu.run_processor(IPRT1=2,verbose=True)

#%%
# simu.run_preprocessor(verbose=False)
# # simu.run_processor(IPRT1=3,verbose=True)

# # # simu.read_inputs('atmbc')
# simu.update_parm(TIMPRTi=[1800,7200])

grid3d = simu.read_outputs('grid3d')

# t_atmbc = [0,86400,86400*4]
# v_atmbc_t_rain = 1e-7
# netValue = np.zeros(len(t_atmbc))
# netValue[0] = v_atmbc_t_rain
# # fig, ax = plt.subplots()
# # ax.imshow(v_atmbc)


# #% Create an empty dataframe of SPP and set default SPP properties 
# df_SPP_map = simu.init_soil_SPP_map_df(nzones=1,nstr=15)
# SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)


# SPP_map['PERMX'] = 6.88e-4
# SPP_map['PERMY'] = 6.88e-4
# SPP_map['PERMZ'] = 6.88e-4

# #% Update soil file
# simu.update_soil(SPP_map=SPP_map)


# simu.update_ic(INDP=0,
#                IPOND=0,
#                pressure_head_ini=-15
#                   )


# simu.update_atmbc(
#                     HSPATM=1,
#                     IETO=1,
#                     time=t_atmbc,
#                     netValue=netValue
#                   )

# simu.update_parm(
#                         IPRT=4,
#                         VTKF=2, # dont write vtk files
#                         )


#%%

# # ---- 3D mesh nodes scatterplot ----
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection="3d")
# ax.set_title("3D mesh nodes")


# ax.scatter(x, y, z, s=2, c=z, cmap="viridis")

# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")

# plt.show()



#%%
# simu.grid3d

# simu.hapin
# simu.hapin
# Out[8]: 
# {'delta_x': 5,
#  'delta_y': 5,
#  'N': 61,
#  'M': 35,
#  'N_celle': 400,
#  'xllcorner': 0.0,
#  'yllcorner': 0.0,

simu.grid3d
# Out[2]: 
# {'nnod': np.float64(1271.0),
#  'nnod3': np.float64(20336.0),
#  'nel': np.float64(105660.0),
#  'mesh3d_nodes': array([[245.  , 175.  , 622.04],
#         [250.  , 175.  , 623.28],
#         [255.  , 175.  , 623.9 ],
#         ...,
#         [ 85.  ,   0.  , 564.96],
#         [ 90.  ,   0.  , 564.96],
#         [ 95.  ,   0.  , 564.96]]),
#  'mesh_tetra': array([[1.0000e+00, 5.0000e+00, 6.0000e+00, 1.2720e+03, 1.0000e+00],
#         [1.2720e+03, 1.2760e+03, 1.2770e+03, 6.0000e+00, 1.0000e+00],
#         [5.0000e+00, 6.0000e+00, 1.2760e+03, 1.2720e+03, 1.0000e+00],
#         ...,
#         [1.9054e+04, 1.9055e+04, 1.9065e+04, 2.0325e+04, 1.0000e+00],
#         [2.0325e+04, 2.0326e+04, 2.0336e+04, 1.9065e+04, 1.0000e+00],
#         [1.9055e+04, 1.9065e+04, 2.0326e+04, 2.0325e+04, 1.0000e+00]])}


# #%%


#%%
# import xarray as xr
# # Build coordinate vectors from hapin
# x = simu.hapin['xllcorner'] + np.arange(simu.hapin['N']) * simu.hapin['delta_x']
# y = simu.hapin['yllcorner'] + np.arange(simu.hapin['M']) * simu.hapin['delta_y']
# grid3d = simu.grid3d

# # Wrap into xarray
# ds_hapin = xr.Dataset(
#     data_vars={
#         "mesh3d_nodes": (("node", "coord"), grid3d["mesh3d_nodes"]),
#         "mesh_tetra": (("element", "nodes_per_element"), grid3d["mesh_tetra"]),
#     },
#     coords={
#         "x": ("x", x),
#         "y": ("y", y),
#         "node": np.arange(grid3d["mesh3d_nodes"].shape[0]),
#         "coord": ["X", "Y", "Z"],
#         "element": np.arange(grid3d["mesh_tetra"].shape[0]),
#         "nodes_per_element": np.arange(grid3d["mesh_tetra"].shape[1]),
#     },
#     attrs={
#         "delta_x": simu.hapin["delta_x"],
#         "delta_y": simu.hapin["delta_y"],
#         "N": simu.hapin["N"],
#         "M": simu.hapin["M"],
#         "N_celle": simu.hapin["N_celle"],
#         "nnod": grid3d["nnod"],
#         "nnod3": grid3d["nnod3"],
#         "nel": grid3d["nel"],
#     }
# )

# print(ds_hapin)
# ds_hapin["mask"] = (("y", "x"), raster_DEM_masked)



#%%


#%%

# import matplotlib.pyplot as plt

# # ---- 2D structured grid from hapin ----
# fig, ax = plt.subplots(figsize=(6, 5))
# ax.set_title("Structured 2D grid (hapin)")
# ax.set_aspect("equal")

# # Plot grid lines
# for xval in ds_hapin.x.values:
#     ax.plot([xval, xval], [ds_hapin.y.values.min(), ds_hapin.y.values.max()], color="lightgrey", lw=0.5)
# for yval in ds_hapin.y.values:
#     ax.plot([ds_hapin.x.values.min(), ds_hapin.x.values.max()], [yval, yval], color="lightgrey", lw=0.5)

# ax.set_xlabel("X")
# ax.set_ylabel("Y")

# plt.show()


# # ---- 3D mesh nodes scatterplot ----
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection="3d")
# ax.set_title("3D mesh nodes")

# X = ds_hapin.mesh3d_nodes.sel(coord="X").values
# Y = ds_hapin.mesh3d_nodes.sel(coord="Y").values
# Z = ds_hapin.mesh3d_nodes.sel(coord="Z").values

# ax.scatter(X, Y, Z, s=2, c=Z, cmap="viridis")

# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")

# plt.show()



#%%

# ds_hapin.mask.plot.imshow()

# # 2D mask from DEM
# mask2d = raster_DEM_masked != -9999  # shape (M, N)
# ds_hapin["mask2d"] = (("y", "x"), mask2d)

# # Extract parameters from hapin
# dx = simu.hapin["delta_x"]
# dy = simu.hapin["delta_y"]
# x0 = simu.hapin["xllcorner"]
# y0 = simu.hapin["yllcorner"]

# # Node coordinates
# X = ds_hapin.mesh3d_nodes.sel(coord="X").values
# Y = ds_hapin.mesh3d_nodes.sel(coord="Y").values

# # Convert to indices
# ix = ((X - x0) / dx).astype(int)
# iy = ((Y - y0) / dy).astype(int)

# # Clip indices inside raster bounds
# ix = np.clip(ix, 0, simu.hapin["N"] - 1)
# iy = np.clip(iy, 0, simu.hapin["M"] - 1)

# # Build mask for nodes
# mask3d = mask2d[iy, ix]   # boolean mask for each node
# ds_hapin["mask3d"] = (("node",), mask3d)

# #%%

# # 2D mask
# ds_hapin.mask2d.plot.imshow()

# # 3D nodes, masked
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection="3d")
# ax.set_title("3D mesh nodes (masked)")

# X = ds_hapin.mesh3d_nodes.sel(coord="X").values
# Y = ds_hapin.mesh3d_nodes.sel(coord="Y").values
# Z = ds_hapin.mesh3d_nodes.sel(coord="Z").values

# valid = ds_hapin.mask3d.values
# ax.scatter(X[valid], Y[valid], Z[valid], s=2, c=Z[valid], cmap="viridis")

# plt.show()


# ds_et.ETa.plot.imshow()

#%%
ds_et.ETp.plot.imshow()


# ds_hapin
# # Use the time axis from ds_et
# ds_hapin = ds_hapin.assign_coords(time=ds_et.time)

# print(ds_hapin)


#%% Coarsening since EO values are often more coarse than hydrological model cell resolution

# import numpy as np
# import xarray as xr

# hapin = simu.hapin

# # Original grid info
# N = hapin["N"]
# M = hapin["M"]
# dx = hapin["delta_x"]
# dy = hapin["delta_y"]
# x0 = hapin["xllcorner"]
# y0 = hapin["yllcorner"]

# # Coarsening factor
# factor = 3

# # Coarse grid sizes
# N_coarse = N // factor
# M_coarse = M // factor
# dx_coarse = dx * factor
# dy_coarse = dy * factor

# # Coarse cell centers
# x_coarse = x0 + (np.arange(N_coarse) + 0.5) * dx_coarse
# y_coarse = y0 + (np.arange(M_coarse) + 0.5) * dy_coarse

# # Meshgrid
# Xc, Yc = np.meshgrid(x_coarse, y_coarse)

# # Smooth ETa / ETp patterns on coarse grid
# ETa_coarse = 3 + 2 * np.sin(np.pi * Xc / Xc.max()) * np.cos(np.pi * Yc / Yc.max())
# ETp_coarse = 4 + 1.5 * np.sin(np.pi * Xc / Xc.max()) * np.cos(np.pi * Yc / Yc.max())

# # Create xarray Dataset
# ds_et_coarse = xr.Dataset(
#     data_vars={
#         "ETa": (("y", "x"), ETa_coarse),
#         "ETp": (("y", "x"), ETp_coarse),
#     },
#     coords={
#         "x": ("x", x_coarse),
#         "y": ("y", y_coarse),
#     },
#     attrs={
#         "resolution_x": dx_coarse,
#         "resolution_y": dy_coarse,
#         "xllcorner": x0,
#         "yllcorner": y0,
#     }
# )

# print(ds_et_coarse)


#%%

# import numpy as np

# # 3D node coordinates
# X_nodes = ds_hapin.mesh3d_nodes.sel(coord="X").values
# Y_nodes = ds_hapin.mesh3d_nodes.sel(coord="Y").values

# # ET 2D grid and coordinates
# ETa_2d = ds_et.ETa.values
# ETp_2d = ds_et.ETp.values
# x_grid = ds_et.x.values
# y_grid = ds_et.y.values

# # Convert node XY to raster indices
# dx = hapin["delta_x"]
# dy = hapin["delta_y"]
# x0 = hapin["xllcorner"]
# y0 = hapin["yllcorner"]

# ix = ((X_nodes - x0) / dx).astype(int)
# iy = ((Y_nodes - y0) / dy).astype(int)

# # Clip indices to valid range
# ix = np.clip(ix, 0, hapin["N"] - 1)
# iy = np.clip(iy, 0, hapin["M"] - 1)

# # Apply 2D mask
# mask2d = ds_hapin.mask2d.values
# mask_nodes = mask2d[iy, ix]
# mask_nodes_time = np.tile(mask_nodes[None, :], (len(time), 1))  # shape (time, node)

# # Map ETa/ETp to nodes, set NaN where masked
# ETa_nodes = np.where(mask_nodes, ETa_2d[:, iy, ix], np.nan)
# ETp_nodes = np.where(mask_nodes, ETp_2d[:, iy, ix], np.nan)

# # Add to 3D mesh dataset
# ds_hapin["ETa"] = (("time","node"), ETa_nodes)
# ds_hapin["ETp"] =  (("time","node"), ETp_nodes)
# ds_hapin["mask3d"] =  (("time","node"), mask_nodes_time)

# print(ds_hapin)

#%%

# import numpy as np

# # 3D node coordinates
# X_nodes = ds_hapin.mesh3d_nodes.sel(coord="X").values
# Y_nodes = ds_hapin.mesh3d_nodes.sel(coord="Y").values

# # Coarse ET 2D grid and coordinates
# ETa_2d = ds_et_coarse.ETa.values
# ETp_2d = ds_et_coarse.ETp.values
# x_grid = ds_et_coarse.x.values
# y_grid = ds_et_coarse.y.values

# # Coarse grid resolution
# dx = ds_et_coarse.attrs["resolution_x"]
# dy = ds_et_coarse.attrs["resolution_y"]
# x0 = ds_et_coarse.attrs["xllcorner"]
# y0 = ds_et_coarse.attrs["yllcorner"]

# # Convert node XY to raster indices on coarse grid
# ix = ((X_nodes - x0) / dx).astype(int)
# iy = ((Y_nodes - y0) / dy).astype(int)

# # Clip indices to valid range
# ix = np.clip(ix, 0, len(x_grid) - 1)
# iy = np.clip(iy, 0, len(y_grid) - 1)

# # Apply 2D mask from fine grid (optional: could coarsen mask too)
# mask2d = ds_hapin.mask2d.values
# mask_nodes = mask2d[np.clip(iy * 3, 0, hapin["M"]-1), np.clip(ix * 3, 0, hapin["N"]-1)]

# # Map ETa/ETp to nodes, set NaN where masked
# ETa_nodes = np.where(mask_nodes, ETa_2d[iy, ix], np.nan)
# ETp_nodes = np.where(mask_nodes, ETp_2d[iy, ix], np.nan)

# # Add to 3D mesh dataset
# ds_hapin["ETa_coarse"] = (("node",), ETa_nodes)
# ds_hapin["ETp_coarse"] = (("node",), ETp_nodes)
# ds_hapin["mask3d_coarse"] = (("node",), mask_nodes)

# print(ds_hapin)

# simu.grid3d

#%%

# ds_hapin['ETa'].isel(time=0,node=slice(0,int(simu.grid3d['nnod']))).values


#%%


# import matplotlib.pyplot as plt

# # Extract ETa for the first time step
# ETa_nodes = ds_hapin["ETa"].isel(time=0).values  # shape (1271,)
# nnod = 1271

# # Plot ETa vs node index
# plt.figure(figsize=(12,4))
# plt.plot(ds_hapin.node[:nnod], ETa_nodes[0:nnod], marker='o', linestyle='-', markersize=3)
# plt.xlabel('Node index')
# plt.ylabel('ETa')
# plt.title('ETa vs Node')
# plt.grid(True)
# plt.show()


#%%

# import matplotlib.pyplot as plt
# import numpy as np

# # Apply mask: set values outside DEM to NaN
# ETa_masked = np.where(ds_hapin.mask2d.values, ds_et.ETa.isel(time=0).values, np.nan)
# ETp_masked = np.where(ds_hapin.mask2d.values, ds_et.ETp.isel(time=0).values, np.nan)

# # Plot ETa
# plt.figure(figsize=(8, 6))
# plt.imshow(ETa_masked, origin="lower", 
#            extent=[ds_et.x.min(), ds_et.x.max(), ds_et.y.min(), ds_et.y.max()],
#            cmap="viridis", aspect="auto")
# plt.colorbar(label="ETa (mm/day)")
# plt.title("ETa on Structured 2D Grid")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.show()

# # Plot ETp
# plt.figure(figsize=(8, 6))
# plt.imshow(ETp_masked, origin="lower", 
#            extent=[ds_et.x.min(), ds_et.x.max(), ds_et.y.min(), ds_et.y.max()],
#            cmap="plasma", aspect="auto")
# plt.colorbar(label="ETp (mm/day)")
# plt.title("ETp on Structured 2D Grid")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.show()

#%%
# import numpy as np
# from scipy.interpolate import RegularGridInterpolator

# # 3D node coordinates
# X_nodes = ds_hapin.mesh3d_nodes.sel(coord="X").values
# Y_nodes = ds_hapin.mesh3d_nodes.sel(coord="Y").values

# # Coarse ET 2D grid
# ETa_coarse = ds_et_coarse.ETa.values
# ETp_coarse = ds_et_coarse.ETp.values
# x_coarse = ds_et_coarse.x.values
# y_coarse = ds_et_coarse.y.values

# # Build interpolators
# interp_ETa = RegularGridInterpolator((y_coarse, x_coarse), 
#                                      ETa_coarse,
#                                      bounds_error=False, 
#                                      fill_value=np.nan)

# interp_ETp = RegularGridInterpolator((y_coarse, x_coarse), 
#                                      ETp_coarse, 
#                                      bounds_error=False, 
#                                      fill_value=np.nan)

# # Interpolate ET values at 3D node XY positions
# points = np.column_stack([Y_nodes, X_nodes])  # (y, x) for each node
# ETa_nodes = interp_ETa(points)
# ETp_nodes = interp_ETp(points)

# # Apply 3D mask (from fine DEM)
# mask_nodes = ds_hapin.mask3d.values
# ETa_nodes = np.where(mask_nodes, ETa_nodes, np.nan)
# ETp_nodes = np.where(mask_nodes, ETp_nodes, np.nan)

# # Add to 3D mesh dataset
# ds_hapin["ETa_coarse_interp"] = (("node",), ETa_nodes)
# ds_hapin["ETp_coarse_interp"] = (("node",), ETp_nodes)

# print(ds_hapin)

#%%

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np

# # 1️⃣ Fine ET on 2D grid
# plt.figure(figsize=(8, 6))
# ETa_masked_fine = np.where(ds_hapin.mask2d.values, ds_et.isel(time=0).ETa.values, np.nan)
# plt.imshow(ETa_masked_fine, origin="lower",
#            extent=[ds_et.x.min(), ds_et.x.max(), ds_et.y.min(), ds_et.y.max()],
#            cmap="viridis", aspect="auto")
# plt.colorbar(label="ETa (mm/day)")
# plt.title("Fine ETa on 2D Grid")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.show()


# # 2️⃣ Coarse ET on 2D grid
# plt.figure(figsize=(8, 6))
# # Downsample mask to coarse grid
# M_coarse = len(ds_et_coarse.y)
# N_coarse = len(ds_et_coarse.x)
# mask_coarse = ds_hapin.mask2d.values[:M_coarse*3:3, :N_coarse*3:3]
# ETa_masked_coarse = np.where(mask_coarse, ds_et_coarse.ETa.values, np.nan)
# plt.imshow(ETa_masked_coarse, origin="lower",
#            extent=[ds_et_coarse.x.min(), ds_et_coarse.x.max(), ds_et_coarse.y.min(), ds_et_coarse.y.max()],
#            cmap="viridis", aspect="auto")
# plt.colorbar(label="ETa (mm/day)")
# plt.title("Coarse ETa on 2D Grid")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.show()


# # 3️⃣ Interpolated coarse ET on 3D mesh nodes
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection="3d")
# valid = ~np.isnan(ds.ETa_coarse_interp.values)
# sc = ax.scatter(
#     ds.mesh3d_nodes.sel(coord="X").values[valid],
#     ds.mesh3d_nodes.sel(coord="Y").values[valid],
#     ds.mesh3d_nodes.sel(coord="Z").values[valid],
#     c=ds.ETa_coarse_interp.values[valid],
#     cmap="viridis",
#     s=2
# )
# plt.colorbar(sc, label="ETa (mm/day)")
# ax.set_title("Coarse ETa Interpolated on 3D Mesh Nodes")
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")
# plt.show()


#%%

# simu.create_mesh_vtk()

#%%

#%%
