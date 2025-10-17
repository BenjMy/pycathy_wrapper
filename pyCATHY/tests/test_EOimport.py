#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

# Import necessary libraries
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv

# Import pyCATHY modules for handling mesh, inputs, and outputs
import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as cplt

#%% Initialize CATHY model
# Define the project directory and model name. This example uses 'atmbc_spatially_from_weill'.
path2prj = "../SSHydro/"  # Replace with your local project path
simu = cathy_tools.CATHY(dirName=path2prj, prj_name="atmbc_spatially_from_weill_withnodata")
figpath = "../results/DA_ET_test/"  # Path to store figures/results

#%% Load and manipulate the DEM (Digital Elevation Model)
# Read the DEM input file
DEM, dem_header = simu.read_inputs('dem')

# Create a new DEM array filled with ones and add irregular boundary and invalid values (-9999)
DEM_new = np.ones(np.shape(DEM))  # Initialize new DEM with ones
DEM_new[-1, -1] = 1 - 1e-3  # Adjust a specific corner value
DEM_new[10:20, 0:10] = -9999  # Add an interior block of invalid values to simulate an irregular boundary
DEM_new[0:3, 15:20] = -9999  # Add an interior block of invalid values to simulate an irregular boundary

# Update the CATHY inputs with the modified DEM
simu.update_prepo_inputs(DEM_new)

# Visualize the updated DEM
simu.show_input('dem')

#%% Preprocess and mesh generation
# Run the preprocessor to handle inputs and generate the mesh
simu.run_preprocessor()

# Create a 3D mesh visualization (VTK format)
simu.create_mesh_vtk(verbose=True)

# Load the 3D grid output
grid3d = simu.read_outputs('grid3d')

# Set parameters for elevation
simu.dem_parameters
elevation_increment = 0.5 / 21  # Define elevation increment per row
elevation_matrix = np.ones([21, 21])  # Initialize the elevation matrix

# Populate elevation_matrix with incremental values based on row index
for row in range(21):
    elevation_matrix[row, :] += row * elevation_increment

#%% Define spatially variable atmospheric boundary condition inputs

# # Set up time intervals and cycles for the boundary condition
# interval = 5  # Number of intervals
# ncycles = 7   # Number of cycles
# t_atmbc = np.linspace(1e-3, 36e3 * ncycles, interval * ncycles)  # Time vector

# # Atmospheric boundary condition value
# v_atmbc_value = -2e-7  # Set the boundary condition value

# # Check if the number of nodes matches the flattened elevation matrix
# if int(grid3d['nnod']) == len(np.ravel(elevation_matrix)):
#     # Calculate the atmospheric boundary condition for each node based on elevation
#     v_atmbc = np.ones(int(grid3d['nnod'])) * v_atmbc_value * np.ravel(elevation_matrix)
# else:
#     # For cases where the number of nodes doesn't match, calculate for all nodes
#     v_atmbc_all_nodes = np.ones(len(np.ravel(elevation_matrix))) * v_atmbc_value * np.ravel(np.exp(elevation_matrix**2))

#     # Reshape the boundary condition values to match the DEM shape
#     v_atmbc_mat = np.reshape(v_atmbc_all_nodes, [np.shape(simu.DEM)[0] + 1, np.shape(simu.DEM)[0] + 1])

#     # Mask invalid values in the DEM (-9999) by setting them to NaN
#     maskDEM_novalid = np.where(DEM_new == -9999)
#     v_atmbc_mat[maskDEM_novalid] = np.nan

#     # Flatten the masked matrix and remove NaN values
#     v_atmbc = np.ravel(v_atmbc_mat)
#     v_atmbc = v_atmbc[~np.isnan(v_atmbc)]  # Use ~np.isnan to filter out NaN values

# # Visualize the spatial variation of the atmospheric boundary condition
# fig, ax = plt.subplots()
# img = ax.imshow(v_atmbc_mat)
# plt.colorbar(img)

# # Update the atmospheric boundary condition (ATMB) parameters in CATHY
# simu.update_atmbc(
#     HSPATM=0,
#     IETO=0,
#     time=t_atmbc,
#     netValue=[v_atmbc] * len(t_atmbc)  # Apply the same boundary condition at all times
# )

# Update the model parameters (time control) in CATHY
# simu.update_parm(TIMPRTi=t_atmbc)


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

# ---------------------------------------------------------------------
# 1. Create ET dataset on structured grid
# ---------------------------------------------------------------------
def create_et_dataset(hapin, coarse_factor=1):
    """
    Create synthetic ETa and ETp dataset on a structured grid.
    Can create at full resolution or coarser (by coarse_factor).
    """
    N = hapin["N"] // coarse_factor
    M = hapin["M"] // coarse_factor
    dx = hapin["delta_x"] * coarse_factor
    dy = hapin["delta_y"] * coarse_factor
    x0 = hapin["xllcorner"]
    y0 = hapin["yllcorner"]

    # Cell centers
    x = x0 + (np.arange(N) + 0.5) * dx
    y = y0 + (np.arange(M) + 0.5) * dy
    X, Y = np.meshgrid(x, y)

    # Example patterns (replace with your data later)
    ETa = 3 + 2 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())
    ETp = 4 + 1.5 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())

    return xr.Dataset(
        data_vars={"ETa": (("y", "x"), ETa),
                   "ETp": (("y", "x"), ETp)},
        coords={"x": ("x", x), "y": ("y", y)},
        attrs={"resolution_x": dx, "resolution_y": dy,
               "xllcorner": x0, "yllcorner": y0}
    )


# ---------------------------------------------------------------------
# 2. Interpolate ET values from structured grid to 3D mesh nodes
# ---------------------------------------------------------------------
import numpy as np

def map_et_to_mesh_from_grid(ds_et, ds_mesh, fine_mask2d, hapin):
    """
    Map ETa/ETp from any-resolution 2D grid to 3D mesh nodes using direct index mapping.
    
    Parameters
    ----------
    ds_et : xarray.Dataset
        ET dataset with dims ('y','x') and variables ETa, ETp
    ds_mesh : xarray.Dataset
        Mesh dataset with mesh3d_nodes (dims: node, coord)
    fine_mask2d : np.ndarray
        Boolean mask of fine resolution (dims MxN)
    hapin : dict
        Fine grid info (delta_x, delta_y, N, M, xllcorner, yllcorner)
    
    Returns
    -------
    ds_mesh : xarray.Dataset
        Mesh dataset updated with ETa/ETp and mask3d
    """
    # Node coordinates
    X_nodes = ds_mesh.mesh3d_nodes.sel(coord="X").values
    Y_nodes = ds_mesh.mesh3d_nodes.sel(coord="Y").values

    # ET values and grid
    ETa_2d = ds_et.ETa.values
    ETp_2d = ds_et.ETp.values
    x_grid = ds_et.x.values
    y_grid = ds_et.y.values

    # Grid resolution and origin
    dx = ds_et.attrs["resolution_x"]
    dy = ds_et.attrs["resolution_y"]
    x0 = ds_et.attrs["xllcorner"]
    y0 = ds_et.attrs["yllcorner"]

    # Convert node XY to raster indices
    ix = ((X_nodes - x0) / dx).astype(int)
    iy = ((Y_nodes - y0) / dy).astype(int)

    # Clip indices
    ix = np.clip(ix, 0, len(x_grid) - 1)
    iy = np.clip(iy, 0, len(y_grid) - 1)

    # Map fine mask to node positions
    scale_x = hapin["N"] // len(x_grid)
    scale_y = hapin["M"] // len(y_grid)
    mask_nodes = fine_mask2d[np.clip(iy * scale_y, 0, hapin["M"]-1),
                              np.clip(ix * scale_x, 0, hapin["N"]-1)]

    # Map ET values to nodes, NaN where masked
    ETa_nodes = np.where(mask_nodes, ETa_2d[iy, ix], np.nan)
    ETp_nodes = np.where(mask_nodes, ETp_2d[iy, ix], np.nan)

    # Assign to dataset
    ds_mesh["ETa"] = (("node",), ETa_nodes)
    ds_mesh["ETp"] = (("node",), ETp_nodes)
    ds_mesh["mask3d"] = (("node",), mask_nodes)

    return ds_mesh


import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

import numpy as np
import xarray as xr

def build_mesh_dataset(hapin, grid3d):
    """
    Build an xarray Dataset containing:
      - structured grid coordinates from hapin
      - unstructured mesh (nodes + tetrahedra) from grid3d
    
    Parameters
    ----------
    hapin : dict
        Dictionary with structured grid info:
          {"delta_x", "delta_y", "N", "M", "N_celle", "xllcorner", "yllcorner"}
    grid3d : dict
        Dictionary with 3D mesh info:
          {"nnod", "nnod3", "nel", "mesh3d_nodes", "mesh_tetra"}
    
    Returns
    -------
    xr.Dataset
        Dataset containing mesh nodes, mesh elements, structured coordinates.
    """
    # Structured grid coordinates
    x = hapin['xllcorner'] + np.arange(hapin['N']) * hapin['delta_x']
    y = hapin['yllcorner'] + np.arange(hapin['M']) * hapin['delta_y']

    # Wrap into xarray
    ds = xr.Dataset(
        data_vars={
            "mesh3d_nodes": (("node", "coord"), grid3d["mesh3d_nodes"]),
            "mesh_tetra": (("element", "nodes_per_element"), grid3d["mesh_tetra"]),
        },
        coords={
            "x": ("x", x),
            "y": ("y", y),
            "node": np.arange(grid3d["mesh3d_nodes"].shape[0]),
            "coord": ["X", "Y", "Z"],
            "element": np.arange(grid3d["mesh_tetra"].shape[0]),
            "nodes_per_element": np.arange(grid3d["mesh_tetra"].shape[1]),
        },
        attrs={
            "delta_x": hapin["delta_x"],
            "delta_y": hapin["delta_y"],
            "N": hapin["N"],
            "M": hapin["M"],
            "N_celle": hapin["N_celle"],
            "nnod": grid3d["nnod"],
            "nnod3": grid3d["nnod3"],
            "nel": grid3d["nel"],
        }
    )
    return ds



# ---------- Mask creation: 2D ----------
def make_mask2d_from_raster(raster_DEM_masked, hapin, nodata=-9999):
    """
    Create a boolean 2D mask (y,x) from raster_DEM_masked.
    - raster_DEM_masked: 2D numpy array with shape (M, N) or compatible.
    - hapin: dictionary with 'M','N','delta_x','delta_y','xllcorner','yllcorner'.
    Returns: xr.DataArray mask2d (dims: y,x) with cell-center coords.
    """
    M = hapin["M"]
    N = hapin["N"]
    if raster_DEM_masked.shape != (M, N):
        raise ValueError(f"raster_DEM_masked.shape {raster_DEM_masked.shape} != (M,N)=({M},{N}). "
                         "Make sure your raster matches hapin grid dimensions.")
    mask = raster_DEM_masked != nodata
    # Use cell centers for coords (consistent with ds_et creation)
    dx = hapin["delta_x"]
    dy = hapin["delta_y"]
    x0 = hapin["xllcorner"]
    y0 = hapin["yllcorner"]
    x_centers = x0 + (np.arange(N) + 0.5) * dx
    y_centers = y0 + (np.arange(M) + 0.5) * dy

    mask_da = xr.DataArray(mask.astype(bool),
                           coords={"y": ("y", y_centers), "x": ("x", x_centers)},
                           dims=("y", "x"),
                           name="mask2d")
    return mask_da


# ---------- Mask creation: 3D from 2D ----------
def make_mask3d_from_mask2d(ds_mesh, mask2d_da, hapin=None, method="interp", threshold=0.5):
    """
    Build a boolean mask for mesh nodes based on 2D mask.
    Parameters:
      - ds_mesh: xarray Dataset containing mesh3d_nodes with coord labels ['X','Y','Z'] (dim 'node').
      - mask2d_da: xr.DataArray boolean mask with dims ('y','x') and coords y,x (cell centers).
      - hapin: needed only for method='nearest' to compute indices. If method='interp', hapin optional.
      - method: 'interp' (interpolate mask, robust) or 'nearest' (nearest-cell lookup using hapin).
      - threshold: for interp, values >= threshold considered True.
    Returns:
      - xr.DataArray mask3d (dims: node)
    """
    X_nodes = ds_mesh.mesh3d_nodes.sel(coord="X").values
    Y_nodes = ds_mesh.mesh3d_nodes.sel(coord="Y").values
    if method == "interp":
        # interpolate mask (float) then threshold
        interp = RegularGridInterpolator((mask2d_da.y.values, mask2d_da.x.values),
                                         mask2d_da.values.astype(float),
                                         bounds_error=False, fill_value=0.0)
        pts = np.column_stack([Y_nodes, X_nodes])  # (y,x)
        vals = interp(pts)   # floats in [0,1]
        mask_nodes = vals >= float(threshold)
    elif method == "nearest":
        if hapin is None:
            raise ValueError("hapin must be provided for method='nearest'")
        dx = hapin["delta_x"]
        dy = hapin["delta_y"]
        x0 = hapin["xllcorner"]
        y0 = hapin["yllcorner"]
        # mask2d coords are cell centers -> x0 + (i+0.5)*dx
        # compute index via rounding
        ix = np.round((X_nodes - (x0 + 0.5 * dx)) / dx).astype(int)
        iy = np.round((Y_nodes - (y0 + 0.5 * dy)) / dy).astype(int)
        # clip
        ix = np.clip(ix, 0, mask2d_da.sizes["x"] - 1)
        iy = np.clip(iy, 0, mask2d_da.sizes["y"] - 1)
        mask_nodes = mask2d_da.values[iy, ix].astype(bool)
    else:
        raise ValueError("method must be 'interp' or 'nearest'")

    mask3d_da = xr.DataArray(mask_nodes,
                             coords={"node": ("node", ds_mesh.node.values)},
                             dims=("node",),
                             name="mask3d")
    return mask3d_da


# ---------- Attach masks to dataset ----------
def attach_masks_to_mesh_dataset(ds_mesh, mask2d_da, mask3d_da, names=("mask2d", "mask3d")):
    """
    Attach provided mask DataArrays to ds_mesh as variables.
    Returns ds_mesh (modified).
    """
    ds_mesh[names[0]] = (("y", "x"), mask2d_da.values) if isinstance(mask2d_da, xr.DataArray) else mask2d_da
    ds_mesh[names[1]] = (("node",), mask3d_da.values) if isinstance(mask3d_da, xr.DataArray) else mask3d_da
    return ds_mesh



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# -------- Plot ET on 2D structured grid --------
def plot_et2d(ds_map, var="ETa", cmap="viridis"):
    """
    Plot ET variable on 2D structured grid using ds_map.
    - ds_map: xarray Dataset with dims (y, x) and mask2d
    - var: "ETa" or "ETp"
    """
    # Reconstruct 2D field from node values
    # Use average of nodes that belong to each grid cell (quick demo with mask2d)
    data2d = np.where(ds_map["mask2d"].values, np.nan, np.nan)  # placeholder
    # for simplicity we just plot mask2d for now
    # if ETa/ETp per cell needed, must remap nodes->grid cells (optional)

    plt.figure(figsize=(6, 5))
    im = plt.pcolormesh(ds_map.x, ds_map.y, ds_map["mask2d"], 
                        shading="auto", cmap=cmap)
    plt.colorbar(im, label=var)
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.title(f"{var} on 2D grid")
    plt.axis("equal")
    plt.show()


# -------- Plot ET on 3D mesh --------
def plot_et3d(ds_map, var="ETa", cmap="viridis", s=15, alpha=0.7):
    """
    Scatter plot ET variable on 3D mesh nodes.
    - ds_map: xarray Dataset with mesh3d_nodes, ETa/ETp, mask3d
    - var: "ETa" or "ETp"
    """
    coords = ds_map.mesh3d_nodes
    X = coords.sel(coord="X").values
    Y = coords.sel(coord="Y").values
    Z = coords.sel(coord="Z").values
    vals = ds_map[var].values
    mask = ds_map["mask3d"].values

    # Apply mask
    X, Y, Z, vals = X[mask], Y[mask], Z[mask], vals[mask]

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(X, Y, Z, c=vals, cmap=cmap, s=s, alpha=alpha)
    fig.colorbar(sc, ax=ax, label=var)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title(f"{var} mapped on 3D mesh")
    plt.show()


# raster_DEM_masked

raster_DEM_masked = simu.DEM 

ds_meshCATHY = build_mesh_dataset(simu.hapin, grid3d)

# make a mask2d from raster_DEM_masked
mask2d_da = make_mask2d_from_raster(raster_DEM_masked,
                                    simu.hapin, 
                                    nodata=-9999
                                    )

# make mask3d (interpolated)
mask3d_da = make_mask3d_from_mask2d(ds_meshCATHY, mask2d_da, method="interp")

# attach to ds_mesh
ds_meshCATHY = attach_masks_to_mesh_dataset(ds_meshCATHY, mask2d_da, mask3d_da)

# hapin = simu.hapin
# # 1. Create fine and coarse ET datasets
ds_et = create_et_dataset(simu.hapin, coarse_factor=1)
ds_et_coarse = create_et_dataset(simu.hapin, coarse_factor=3)


ds_et['ETa'].plot.imshow()


# # 2. Map coarse ET dataset to 3D mesh nodes
ds_map = map_et_to_mesh_from_grid(ds_et,
                                  ds_meshCATHY, 
                                  ds_meshCATHY.mask2d.values,
                                  simu.hapin
                                  )


# # 3. Plot results
# plot_et_fields(ds, ds_et, ds_et_coarse)

# 2D plot
plot_et2d(ds_map, 
          var="ETa"
          )

# 3D plot
# plot_et3d(ds_map, var="ETa")

#%%

import numpy as np
import xarray as xr

# Example input
hapin = simu.hapin

grid3d = simu.grid3d

# Build coordinate vectors from hapin
x = hapin['xllcorner'] + np.arange(hapin['N']) * hapin['delta_x']
y = hapin['yllcorner'] + np.arange(hapin['M']) * hapin['delta_y']

# Wrap into xarray
ds = xr.Dataset(
    data_vars={
        "mesh3d_nodes": (("node", "coord"), grid3d["mesh3d_nodes"]),
        "mesh_tetra": (("element", "nodes_per_element"), grid3d["mesh_tetra"]),
    },
    coords={
        "x": ("x", x),
        "y": ("y", y),
        "node": np.arange(grid3d["mesh3d_nodes"].shape[0]),
        "coord": ["X", "Y", "Z"],
        "element": np.arange(grid3d["mesh_tetra"].shape[0]),
        "nodes_per_element": np.arange(grid3d["mesh_tetra"].shape[1]),
    },
    attrs={
        "delta_x": hapin["delta_x"],
        "delta_y": hapin["delta_y"],
        "N": hapin["N"],
        "M": hapin["M"],
        "N_celle": hapin["N_celle"],
        "nnod": grid3d["nnod"],
        "nnod3": grid3d["nnod3"],
        "nel": grid3d["nel"],
    }
)

print(ds)


import matplotlib.pyplot as plt

# ---- 2D structured grid from hapin ----
fig, ax = plt.subplots(figsize=(6, 5))
ax.set_title("Structured 2D grid (hapin)")
ax.set_aspect("equal")

# Plot grid lines
for xval in ds.x.values:
    ax.plot([xval, xval], [ds.y.values.min(), ds.y.values.max()], color="lightgrey", lw=0.5)
for yval in ds.y.values:
    ax.plot([ds.x.values.min(), ds.x.values.max()], [yval, yval], color="lightgrey", lw=0.5)

ax.set_xlabel("X")
ax.set_ylabel("Y")

plt.show()


# ---- 3D mesh nodes scatterplot ----
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.set_title("3D mesh nodes")

X = ds.mesh3d_nodes.sel(coord="X").values
Y = ds.mesh3d_nodes.sel(coord="Y").values
Z = ds.mesh3d_nodes.sel(coord="Z").values

ax.scatter(X, Y, Z, s=2, c=Z, cmap="viridis")

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()


# raster_DEM_masked : shape (M, N) with -9999 for no data
mask = raster_DEM_masked != -9999  # True where valid

# Add to dataset
ds["mask"] = (("y", "x"), mask)

ds.mask.plot.imshow()

# 2D mask from DEM
mask2d = raster_DEM_masked != -9999  # shape (M, N)
ds["mask2d"] = (("y", "x"), mask2d)


hapin = simu.hapin
# Extract parameters from hapin
dx = hapin["delta_x"]
dy = hapin["delta_y"]
x0 = hapin["xllcorner"]
y0 = hapin["yllcorner"]

# Node coordinates
X = ds.mesh3d_nodes.sel(coord="X").values
Y = ds.mesh3d_nodes.sel(coord="Y").values

# Convert to indices
ix = ((X - x0) / dx).astype(int)
iy = ((Y - y0) / dy).astype(int)

# Clip indices inside raster bounds
ix = np.clip(ix, 0, hapin["N"] - 1)
iy = np.clip(iy, 0, hapin["M"] - 1)

# Build mask for nodes
mask3d = mask2d[iy, ix]   # boolean mask for each node
ds["mask3d"] = (("node",), mask3d)

# 2D mask
ds.mask2d.plot.imshow()

# 3D nodes, masked
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.set_title("3D mesh nodes (masked)")

X = ds.mesh3d_nodes.sel(coord="X").values
Y = ds.mesh3d_nodes.sel(coord="Y").values
Z = ds.mesh3d_nodes.sel(coord="Z").values

valid = ds.mask3d.values
ax.scatter(X[valid], Y[valid], Z[valid], s=2, c=Z[valid], cmap="viridis")

plt.show()

import numpy as np
import xarray as xr

# Grid info
N = hapin["N"]
M = hapin["M"]
dx = hapin["delta_x"]
dy = hapin["delta_y"]
x0 = hapin["xllcorner"]
y0 = hapin["yllcorner"]



import numpy as np
import xarray as xr

hapin = simu.hapin


# Cell centers
x = x0 + (np.arange(N) + 0.5) * dx
y = y0 + (np.arange(M) + 0.5) * dy

# Create meshgrid for smooth patterns
X, Y = np.meshgrid(x, y)

# Smooth ETa / ETp patterns (replace with real data if available)
ETa = 3 + 2 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())  # 2D gradient
ETp = 4 + 1.5 * np.sin(np.pi * X / X.max()) * np.cos(np.pi * Y / Y.max())

# Create xarray Dataset
ds_et = xr.Dataset(
    data_vars={
        "ETa": (("y", "x"), ETa),
        "ETp": (("y", "x"), ETp),
    },
    coords={
        "x": ("x", x),
        "y": ("y", y),
    },
    attrs={
        "resolution_x": dx,
        "resolution_y": dy,
        "xllcorner": x0,
        "yllcorner": y0,
    }
)

# Quick visualization
ds_et.ETa.plot.imshow(cmap="viridis")
ds_et.ETp.plot.imshow(cmap="plasma")

import numpy as np

# 3D node coordinates
X_nodes = ds.mesh3d_nodes.sel(coord="X").values
Y_nodes = ds.mesh3d_nodes.sel(coord="Y").values

# ET 2D grid and coordinates
ETa_2d = ds_et.ETa.values
ETp_2d = ds_et.ETp.values
x_grid = ds_et.x.values
y_grid = ds_et.y.values

# Convert node XY to raster indices
dx = hapin["delta_x"]
dy = hapin["delta_y"]
x0 = hapin["xllcorner"]
y0 = hapin["yllcorner"]

ix = ((X_nodes - x0) / dx).astype(int)
iy = ((Y_nodes - y0) / dy).astype(int)

# Clip indices to valid range
ix = np.clip(ix, 0, hapin["N"] - 1)
iy = np.clip(iy, 0, hapin["M"] - 1)

# Apply 2D mask
mask2d = ds.mask2d.values
mask_nodes = mask2d[iy, ix]

# Map ETa/ETp to nodes, set NaN where masked
ETa_nodes = np.where(mask_nodes, ETa_2d[iy, ix], np.nan)
ETp_nodes = np.where(mask_nodes, ETp_2d[iy, ix], np.nan)

# Add to 3D mesh dataset
ds["ETa"] = (("node",), ETa_nodes)
ds["ETp"] = (("node",), ETp_nodes)
ds["mask3d"] = (("node",), mask_nodes)

print(ds)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
valid = ~np.isnan(ds.ETa.values)

sc = ax.scatter(
    ds.mesh3d_nodes.sel(coord="X").values[valid],
    ds.mesh3d_nodes.sel(coord="Y").values[valid],
    ds.mesh3d_nodes.sel(coord="Z").values[valid],
    c=ds.ETa.values[valid],
    cmap="viridis",
    s=2
)
plt.colorbar(sc, label="ETa (mm/day)")
plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Apply mask: set values outside DEM to NaN
ETa_masked = np.where(ds.mask2d.values, ds_et.ETa.values, np.nan)
ETp_masked = np.where(ds.mask2d.values, ds_et.ETp.values, np.nan)

# Plot ETa
plt.figure(figsize=(8, 6))
plt.imshow(ETa_masked, origin="lower", 
           extent=[ds_et.x.min(), ds_et.x.max(), ds_et.y.min(), ds_et.y.max()],
           cmap="viridis", aspect="auto")
plt.colorbar(label="ETa (mm/day)")
plt.title("ETa on Structured 2D Grid")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

# Plot ETp
plt.figure(figsize=(8, 6))
plt.imshow(ETp_masked, origin="lower", 
           extent=[ds_et.x.min(), ds_et.x.max(), ds_et.y.min(), ds_et.y.max()],
           cmap="plasma", aspect="auto")
plt.colorbar(label="ETp (mm/day)")
plt.title("ETp on Structured 2D Grid")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()


#%%

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def process_grid_and_mesh(simu, raster_DEM_masked, generate_et_data=True, plot_results=True):
    """
    Comprehensive function to process 2D structured grids and 3D unstructured meshes.
    
    Parameters:
    -----------
    simu : object
        Simulation object containing hapin and grid3d attributes
    raster_DEM_masked : ndarray
        2D array with DEM data, -9999 for no data values
    generate_et_data : bool, default True
        Whether to generate synthetic ET data or use existing data
    plot_results : bool, default True
        Whether to generate visualization plots
        
    Returns:
    --------
    ds : xarray.Dataset
        Dataset containing 3D mesh data with interpolated values
    ds_et : xarray.Dataset  
        Dataset containing 2D structured grid ET data
    """
    
    # Extract input data
    hapin = simu.hapin
    grid3d = simu.grid3d
    
    # ========== BUILD 2D STRUCTURED GRID ==========
    
    # Grid parameters
    N = hapin["N"]
    M = hapin["M"]
    dx = hapin["delta_x"]
    dy = hapin["delta_y"]
    x0 = hapin["xllcorner"]
    y0 = hapin["yllcorner"]
    
    # Build coordinate vectors (cell centers for ET data)
    x_centers = x0 + (np.arange(N) + 0.5) * dx
    y_centers = y0 + (np.arange(M) + 0.5) * dy
    
    # Build coordinate vectors (grid lines for visualization)
    x_lines = x0 + np.arange(N) * dx
    y_lines = y0 + np.arange(M) * dy
    
    # ========== CREATE MASKS ==========
    
    # 2D mask from DEM (True where valid data exists)
    mask2d = raster_DEM_masked != -9999  # shape (M, N)
    
    # ========== BUILD 3D MESH DATASET ==========
    
    # Create main dataset with 3D mesh
    ds = xr.Dataset(
        data_vars={
            "mesh3d_nodes": (("node", "coord"), grid3d["mesh3d_nodes"]),
            "mesh_tetra": (("element", "nodes_per_element"), grid3d["mesh_tetra"]),
            "mask2d": (("y", "x"), mask2d),
        },
        coords={
            "x": ("x", x_centers),
            "y": ("y", y_centers),
            "x_lines": ("x_lines", x_lines),
            "y_lines": ("y_lines", y_lines),
            "node": np.arange(grid3d["mesh3d_nodes"].shape[0]),
            "coord": ["X", "Y", "Z"],
            "element": np.arange(grid3d["mesh_tetra"].shape[0]),
            "nodes_per_element": np.arange(grid3d["mesh_tetra"].shape[1]),
        },
        attrs={
            "delta_x": dx,
            "delta_y": dy,
            "N": N,
            "M": M,
            "N_celle": hapin["N_celle"],
            "nnod": grid3d["nnod"],
            "nnod3": grid3d["nnod3"],
            "nel": grid3d["nel"],
            "xllcorner": x0,
            "yllcorner": y0,
        }
    )
    
    # ========== CREATE 3D NODE MASK ==========
    
    # Extract 3D node coordinates
    X_nodes = ds.mesh3d_nodes.sel(coord="X").values
    Y_nodes = ds.mesh3d_nodes.sel(coord="Y").values
    Z_nodes = ds.mesh3d_nodes.sel(coord="Z").values
    
    # Convert node coordinates to grid indices
    ix = ((X_nodes - x0) / dx).astype(int)
    iy = ((Y_nodes - y0) / dy).astype(int)
    
    # Clip indices to valid range
    ix = np.clip(ix, 0, N - 1)
    iy = np.clip(iy, 0, M - 1)
    
    # Create 3D node mask
    mask3d = mask2d[iy, ix]
    ds["mask3d"] = (("node",), mask3d)
    
    # ========== GENERATE OR PROCESS ET DATA ==========
    
    if generate_et_data:
        # Create meshgrid for synthetic patterns
        X_grid, Y_grid = np.meshgrid(x_centers, y_centers)
        
        # Generate synthetic ET patterns
        ETa = 3 + 2 * np.sin(np.pi * X_grid / X_grid.max()) * np.cos(np.pi * Y_grid / Y_grid.max())
        ETp = 4 + 1.5 * np.sin(np.pi * X_grid / X_grid.max()) * np.cos(np.pi * Y_grid / Y_grid.max())
    else:
        # Placeholder for real ET data - replace with actual data loading
        ETa = np.ones((M, N)) * 2.5  # Default values
        ETp = np.ones((M, N)) * 3.5
    
    # Create ET dataset
    ds_et = xr.Dataset(
        data_vars={
            "ETa": (("y", "x"), ETa),
            "ETp": (("y", "x"), ETp),
        },
        coords={
            "x": ("x", x_centers),
            "y": ("y", y_centers),
        },
        attrs={
            "resolution_x": dx,
            "resolution_y": dy,
            "xllcorner": x0,
            "yllcorner": y0,
        }
    )
    
    # ========== INTERPOLATE ET DATA TO 3D NODES ==========
    
    # Map 2D ET values to 3D nodes
    ETa_nodes = np.where(mask3d, ETa[iy, ix], np.nan)
    ETp_nodes = np.where(mask3d, ETp[iy, ix], np.nan)
    
    # Add ET data to 3D mesh dataset
    ds["ETa"] = (("node",), ETa_nodes)
    ds["ETp"] = (("node",), ETp_nodes)
    
    # ========== VISUALIZATION ==========
    
    if plot_results:
        # Create subplot layout
        fig = plt.figure(figsize=(16, 12))
        
        # 1. 2D structured grid
        ax1 = plt.subplot(2, 3, 1)
        ax1.set_title("Structured 2D Grid")
        ax1.set_aspect("equal")
        
        # Plot grid lines
        for xval in x_lines:
            ax1.plot([xval, xval], [y_lines.min(), y_lines.max()], color="lightgrey", lw=0.5)
        for yval in y_lines:
            ax1.plot([x_lines.min(), x_lines.max()], [yval, yval], color="lightgrey", lw=0.5)
        
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        
        # 2. 2D mask
        ax2 = plt.subplot(2, 3, 2)
        ax2.set_title("2D Mask from DEM")
        im2 = ax2.imshow(mask2d, origin="lower", cmap="viridis", aspect="auto",
                        extent=[x_centers.min(), x_centers.max(), y_centers.min(), y_centers.max()])
        plt.colorbar(im2, ax=ax2)
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        
        # 3. ETa on 2D grid (masked)
        ax3 = plt.subplot(2, 3, 3)
        ax3.set_title("ETa on 2D Grid")
        ETa_masked = np.where(mask2d, ETa, np.nan)
        im3 = ax3.imshow(ETa_masked, origin="lower", cmap="viridis", aspect="auto",
                        extent=[x_centers.min(), x_centers.max(), y_centers.min(), y_centers.max()])
        plt.colorbar(im3, ax=ax3, label="ETa (mm/day)")
        ax3.set_xlabel("X")
        ax3.set_ylabel("Y")
        
        # 4. ETp on 2D grid (masked)
        ax4 = plt.subplot(2, 3, 4)
        ax4.set_title("ETp on 2D Grid")
        ETp_masked = np.where(mask2d, ETp, np.nan)
        im4 = ax4.imshow(ETp_masked, origin="lower", cmap="plasma", aspect="auto",
                        extent=[x_centers.min(), x_centers.max(), y_centers.min(), y_centers.max()])
        plt.colorbar(im4, ax=ax4, label="ETp (mm/day)")
        ax4.set_xlabel("X")
        ax4.set_ylabel("Y")
        
        # 5. 3D mesh nodes (all)
        ax5 = plt.subplot(2, 3, 5, projection="3d")
        ax5.set_title("3D Mesh Nodes (All)")
        sc5 = ax5.scatter(X_nodes, Y_nodes, Z_nodes, s=2, c=Z_nodes, cmap="viridis")
        ax5.set_xlabel("X")
        ax5.set_ylabel("Y")
        ax5.set_zlabel("Z")
        
        # 6. 3D mesh nodes with ETa (masked)
        ax6 = plt.subplot(2, 3, 6, projection="3d")
        ax6.set_title("3D Mesh Nodes with ETa")
        valid_nodes = ~np.isnan(ETa_nodes)
        sc6 = ax6.scatter(X_nodes[valid_nodes], Y_nodes[valid_nodes], Z_nodes[valid_nodes],
                         c=ETa_nodes[valid_nodes], cmap="viridis", s=2)
        plt.colorbar(sc6, ax=ax6, label="ETa (mm/day)", shrink=0.6)
        ax6.set_xlabel("X")
        ax6.set_ylabel("Y")
        ax6.set_zlabel("Z")
        
        plt.tight_layout()
        plt.show()
        
        # Summary statistics
        print("\n" + "="*50)
        print("PROCESSING SUMMARY")
        print("="*50)
        print(f"2D Grid dimensions: {M} x {N}")
        print(f"Grid resolution: dx={dx}, dy={dy}")
        print(f"3D mesh nodes: {len(X_nodes)}")
        print(f"3D mesh elements: {grid3d['mesh_tetra'].shape[0]}")
        print(f"Valid 2D cells: {mask2d.sum()} / {mask2d.size} ({100*mask2d.sum()/mask2d.size:.1f}%)")
        print(f"Valid 3D nodes: {mask3d.sum()} / {len(mask3d)} ({100*mask3d.sum()/len(mask3d):.1f}%)")
        print(f"ETa range (2D): {np.nanmin(ETa_masked):.2f} - {np.nanmax(ETa_masked):.2f}")
        print(f"ETp range (2D): {np.nanmin(ETp_masked):.2f} - {np.nanmax(ETp_masked):.2f}")
        print(f"ETa range (3D nodes): {np.nanmin(ETa_nodes):.2f} - {np.nanmax(ETa_nodes):.2f}")
        print(f"ETp range (3D nodes): {np.nanmin(ETp_nodes):.2f} - {np.nanmax(ETp_nodes):.2f}")
    
    return ds, ds_et


# Example usage:
ds, ds_et = process_grid_and_mesh(simu, raster_DEM_masked, generate_et_data=True, plot_results=True)
# print(ds)
# print(d

#%%

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RegularGridInterpolator

def process_misaligned_grid_and_mesh(simu, raster_DEM_masked, et_params=None, plot_results=True):
    """
    Process 2D structured grids and 3D unstructured meshes with misaligned ET data.
    
    The ET data can have:
    - Different spatial resolution (pixel size)
    - Different coordinate system origin
    - Different spatial extent
    - Rotated or shifted grid alignment
    
    Parameters:
    -----------
    simu : object
        Simulation object containing hapin and grid3d attributes
    raster_DEM_masked : ndarray
        2D array with DEM data, -9999 for no data values
    et_params : dict, optional
        ET grid parameters. If None, uses default misaligned parameters:
        {
            'resolution': 50.0,  # Different pixel size
            'x_offset': 1200,    # X coordinate shift
            'y_offset': -800,    # Y coordinate shift  
            'nx': 80,           # Different number of pixels
            'ny': 60,
            'rotation': 5.0     # Rotation in degrees (future enhancement)
        }
    plot_results : bool, default True
        Whether to generate visualization plots
        
    Returns:
    --------
    ds : xarray.Dataset
        Dataset containing 3D mesh data with interpolated ET values
    ds_et : xarray.Dataset  
        Dataset containing misaligned 2D ET data
    interpolation_info : dict
        Information about the interpolation process
    """
    
    # Extract input data
    hapin = simu.hapin
    grid3d = simu.grid3d
    
    # ========== BUILD ORIGINAL 3D MESH DATASET ==========
    
    # Grid parameters for the main simulation
    N = hapin["N"]
    M = hapin["M"]
    dx = hapin["delta_x"]
    dy = hapin["delta_y"]
    x0 = hapin["xllcorner"]
    y0 = hapin["yllcorner"]
    
    # Build coordinate vectors for main grid
    x_main = x0 + (np.arange(N) + 0.5) * dx
    y_main = y0 + (np.arange(M) + 0.5) * dy
    
    # 2D mask from DEM
    mask2d = raster_DEM_masked != -9999
    
    # Create main 3D mesh dataset
    ds = xr.Dataset(
        data_vars={
            "mesh3d_nodes": (("node", "coord"), grid3d["mesh3d_nodes"]),
            "mesh_tetra": (("element", "nodes_per_element"), grid3d["mesh_tetra"]),
            "mask2d": (("y", "x"), mask2d),
        },
        coords={
            "x": ("x", x_main),
            "y": ("y", y_main),
            "node": np.arange(grid3d["mesh3d_nodes"].shape[0]),
            "coord": ["X", "Y", "Z"],
            "element": np.arange(grid3d["mesh_tetra"].shape[0]),
            "nodes_per_element": np.arange(grid3d["mesh_tetra"].shape[1]),
        },
        attrs={
            "delta_x": dx,
            "delta_y": dy,
            "N": N,
            "M": M,
            "xllcorner": x0,
            "yllcorner": y0,
        }
    )
    
    # ========== SET UP MISALIGNED ET GRID PARAMETERS ==========
    
    # Calculate main grid extent
    main_x_extent = x_main.max() - x_main.min()
    main_y_extent = y_main.max() - y_main.min()
    main_x_center = (x_main.min() + x_main.max()) / 2
    main_y_center = (y_main.min() + y_main.max()) / 2
    
    if et_params is None:
        # Default parameters - ET grid always bigger than main grid
        et_params = {
            'resolution': max(dx, dy) * 2.0,     # 2x coarser than main grid
            'x_offset': main_x_extent * 0.1,     # 10% shift right
            'y_offset': -main_y_extent * 0.05,   # 5% shift down
            'coverage_factor': 1.5,              # ET grid 1.5x larger than main grid
            'buffer_cells': 2                    # Minimum buffer around main grid
        }
    
    # ET grid parameters
    et_dx = et_params['resolution']
    et_dy = et_params['resolution']
    
    # Calculate required ET grid size to fully cover main grid + buffer
    coverage_factor = et_params.get('coverage_factor', 1.5)
    buffer_cells = et_params.get('buffer_cells', 2)
    
    # Minimum required ET grid extent to cover main grid
    required_x_extent = main_x_extent * coverage_factor
    required_y_extent = main_y_extent * coverage_factor
    
    # ET grid dimensions (ensure minimum buffer)
    et_nx = max(int(np.ceil(required_x_extent / et_dx)), 
                int(np.ceil(main_x_extent / et_dx)) + 2 * buffer_cells)
    et_ny = max(int(np.ceil(required_y_extent / et_dy)), 
                int(np.ceil(main_y_extent / et_dy)) + 2 * buffer_cells)
    
    # Override with user-specified dimensions if provided (but warn if too small)
    if 'nx' in et_params:
        min_nx = int(np.ceil(main_x_extent / et_dx)) + 2 * buffer_cells
        if et_params['nx'] < min_nx:
            print(f"Warning: Specified nx={et_params['nx']} may be too small. Minimum recommended: {min_nx}")
        et_nx = et_params['nx']
    
    if 'ny' in et_params:
        min_ny = int(np.ceil(main_y_extent / et_dy)) + 2 * buffer_cells
        if et_params['ny'] < min_ny:
            print(f"Warning: Specified ny={et_params['ny']} may be too small. Minimum recommended: {min_ny}")
        et_ny = et_params['ny']
    
    # Calculate actual ET grid extent
    et_x_extent = et_nx * et_dx
    et_y_extent = et_ny * et_dy
    
    # Position ET grid to fully contain main grid with specified offset
    # Start by centering, then apply offset
    et_x0 = main_x_center - (et_x_extent / 2) + et_params['x_offset']
    et_y0 = main_y_center - (et_y_extent / 2) + et_params['y_offset']
    
    # Verify complete coverage (adjust if needed)
    main_x_min, main_x_max = x_main.min(), x_main.max()
    main_y_min, main_y_max = y_main.min(), y_main.max()
    et_x_min, et_x_max = et_x0, et_x0 + et_x_extent
    et_y_min, et_y_max = et_y0, et_y0 + et_y_extent
    
    # Adjust ET grid position if main grid extends outside
    if main_x_min < et_x_min:
        et_x0 = main_x_min - et_dx * buffer_cells
    elif main_x_max > et_x_max:
        et_x0 = main_x_max - et_x_extent + et_dx * buffer_cells
        
    if main_y_min < et_y_min:
        et_y0 = main_y_min - et_dy * buffer_cells
    elif main_y_max > et_y_max:
        et_y0 = main_y_max - et_y_extent + et_dy * buffer_cells
    
    # Build ET coordinate vectors
    x_et = et_x0 + (np.arange(et_nx) + 0.5) * et_dx
    y_et = et_y0 + (np.arange(et_ny) + 0.5) * et_dy
    
    # ========== GENERATE ET DATA ON MISALIGNED GRID ==========
    
    # Create meshgrid for ET data
    X_et, Y_et = np.meshgrid(x_et, y_et)
    
    # Generate ET patterns with different characteristics than main grid
    # Use different wavelengths and phases to show misalignment effects
    ETa_et = (3.5 + 2.5 * np.sin(2 * np.pi * X_et / (main_x_extent * 0.7)) * 
              np.cos(2 * np.pi * Y_et / (main_y_extent * 0.8)) +
              1.0 * np.sin(4 * np.pi * X_et / main_x_extent))
    
    ETp_et = (4.2 + 1.8 * np.cos(2 * np.pi * X_et / (main_x_extent * 0.9)) * 
              np.sin(2 * np.pi * Y_et / (main_y_extent * 0.6)) +
              0.8 * np.cos(3 * np.pi * Y_et / main_y_extent))
    
    # Create ET dataset with misaligned coordinates
    ds_et = xr.Dataset(
        data_vars={
            "ETa": (("y_et", "x_et"), ETa_et),
            "ETp": (("y_et", "x_et"), ETp_et),
        },
        coords={
            "x_et": ("x_et", x_et),
            "y_et": ("y_et", y_et),
        },
        attrs={
            "resolution_x": et_dx,
            "resolution_y": et_dy,
            "xllcorner": et_x0,
            "yllcorner": et_y0,
            "nx": et_nx,
            "ny": et_ny,
            "alignment": "misaligned_with_main_grid",
        }
    )
    
    # ========== INTERPOLATE ET DATA TO 3D MESH NODES ==========
    
    # Extract 3D node coordinates
    X_nodes = ds.mesh3d_nodes.sel(coord="X").values
    Y_nodes = ds.mesh3d_nodes.sel(coord="Y").values
    Z_nodes = ds.mesh3d_nodes.sel(coord="Z").values
    
    # Create interpolators for ET data
    # Note: RegularGridInterpolator expects (y, x) order for coordinates
    eta_interpolator = RegularGridInterpolator(
        points=(y_et, x_et),
        values=ETa_et,
        method='linear',
        bounds_error=False,
        fill_value=np.nan
    )
    
    etp_interpolator = RegularGridInterpolator(
        points=(y_et, x_et),
        values=ETp_et,
        method='linear',
        bounds_error=False,
        fill_value=np.nan
    )
    
    # Interpolate ET values at 3D node locations
    # Note: RegularGridInterpolator expects points as (y, x) pairs
    node_points = np.column_stack((Y_nodes, X_nodes))
    
    ETa_nodes_interp = eta_interpolator(node_points)
    ETp_nodes_interp = etp_interpolator(node_points)
    
    # Apply original mask (nodes outside DEM area)
    # Convert node coordinates to main grid indices for masking
    ix_main = ((X_nodes - x0) / dx).astype(int)
    iy_main = ((Y_nodes - y0) / dy).astype(int)
    ix_main = np.clip(ix_main, 0, N - 1)
    iy_main = np.clip(iy_main, 0, M - 1)
    mask3d = mask2d[iy_main, ix_main]
    
    # Apply mask to interpolated values
    ETa_nodes_final = np.where(mask3d, ETa_nodes_interp, np.nan)
    ETp_nodes_final = np.where(mask3d, ETp_nodes_interp, np.nan)
    
    # Add interpolated data to 3D mesh dataset
    ds["ETa"] = (("node",), ETa_nodes_final)
    ds["ETp"] = (("node",), ETp_nodes_final)
    ds["mask3d"] = (("node",), mask3d)
    
    # ========== INTERPOLATION INFORMATION ==========
    
    interpolation_info = {
        'main_grid_resolution': (dx, dy),
        'et_grid_resolution': (et_dx, et_dy),
        'resolution_ratio': (et_dx/dx, et_dy/dy),
        'spatial_offset': (et_params['x_offset'], et_params['y_offset']),
        'main_grid_extent': (x_main.min(), x_main.max(), y_main.min(), y_main.max()),
        'et_grid_extent': (x_et.min(), x_et.max(), y_et.min(), y_et.max()),
        'overlap_extent': (
            max(x_main.min(), x_et.min()), 
            min(x_main.max(), x_et.max()),
            max(y_main.min(), y_et.min()), 
            min(y_main.max(), y_et.max())
        ),
        'nodes_with_et_data': np.sum(~np.isnan(ETa_nodes_final)),
        'total_nodes': len(X_nodes),
        'interpolation_method': 'linear'
    }
    
    # ========== VISUALIZATION ==========
    
    if plot_results:
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Main grid with mask
        ax1 = plt.subplot(3, 4, 1)
        ax1.set_title("Main Grid Mask")
        im1 = ax1.imshow(mask2d, origin="lower", cmap="viridis", aspect="auto",
                        extent=[x_main.min(), x_main.max(), y_main.min(), y_main.max()])
        plt.colorbar(im1, ax=ax1)
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        
        # 2. ET grid ETa (misaligned)
        ax2 = plt.subplot(3, 4, 2)
        ax2.set_title("ET Grid - ETa (Misaligned)")
        im2 = ax2.imshow(ETa_et, origin="lower", cmap="viridis", aspect="auto",
                        extent=[x_et.min(), x_et.max(), y_et.min(), y_et.max()])
        plt.colorbar(im2, ax=ax2, label="ETa (mm/day)")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        
        # 3. ET grid ETp (misaligned)
        ax3 = plt.subplot(3, 4, 3)
        ax3.set_title("ET Grid - ETp (Misaligned)")
        im3 = ax3.imshow(ETp_et, origin="lower", cmap="plasma", aspect="auto",
                        extent=[x_et.min(), x_et.max(), y_et.min(), y_et.max()])
        plt.colorbar(im3, ax=ax3, label="ETp (mm/day)")
        ax3.set_xlabel("X")
        ax3.set_ylabel("Y")
        
        # 4. Grid comparison (overlay)
        ax4 = plt.subplot(3, 4, 4)
        ax4.set_title("Grid Comparison")
        # Main grid outline
        ax4.plot([x_main.min(), x_main.max(), x_main.max(), x_main.min(), x_main.min()],
                [y_main.min(), y_main.min(), y_main.max(), y_main.max(), y_main.min()],
                'b-', linewidth=2, label='Main Grid')
        # ET grid outline
        ax4.plot([x_et.min(), x_et.max(), x_et.max(), x_et.min(), x_et.min()],
                [y_et.min(), y_et.min(), y_et.max(), y_et.max(), y_et.min()],
                'r--', linewidth=2, label='ET Grid')
        ax4.set_xlabel("X")
        ax4.set_ylabel("Y")
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        # 5. 3D mesh nodes (all)
        ax5 = plt.subplot(3, 4, 5, projection="3d")
        ax5.set_title("3D Mesh Nodes (All)")
        sc5 = ax5.scatter(X_nodes, Y_nodes, Z_nodes, s=1, c=Z_nodes, cmap="viridis", alpha=0.6)
        ax5.set_xlabel("X")
        ax5.set_ylabel("Y")
        ax5.set_zlabel("Z")
        
        # 6. 3D mesh nodes with interpolated ETa
        ax6 = plt.subplot(3, 4, 6, projection="3d")
        ax6.set_title("3D Nodes - Interpolated ETa")
        valid_eta = ~np.isnan(ETa_nodes_final)
        if np.any(valid_eta):
            sc6 = ax6.scatter(X_nodes[valid_eta], Y_nodes[valid_eta], Z_nodes[valid_eta],
                             c=ETa_nodes_final[valid_eta], cmap="viridis", s=2)
            plt.colorbar(sc6, ax=ax6, label="ETa (mm/day)", shrink=0.6)
        ax6.set_xlabel("X")
        ax6.set_ylabel("Y")
        ax6.set_zlabel("Z")
        
        # 7. 3D mesh nodes with interpolated ETp
        ax7 = plt.subplot(3, 4, 7, projection="3d")
        ax7.set_title("3D Nodes - Interpolated ETp")
        valid_etp = ~np.isnan(ETp_nodes_final)
        if np.any(valid_etp):
            sc7 = ax7.scatter(X_nodes[valid_etp], Y_nodes[valid_etp], Z_nodes[valid_etp],
                             c=ETp_nodes_final[valid_etp], cmap="plasma", s=2)
            plt.colorbar(sc7, ax=ax7, label="ETp (mm/day)", shrink=0.6)
        ax7.set_xlabel("X")
        ax7.set_ylabel("Y")
        ax7.set_zlabel("Z")
        
        # 8. Interpolation coverage map
        ax8 = plt.subplot(3, 4, 8)
        ax8.set_title("Interpolation Coverage")
        # Show which nodes have valid interpolated data
        coverage_map = np.zeros((M, N))
        valid_interp = ~np.isnan(ETa_nodes_interp)  # Before masking
        for i, (x, y, valid) in enumerate(zip(X_nodes, Y_nodes, valid_interp)):
            ix = int(np.clip((x - x0) / dx, 0, N-1))
            iy = int(np.clip((y - y0) / dy, 0, M-1))
            if valid:
                coverage_map[iy, ix] = 1
        im8 = ax8.imshow(coverage_map, origin="lower", cmap="RdYlGn", aspect="auto",
                        extent=[x_main.min(), x_main.max(), y_main.min(), y_main.max()])
        plt.colorbar(im8, ax=ax8, label="Has ET Data")
        ax8.set_xlabel("X")
        ax8.set_ylabel("Y")
        
        
        plt.tight_layout()
        plt.show()
        
        # Print detailed summary
        print("\n" + "="*70)
        print("MISALIGNED GRID PROCESSING SUMMARY")
        print("="*70)
        print(f"Main Grid:")
        print(f"  - Resolution: {dx:.1f} x {dy:.1f}")
        print(f"  - Dimensions: {N} x {M}")
        print(f"  - Extent: X=[{x_main.min():.0f}, {x_main.max():.0f}], Y=[{y_main.min():.0f}, {y_main.max():.0f}]")
        print(f"\nET Grid:")
        print(f"  - Resolution: {et_dx:.1f} x {et_dy:.1f}")
        print(f"  - Dimensions: {et_nx} x {et_ny}")
        print(f"  - Extent: X=[{x_et.min():.0f}, {x_et.max():.0f}], Y=[{y_et.min():.0f}, {y_et.max():.0f}]")
        print(f"  - Offset: X={et_params['x_offset']:.1f}, Y={et_params['y_offset']:.1f}")
        print(f"  - Coverage factor: {et_params.get('coverage_factor', 1.0):.1f}")
        print(f"  - Buffer cells: {et_params.get('buffer_cells', 0)}")
        
        # Verify coverage
        coverage_x = (x_et.max() - x_et.min()) / (x_main.max() - x_main.min())
        coverage_y = (y_et.max() - y_et.min()) / (y_main.max() - y_main.min())
        print(f"  - Actual coverage ratio: {coverage_x:.2f}x (X), {coverage_y:.2f}x (Y)")
        
        # Check if main grid is fully contained
        fully_contained = (x_et.min() <= x_main.min() and x_et.max() >= x_main.max() and
                          y_et.min() <= y_main.min() and y_et.max() >= y_main.max())
        print(f"  - Main grid fully contained: {'✓ Yes' if fully_contained else '✗ No'}")
        print(f"\nInterpolation Results:")
        print(f"  - Total 3D nodes: {len(X_nodes)}")
        print(f"  - Nodes with ET data: {interpolation_info['nodes_with_et_data']} ({100*interpolation_info['nodes_with_et_data']/len(X_nodes):.1f}%)")
        print(f"  - ETa range: {np.nanmin(ETa_nodes_final):.2f} - {np.nanmax(ETa_nodes_final):.2f}")
        print(f"  - ETp range: {np.nanmin(ETp_nodes_final):.2f} - {np.nanmax(ETp_nodes_final):.2f}")
        print(f"  - Resolution ratio (ET/Main): {interpolation_info['resolution_ratio'][0]:.2f} x {interpolation_info['resolution_ratio'][1]:.2f}")
    
    return ds, ds_et, interpolation_info


# Example usage with different scenarios:

def create_example_scenarios(main_grid_info):
    """
    Create example scenarios where ET grid is always larger than main grid
    
    Parameters:
    -----------
    main_grid_info : dict
        Dictionary with keys: 'dx', 'dy', 'x_extent', 'y_extent', 'x_center', 'y_center'
    """
    
    dx, dy = main_grid_info['dx'], main_grid_info['dy']
    x_extent, y_extent = main_grid_info['x_extent'], main_grid_info['y_extent']
    
    scenarios = {
        'fine_shifted': {
            'resolution': max(dx, dy) * 0.8,    # Finer resolution
            'x_offset': x_extent * 0.1,         # 10% shift right
            'y_offset': -y_extent * 0.05,       # 5% shift down
            'coverage_factor': 1.4,             # 40% larger coverage
            'buffer_cells': 3,                  # 3-cell buffer
        },
        
        'coarse_centered': {
            'resolution': max(dx, dy) * 3.0,    # Much coarser resolution
            'x_offset': 0.0,                    # Centered
            'y_offset': 0.0,                    # Centered
            'coverage_factor': 2.0,             # Double coverage
            'buffer_cells': 2,                  # 2-cell buffer
        },
        
        'medium_offset': {
            'resolution': max(dx, dy) * 1.5,    # Moderately coarser
            'x_offset': x_extent * 0.15,        # 15% shift right
            'y_offset': y_extent * 0.1,         # 10% shift up
            'coverage_factor': 1.6,             # 60% larger coverage
            'buffer_cells': 4,                  # 4-cell buffer
        },
        
        'large_coarse': {
            'resolution': max(dx, dy) * 4.0,    # Very coarse resolution
            'x_offset': -x_extent * 0.05,       # 5% shift left
            'y_offset': -y_extent * 0.08,       # 8% shift down
            'coverage_factor': 2.5,             # 2.5x larger coverage
            'buffer_cells': 1,                  # Minimal buffer
        },
        
        'minimal_coverage': {
            'resolution': max(dx, dy) * 2.0,    # 2x coarser
            'x_offset': x_extent * 0.2,         # 20% shift (will be adjusted)
            'y_offset': -y_extent * 0.15,       # 15% shift (will be adjusted)
            'coverage_factor': 1.1,             # Minimal coverage (just bigger)
            'buffer_cells': 1,                  # Minimal buffer
        }
    }
    
    return scenarios

def run_misalignment_example(simu, raster_DEM_masked, scenario_name='default'):
    """
    Run the misaligned grid processing with proper parameter scaling
    
    Parameters:
    -----------
    simu : object
        Simulation object
    raster_DEM_masked : ndarray
        DEM data
    scenario_name : str
        One of: 'default', 'fine_shifted', 'coarse_offset', 'partial_overlap', 'minimal_overlap'
    """
    
    # Extract main grid info for proper scaling
    hapin = simu.hapin
    dx, dy = hapin["delta_x"], hapin["delta_y"]
    x0, y0 = hapin["xllcorner"], hapin["yllcorner"]
    N, M = hapin["N"], hapin["M"]
    
    x_extent = N * dx
    y_extent = M * dy
    x_center = x0 + x_extent / 2
    y_center = y0 + y_extent / 2
    
    main_grid_info = {
        'dx': dx, 'dy': dy,
        'x_extent': x_extent, 'y_extent': y_extent,
        'x_center': x_center, 'y_center': y_center
    }
    
    # Get scenario parameters
    if scenario_name == 'default':
        et_params = None  # Use automatic defaults
    else:
        scenarios = create_example_scenarios(main_grid_info)
        if scenario_name not in scenarios:
            print(f"Unknown scenario '{scenario_name}'. Available: {list(scenarios.keys())}")
            et_params = None
        else:
            et_params = scenarios[scenario_name]
    
    print(f"\nRunning scenario: {scenario_name}")
    print(f"Main grid: {N}x{M}, resolution: {dx}x{dy}, extent: {x_extent:.1f}x{y_extent:.1f}")
    
    if et_params:
        print(f"ET params: res={et_params['resolution']:.1f}, offset=({et_params['x_offset']:.1f}, {et_params['y_offset']:.1f})")
    
    # Run the processing
    ds, ds_et, info = process_misaligned_grid_and_mesh(simu, raster_DEM_masked, 
                                                       et_params=et_params, 
                                                       plot_results=True)
    
    return ds, ds_et, info

# Example usage:

# Simple usage with automatic parameter scaling
ds, ds_et, info = run_misalignment_example(simu, raster_DEM_masked, 'default')

# Try different scenarios
ds, ds_et, info = run_misalignment_example(simu, raster_DEM_masked, 'fine_shifted')
ds, ds_et, info = run_misalignment_example(simu, raster_DEM_masked, 'coarse_offset') 
ds, ds_et, info = run_misalignment_example(simu, raster_DEM_masked, 'partial_overlap')

# Manual parameter specification
custom_params = {
    'resolution': 1.0,      # Specify resolution directly
    'x_offset': 2.0,        # Specify offsets directly  
    'y_offset': -1.5,
    'nx': 15,               # Grid dimensions
    'ny': 12
}
ds, ds_et, info = process_misaligned_grid_and_mesh(simu, raster_DEM_masked,
                                                   et_params=custom_params)