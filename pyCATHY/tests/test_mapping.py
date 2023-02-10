#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:23:58 2022

@author: ben
"""

import os
import sys

import pyCATHY
from pyCATHY import cathy_tools

sys.path.append(
    os.path.join("/home/ben/Documents/GitHub/BenjMy/Re-della-Pietra/", "lib")
)

import numpy as np
import pyvista as pv
import utils_CATHY
import vtk

from pyCATHY.importers import cathy_inputs as in_CT

# import vtkchange

# vtk.__version__

#%%
simu = cathy_tools.CATHY(dirName="./data/", prj_name="mapping")  # args.scenario

#%%

simu.run_processor(IPRT1=3, verbose=True)

simu.grid3d


len(simu.grid3d["nodes_idxyz"])
#%%


#%%

simu.update_zone()
zone_names = ["unique_zone"]
layers_depths = [[0, 1]]
simu.zone

simu.update_prepo_inputs()

simu.dem_parameters

dem, header_dem = simu.read_inputs("dem")
meshpv = pv.read("./mapping/vtk/mapping.vtk")

#%%

simu.hapin
# simu.DEM

# from pyvista import examples
# demex = examples.download_crater_topo()
# terrain = demex.warp_by_scalar()
# terrain.plot()


# dem

grid = pv.UniformGrid()


xu = np.unique(meshpv.points[:, 0])
yu = np.unique(meshpv.points[:, 1])

dempar = simu.dem_parameters["zratio(i),i=1,nstr"].split("\t")
dempar_ratio = [float(d) for d in dempar]

# zu = [-np.cumsum(dempar_ratio[:i+1])[-1]*simu.dem_parameters['base'] for i in range(simu.dem_parameters['nstr'])]
zu = 0

x, y, z = np.meshgrid(xu, yu, zu)


# grid = pv.UniformGrid((simu.hapin['N'],
#                        simu.hapin['M'],
#                        1)
#                       )
grid = pv.StructuredGrid(x, y, z)

# grid.cell_data_to_point_data()
# grid['dem']

len(np.ravel(dem))
# grid.add_field_data(np.ravel(dem),'dem')

grid["dem"] = np.ravel(dem.T) - 1
# grid.set_active_scalars('dem')
# grid.plot(show_edges=(True))

grid = grid.point_data_to_cell_data()
# grid.plot()
grid = grid.cell_data_to_point_data()

terrain = grid.warp_by_scalar()
# terrain.plot()


# z_cells = np.array([25] * 5 + [35] * 3 + [50] * 2 + [75, 100])


dempar = simu.dem_parameters["zratio(i),i=1,nstr"].split("\t")
dempar_ratio = [float(d) for d in dempar]
z_cells = [
    -np.cumsum(dempar_ratio[: i + 1])[-1] * simu.dem_parameters["base"]
    for i in range(simu.dem_parameters["nstr"])
]

xx = np.repeat(terrain.x, len(z_cells), axis=-1)
yy = np.repeat(terrain.y, len(z_cells), axis=-1)
zz = (
    np.repeat(terrain.z, len(z_cells), axis=-1) - z_cells
)  # np.cumsum(z_cells).reshape((1, 1, -1))


min(np.ravel(terrain.z))


mesh = pv.StructuredGrid(xx, yy, -zz)
mesh["Elevation"] = zz.ravel(order="F")
mesh

# mesh.plot(show_edges=True, lighting=False)


# p = pv.Plotter(notebook=False)
# p.add_mesh(mesh,show_edges=True)
# p.add_bounding_box()
# p.show_grid()
# p.show()


# mesh_unstr = mesh.cast_to_unstructured_grid()


# p = pv.Plotter(notebook=False)
# p.add_mesh(mesh_unstr,show_edges=True)
# p.add_bounding_box()
# p.show_grid()
# p.show()


p = pv.Plotter(notebook=False)
p.add_mesh(meshpv, show_edges=True)
p.add_mesh(mesh)
p.add_bounding_box()
p.show_grid()
p.show()


#%%
# extend to 3d the zone raster file
# -----------------------------------
zones3d = utils_CATHY.zone3d(simu)


#%%

# insert layers flag into the 3d the zone raster file
# -----------------------------------
zones3d_layered = utils_CATHY.create_layers_inzones3d(
    simu, zones3d, zone_names, layers_depths
)

#%%

# plot 3d zones files layered
# ------------------------------------------
utils_CATHY.plot_zones3d_layered(simu, zones3d_layered)

#%%


# meshpv2.plot()
# index = meshpv.find_cells_within_bounds([-2.0, 2.0, -2.0, 2.0, -2.0, 2.0])

#%%

# test = zones3d_layered[1:,:,:]
# np.shape(test)

# mesh_unstr['zones'] =  np.ravel(zones3d_layered[:,:,:])
# # grid.cell_data['my_array'] = np.ravel(zones3d_layered[1:,:,:])

# mesh_unstr.set_active_scalars('zones')
# mesh_unstr = mesh.point_data_to_cell_data()

# mesh_unstr.plot(show_edges=True)
# result = grid.interpolate(meshpv)
# result = meshpv.interpolate(mesh_unstr,
#                             )

# result.plot(show_edges=True)


#%%
# meshpv.plot(show_edges=1)

xu = np.unique(meshpv.points[:, 0])
yu = np.unique(meshpv.points[:, 1])

dempar = simu.dem_parameters["zratio(i),i=1,nstr"].split("\t")
dempar_ratio = [float(d) for d in dempar]

zu = [
    -np.cumsum(dempar_ratio[: i + 1])[-1] * simu.dem_parameters["base"]
    for i in range(simu.dem_parameters["nstr"])
]


len(xu)
len(yu)
np.shape(zones3d_layered)

x, y, z = np.meshgrid(xu, yu, zu)


grid = pv.StructuredGrid(x, y, z)

grid_cell_centers = grid.cell_centers()

len(np.unique(grid_cell_centers.points[:, 0]))
len(np.unique(grid_cell_centers.points[:, 1]))
len(np.unique(grid_cell_centers.points[:, 2]))


np.random.random(grid.n_cells)

grid.cell_data["my_array"] = np.ravel(zones3d_layered[1:, :, :])

# result = grid.interpolate(meshpv)
result = meshpv.interpolate(grid)
# result.plot()

# result.set_active_scalars('my_array')
# result.plot(show_edges=1)

result.save("test.vtk", binary=False)


grid["flagZ"] = zones3d_layered[1:, :, :]
pl = pv.Plotter()
actor = pl.add_mesh(grid, show_edges=True)
actor = pl.add_points(
    grid_cell_centers, render_points_as_spheres=True, color="red", point_size=20
)
pl.show()


grid.plot(show_edges=1)
p = pv.Plotter(notebook=True)
p.add_mesh(grid)


# meshpv_ptsdata = meshpv.cell_data_to_point_data()
#%%

from pyvista import examples

grid = examples.load_explicit_structured()
grid["zones"] = np.r_[np.ones(4), np.ones(206) * 2]
grid.set_active_scalars("zones")
grid = grid.point_data_to_cell_data()
grid.plot(show_edges=True, show_bounds=True)


#%%

# Starting from the generated mesh
# -------------------------------------

# test = zones3d_layered[1:,:,:]
# meshpv_struc = meshpv.cast_to_explicit_structured_grid()
# meshpv_struc.plot(show_edges=1)
# cell_centers = meshpv.cell_centers()


# pl = pv.Plotter()
# actor = pl.add_mesh(meshpv, show_edges=True)
# actor = pl.add_points(cell_centers, render_points_as_spheres=True,
#                       color='red', point_size=20)
# pl.show()


#%%
# utils_CATHY.get_layers_depth(simu)

# simu.update_ic(INDP=0,
#                pressure_head_ini=-5
#                )
