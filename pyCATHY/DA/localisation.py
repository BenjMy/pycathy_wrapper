#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Managing Data Assimilation process localisation. 
    Prepare for DA class
"""
from pyCATHY.cathy_tools import CATHY
import numpy as np 
import pyCATHY.meshtools as mt
from scipy.spatial.distance import cdist


def create_mask_localisation(localisation,veg_map,zones,hapin,grid3d,DEM):
    # Determine the map nodes based on the localisation type
    if localisation == 'veg_map':
        map_surface_nodes = mt.map_cells_to_nodes(veg_map, (hapin['M'] + 1, 
                                                           hapin['N'] + 1)
                                                  )
    elif localisation == 'zones':
        map_surface_nodes = mt.map_cells_to_nodes(zones, (hapin['M'] + 1, 
                                                           hapin['N'] + 1)
                                                  )
    elif localisation == 'nodes':
        try:
            map_surface_nodes = np.arange(grid3d['nnod']).reshape([hapin['M'] + 1,
                                                                    hapin['N'] + 1]
                                                                )  # Treat each node individually   
        except:
            print('Issue here when the mesh contains nan values (boundaries)')
            # # Mask out-of-bound DEM cells
            # mask_out_bound = DEM == -9999
            
            # # Map cells to nodes
            # # mask_out_bound_nodes = ~np.bool_(
            # #     mt.map_cells_to_nodes(mask_out_bound, (hapin['M'] + 1, hapin['N'] + 1))
            # # )
            
            # # Assign unique indices to out-of-bound surface nodes, keep -9999 elsewhere
            # map_surface_nodes = np.full(mask_out_bound.shape, -9999, dtype=int)
            # map_surface_nodes[mask_out_bound] = np.arange(np.sum(mask_out_bound))
            
            map_surface_nodes = np.arange(grid3d['nnod'])


    map_surface_nodes_t = np.hstack(map_surface_nodes)    
    mesh_nodes_valid = mesh_nodes_local_valid(map_surface_nodes,grid3d)    
    return mesh_nodes_valid, map_surface_nodes


    
def mesh_nodes_local_valid(raster_local, grid3d):
    """
    Extrude 2D or 1D valid surface nodes to 3D (depth inclusion).

    Parameters
    ----------
    raster_local : np.ndarray
        2D raster of surface nodes (with -9999 for invalid), OR 1D list of valid node IDs.
    grid3d : dict
        Must contain 'mesh3d_nodes' (N x 3 array of y, x, z coordinates).

    Returns
    -------
    mesh_nodes_valid : list of np.ndarray
        Each element is the list of 3D node IDs corresponding to one surface node.
    """
    print("""[b]
    Extrude surface valid nodes to 3D == Assuming an observation at the surface 
    will influence all the depths!
    [/b]""")

    # Extract **3D grid nodes** from the mesh data
    grid3d = grid3d['mesh3d_nodes']
    num_nodes = grid3d.shape[0]
    node_ids = np.arange(num_nodes)
    x_coords = grid3d[:, 1]  # x values
    y_coords = grid3d[:, 0]  # y values

    mesh_nodes_valid = []

    # Case 1: raster_local is 2D (raster/grid)
    if raster_local.ndim == 2:
        ix = np.unique(x_coords)
        iy = np.unique(y_coords)
        for msirfi in np.unique(raster_local[raster_local != -9999]):
            mask = (raster_local == msirfi)
            idmxi, idmyi = np.where(mask)
            valid_nodes = []
            for idx, idy in zip(idmxi, idmyi):
                valid_node = node_ids[(x_coords == ix[idx]) & (y_coords == iy[idy])]
                valid_nodes.extend(valid_node)
            mesh_nodes_valid.append(np.array(valid_nodes))

    # Case 2: raster_local is 1D (already list of valid surface nodes)
    elif raster_local.ndim == 1:
        for msirfi in raster_local[raster_local != -9999]:
            # directly map surface node (x, y) to all depths
            valid_nodes = node_ids[(x_coords == x_coords[int(msirfi)]) & (y_coords == y_coords[int(msirfi)])]
            mesh_nodes_valid.append(valid_nodes)
    else:
        raise ValueError("raster_local must be 1D or 2D")

    return mesh_nodes_valid



def build_localization_EM_matrix(grid_coords, ncoils, L, with_coil_covariance=True):
    """
    Build covariance localization matrix for a 2D grid with multiple coils per grid point.
    
    Parameters:
    - grid_coords: (N x 2) array of grid coordinates [(x1, y1), (x2, y2), ...]
    - ncoils: number of coils per grid point
    - L: Gaspari-Cohn localization radius
    - with_coil_covariance: if True, include covariance between coils; else treat coils as independent
    
    Returns:
    - localization_matrix: (N*ncoils) x (N*ncoils) array
    """
    # Observation coordinates including coils
    coil_idx = np.arange(ncoils)
    obs_coords = np.array([[gx, gy, c] for gx, gy in grid_coords for c in coil_idx])
    
    if with_coil_covariance:
        # Full 3D distance including coil index
        dist_matrix = cdist(obs_coords, obs_coords)
        localization_matrix = gaspari_cohn(dist_matrix, L)
    else:
        # Only spatial distance, coils independent
        dist_matrix = cdist(obs_coords[:, :2], obs_coords[:, :2])
        localization_matrix = np.kron(np.eye(ncoils), gaspari_cohn(dist_matrix, L))
    
    return localization_matrix



def gaspari_cohn(r, L):
    """
    Gaspari-Cohn localization function
    r: distance (can be array)
    L: localization radius
    """
    r = np.abs(r) / L
    w = np.zeros_like(r)
    
    mask1 = r <= 1
    mask2 = (r > 1) & (r <= 2)
    
    w[mask1] = (((-0.25 * r[mask1] + 0.5) * r[mask1] + 0.625) * r[mask1] - 5/3) * r[mask1]**2 + 1
    w[mask2] = ((((r[mask2] / 12 - 0.5) * r[mask2] + 0.625) * r[mask2] + 5/3) * r[mask2] - 5) * r[mask2] + 4 - 2/(3*r[mask2])
    w[r > 2] = 0
    return w


def build_localization_matrix(obs_pos, grid_pos, L):
    """
    Build localization matrix using Gaspari-Cohn function.
    
    Parameters
    ----------
    obs_pos : np.array of shape (n_obs, 2)
        Observation positions [[x1, y1], [x2, y2], ...]
    grid_pos : np.array of shape (n_grid, 2)
        Grid positions [[x1, y1], [x2, y2], ...]
    L : float
        Localization radius
        
    Returns
    -------
    np.array of shape (n_grid, n_obs)
        Localization matrix
    """
    n_grid = grid_pos.shape[0]
    n_obs = obs_pos.shape[0]
    
    loc_matrix = np.zeros((n_grid, n_obs))
    
    for i in range(n_grid):
        # distances from grid point i to all observations
        dx = obs_pos[:,0] - grid_pos[i,0]
        dy = obs_pos[:,1] - grid_pos[i,1]
        d = np.sqrt(dx**2 + dy**2)
        
        loc_matrix[i,:] = gaspari_cohn(d, L)
    
    return loc_matrix



