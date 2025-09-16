#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Managing Data Assimilation process localisation. 
    Prepare for DA class
"""
from pyCATHY.cathy_tools import CATHY
import numpy as np 
import pyCATHY.meshtools as mt
from scipy.spatial.distance import cdist


def create_mask_localisation(localisation,veg_map,zones,hapin,grid3d):
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
        map_surface_nodes = np.arange(grid3d['nnod']).reshape([hapin['M'] + 1,
                                                                hapin['N'] + 1]
                                                            )  # Treat each node individually          

    map_surface_nodes_t = np.hstack(map_surface_nodes)    
    
    
    # n = int(np.sqrt(map_surface_nodes.size))  # assume square
    # map_surface_nodes_2d = np.arange(1, n*n+1).reshape(n, n)
         
    mesh_nodes_valid = mesh_nodes_local_valid(map_surface_nodes,grid3d)
    # np.shape(mesh_nodes_valid)
    
    return mesh_nodes_valid, map_surface_nodes


    
def mesh_nodes_local_valid(raster_local,grid3d):
    '''
    Extrude 2d valid node to 3d (depth inclusion)

    Parameters
    ----------
    raster_local : TYPE
        DESCRIPTION.
    grid3d : TYPE
        DESCRIPTION.

    Returns
    -------
    mesh_nodes_valid : TYPE
        DESCRIPTION.

    '''
    print("""[b]
            Extrude 2d valid node to 3d == Assuming that an observation at the surface will influence all the depths!
            [/b]"""
            )
                            
    # Extract **3D grid nodes** from the mesh data
    grid3d = grid3d['mesh3d_nodes']
    num_nodes = grid3d.shape[0]
    node_ids = np.arange(num_nodes)
    x_coords = grid3d[:, 1]  # x values
    y_coords = grid3d[:, 0]  # y values
    ix = np.unique(x_coords)
    iy = np.unique(y_coords)
    mesh_nodes_valid = []
    # Iterate over unique surface node identifiers
    for msirfi in np.unique(raster_local):
        # Get indices where raster_local equals the current identifier
        mask = (raster_local == msirfi)
        # Get the indices of the valid nodes
        idmxi, idmyi = np.where(mask)
        # Filter for valid node IDs based on x and y coordinates
        valid_nodes = []
        for idx, idy in zip(idmxi, idmyi):
            valid_node = node_ids[(x_coords == ix[idx]) & (y_coords == iy[idy])]
            valid_nodes.extend(valid_node)  # Extend list with valid nodes
        mesh_nodes_valid.append(np.array(valid_nodes))
    return mesh_nodes_valid


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



def build_localization_matrix(grid_coords, ncoils, L, with_coil_covariance=True):
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
