#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Managing Data Assimilation process localisation. 
    Prepare for DA class
"""
from pyCATHY.cathy_tools import CATHY
import numpy as np 
import pyCATHY.meshtools as mt


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
    mesh_nodes_valid = mesh_nodes_local_valid(map_surface_nodes,grid3d)
    
    return mesh_nodes_valid, map_surface_nodes


    
def mesh_nodes_local_valid(raster_local,grid3d):
    # Extract 3D grid nodes from the mesh data
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
