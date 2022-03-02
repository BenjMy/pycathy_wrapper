#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 18:03:04 2021

@author: ben
"""

import os 
import numpy as np
import pyvista as pv
import pygimli as pg
import pygimli.meshtools as mt


def trace_mesh_pg(meshIN,meshOUT,method='spline', **kwargs):
    
    # Specify interpolation method 'linear, 'spline', 'harmonic'
    meshIN = pg.load(meshIN)
    meshOUT = pg.load(meshOUT)

    out_data = pg.interpolate(meshIN['ER_converted_CATHYmsh*'], meshOUT, method=method)
    
    return out_data
    

def trace_mesh(meshIN,meshOUT,scalar,threshold=1e-1,**kwargs):
    '''
    Trace meshIN on meshOUT using nearest neigbour interpolation

    Parameters
    ----------
    meshIN : TYPE
        DESCRIPTION.
    meshOUT : TYPE
        DESCRIPTION.
    threshold : TYPE, optional
        DESCRIPTION. The default is 1e-1.

    Returns
    -------
    out_data : TYPE
        DESCRIPTION.
    '''
    
    print(np.mean(meshIN.get_array(scalar)))
    # print(np.mean(meshIN.get_array(scalar)))

    in_nodes_mod = np.array(meshIN.points)
    if 'in_nodes_mod' in kwargs:
        in_nodes_mod = kwargs['in_nodes_mod']

    # out_data = []
    # closest_idx = []

    # if type(meshOUT) is str:
    #     meshOUT = pv.read(meshOUT)


    # cellOUT_centers = meshOUT.cell_centers()

    # for pos in cellOUT_centers.points:

    #     closest_idx, closest = find_nearest_node(pos,in_nodes_mod,
    #                                                   threshold)

    #     if 'nan' in closest:
    #         out_data.append('nan')
    #     else:
    #         out_data.append(meshIN.active_scalars[closest_idx])
    
    # out_data = np.hstack(out_data)
    # len(out_data)   
    
    # print(meshIN.active_scalars_info)
    
    meshIN.set_active_scalars(scalar)
    
    meshIN.points = in_nodes_mod
    
    rd= min([abs(min(np.diff(meshIN.points[:,0]))),
         abs(min(np.diff(meshIN.points[:,1]))),
         abs(min(np.diff(meshIN.points[:,2])))
         ]
        )

    
    result = meshOUT.interpolate(meshIN, radius=rd*25, pass_point_data=True)
    result = result.point_data_to_cell_data()
    out_data = result[scalar]

    # result.save('test.vtk',binary=False)
    
    return out_data
    

# def find_nearest_cellcenter(node_coord,meshIN_nodes_coords,threshold=1e-1, 
#                        **kwargs):
#     '''
#     Find nearest mesh node between two meshes

#     Parameters
#     ----------
#     node_coord : np.array
    
#     meshIN_nodes_coords : np.array

#     threshold : float
#         if distance > threshold --> closest = nan

#     Returns
#     -------
#     closest_idx : list
#         Node indice in the mesh_node.
#     closest : list
#         Node coordinate in the mesh_node.

#     '''
    
#     closest_idx = []
#     closest = []
#     # for i, nc in enumerate(cell_coords):
#         # euclidean distance
        
#     d = ( (meshIN_nodes_coords[:,0] - node_coord[0]) ** 2 + 
#            (meshIN_nodes_coords[:,1]  - node_coord[1]) ** 2 +
#            (meshIN_nodes_coords[:,2]  - node_coord[2]) ** 2
#            ) ** 0.5
    
#     closest_idx.append(np.argmin(d))
#     closest.append(np.vstack(meshIN_nodes_coords[closest_idx,:]))
                
#     if d[np.argmin(d)]>threshold:
#         closest = 'nan'

#     return closest_idx, closest


def find_nearest_node(node_coord,meshIN_nodes_coords,threshold=1e-1, 
                       **kwargs):
    '''
    Find nearest mesh node between two meshes

    Parameters
    ----------
    node_coord : np.array
    
    meshIN_nodes_coords : np.array

    threshold : float
        if distance > threshold --> closest = nan

    Returns
    -------
    closest_idx : list
        Node indice in the mesh_node.
    closest : list
        Node coordinate in the mesh_node.

    '''
    
    closest_idx = []
    closest = []
    # for i, nc in enumerate(cell_coords):
        # euclidean distance
        
    d = ( (meshIN_nodes_coords[:,0] - node_coord[0]) ** 2 + 
           (meshIN_nodes_coords[:,1]  - node_coord[1]) ** 2 +
           (meshIN_nodes_coords[:,2]  - node_coord[2]) ** 2
           ) ** 0.5
    
    closest_idx.append(np.argmin(d))
    closest.append(np.vstack(meshIN_nodes_coords[closest_idx,:]))
                
    if d[np.argmin(d)]>threshold:
        closest = 'nan'

    return closest_idx, closest



def add_attribute_2mesh(data, mesh, name='ER_pred', overwrite=True, 
                        saveMesh=True, **kwargs):
    '''
    add a new mesh attribute to a vtk file

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    mesh : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.
    overwrite : TYPE, optional
        DESCRIPTION. The default is True.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    if type(mesh) is str:
        mesh = pv.read(mesh)
    

    try:
        mesh.point_data[name] = data
    except:
        mesh.cell_data[name] = data

    meshname = name + '.vtk'
    
    if saveMesh:
        path = os.getcwd()
        if 'path' in kwargs:
            path = kwargs['path']
           
        if 'time' in kwargs:
            time = kwargs['time']
            meshname = name  + str(time) +'.vtk' 
            mesh.save(path + meshname, binary=False)
        else:
            mesh.save(path + meshname, binary=False)
                    
        # if overwrite==True:
        #     mesh.save(full_path)
        
    return mesh, name
