#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshing tools
"""

import os 
import numpy as np
import pyvista as pv

try: 
    import pygimli as pg
    import pygimli.meshtools as mt
except ImportError: 
    pygimli = None

import matplotlib.pylab as plt


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
    
    in_nodes_mod = np.array(meshIN.points)
    if 'in_nodes_mod' in kwargs:
        in_nodes_mod = kwargs['in_nodes_mod']    
    meshIN.set_active_scalars(scalar)
    meshIN.points = in_nodes_mod
    
    # set_interpolation_radius()
    
    rd = max(np.diff(meshIN.points[:,0]))/1
    result = meshOUT.interpolate(meshIN, radius=rd, pass_point_data=True)
    
    # plot_2d_interpolation_quality(meshIN,scalar,meshOUT,result)

    result = result.point_data_to_cell_data()
    out_data = result[scalar]
    
    warm_0 = ''
    if len(np.where(out_data == 0))>0:
        warm_0 = 'interpolation created 0 values - replacing them by min value of input CATHY predicted ER mesh'
        
    out_data = np.where(out_data == 0, 1e-3, out_data)
    
    
    return out_data, warm_0

    
def set_interpolation_radius():
    
    # rd= min([abs(min(np.diff(meshIN.points[:,0]))),
    #      abs(min(np.diff(meshIN.points[:,1]))),
    #      abs(min(np.diff(meshIN.points[:,2])))
    #      ]
    #     )
    
    pass
    
    
def plot_2d_interpolation_quality(meshIN,scalar,meshOUT,result):
    
    fig = plt.figure()
    ax1 = plt.subplot(131)
    print(max(meshIN[scalar]))
    print(min(meshIN[scalar]))
    
    meshIN.points[:,0].min()

    meshOUT.points[:,0].min()
    meshOUT.points[:,0].max()
    meshOUT.points[:,1].min()
    meshOUT.points[:,1].max()
    
    cm = plt.cm.get_cmap('RdYlBu')
    # result = meshOUT.interpolate(meshIN, radius=rd, pass_point_data=True)
    sc = ax1.scatter(meshIN.points[:,0],meshIN.points[:,1],c=meshIN[scalar],label='meshIN[scalar]',
                s=35, cmap=cm)
    plt.colorbar(sc)
    cm = plt.cm.get_cmap('RdYlBu')
    # fig = plt.figure()
    ax2 = plt.subplot(132)
    sc = ax2.scatter(meshOUT.points[:,0],meshOUT.points[:,1],c=result[scalar],label='result[scalar]',
                s=35, cmap=cm,vmin=min(meshIN[scalar]), vmax=max(meshIN[scalar]))
    ax2.set_xlim([min(meshIN.points[:,0]),max(meshIN.points[:,0])])
    ax2.set_ylim([min(meshIN.points[:,1]),max(meshIN.points[:,1])])
    
    # fig = plt.figure()
    ax3 = plt.subplot(133)
    sc = ax3.scatter(meshOUT.points[:,0],meshOUT.points[:,1],c=result[scalar],label='result[scalar]',
                s=35, cmap=cm)
    plt.colorbar(sc)
    plt.tight_layout()

    
    def uniquify(path):
        filename, extension = os.path.splitext(path)
        counter = 1
    
        while os.path.exists(path):
            path = filename + str(counter) + extension
            counter += 1

        return path


    savedir = os.getcwd()
    savename_test = os.path.join(savedir, 'interpolation_q.png')
    savename = uniquify(savename_test)
    print(savename)

    plt.savefig(savename)
    # plt.show()
    
    pass

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
    
    # for k in kwargs:
    #     print(k)
    # print(mesh)
    
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
