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


def CATHY_2_Simpeg(mesh_CATHY,ERT_meta_dict,scalar='saturation',show=False,**kwargs):
    pass

def CATHY_2_pg(mesh_CATHY,ERT_meta_dict,scalar='saturation',show=False,**kwargs):
    '''   
    Interpolate CATHY mesh attribute to pygimli mesh.
    Add a new [`scalar`] attribute to the pygimli mesh (create a new mesh)
    
    .. Note:
        Need to flip axis because convention for CATHY and pygimli are different

    Parameters
    ----------
    mesh_CATHY : pvmesh
        CATHY mesh to transform to pygimli.
    ERT_meta_dict :dict
        Dictionnary containing ERT metadata (mesh, format, ..).
    scalar : str, optional
        scalar attribute to interpolate. The default is 'saturation'.
    show : bool, optional
        show the result of the interpolation using pyvista. The default is False.
    **kwargs : TYPE
        path : path of the mesh to overwrite 

    Returns
    -------
    mesh_new_attr : TYPE
        DESCRIPTION.
    scalar_new : TYPE
        DESCRIPTION.

    '''
            
    
    if type(ERT_meta_dict['forward_mesh_vtk_file']) is str:
        mesh_OUT = pv.read(ERT_meta_dict['forward_mesh_vtk_file'])
    
    # flip y and z axis as CATHY and pg have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    in_nodes_mod_pg = np.array(mesh_OUT.points)
    
    
    idx = np.array([0, 2, 1]) # reorder xyz to xzy
    in_nodes_mod_m = in_nodes_mod[:, idx]
    # in_nodes_mod_m = in_nodes_mod[:, idx]
    in_nodes_mod_m[:,2] = -np.flipud(in_nodes_mod_m[:,2])
    in_nodes_mod_m[:,1] = -np.flipud(in_nodes_mod_m[:,1])

    # for i, axesi in enumerate(['x','y','z']):
    #     print('check {} consistency between the 2 meshes'.format(axesi))
    #     print(max(in_nodes_mod_m[:,i]),max(in_nodes_mod_pg[:,i]))
    #     print(min(in_nodes_mod_m[:,i]),min(in_nodes_mod_pg[:,i]))


    if 'mesh_nodes_modif' in ERT_meta_dict.keys():
        in_nodes_mod_m = ERT_meta_dict['mesh_nodes_modif']
        
           
    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']


        
        
    # mesh_OUT, meshOUT = trace_mesh(mesh_CATHY,mesh_OUT,
    #                     scalar=scalar,
    #                     threshold=1e-1,
    #                     in_nodes_mod=in_nodes_mod_m
    #                     )
    
    data_OUT, warm_0 = trace_mesh(mesh_CATHY,mesh_OUT,
                                  scalar=scalar,
                                  threshold=1e-1,
                                  in_nodes_mod=in_nodes_mod_m
                                  )
           
    if len(warm_0)>0:
        print(warm_0)
        
    
    scalar_new = scalar + '_nearIntrp2_pg_msh' 
    if 'time' in kwargs:
        time = kwargs['time']
        mesh_new_attr, name_new_attr = add_attribute_2mesh(data_OUT, 
                                                            mesh_OUT, 
                                                            scalar_new, 
                                                            overwrite=True,
                                                            time=time,
                                                            path=path)
    else:
        mesh_new_attr, name_new_attr = add_attribute_2mesh(data_OUT, 
                                                            mesh_OUT, 
                                                            scalar_new, 
                                                            overwrite=True,
                                                            path=path)
    
    if show:

        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=False)
        p.add_mesh(mesh_new_attr,scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
        
        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
            
    
    return mesh_new_attr, scalar_new


def CATHY_2_Resipy(mesh_CATHY,mesh_Resipy,scalar='saturation',show=False,**kwargs):
 

    # flip y and z axis as CATHY and Resipy have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    in_nodes_mod[:,2] = -np.flipud(in_nodes_mod[:,2])
    in_nodes_mod[:,1] = -np.flipud(in_nodes_mod[:,1])


    # check with a plot positon of the nodes for both meshes
    # ------------------------------------------------------------------------
    # p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
    # p.add_mesh(mesh_CATHY)
    # _ = p.add_points(np.array(mesh_CATHY.points), render_points_as_spheres=True,
    #                         color='red', point_size=20)
    # _ = p.add_points(in_nodes_mod, render_points_as_spheres=True,
    #                         color='blue', point_size=20)
    # _ = p.show_bounds(grid='front', all_edges=True, font_size=50)
    # cpos = p.show(True)
    
    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']

    data_OUT = trace_mesh(mesh_CATHY,mesh_Resipy,
                            scalar=scalar,
                            threshold=1e-1,
                            in_nodes_mod=in_nodes_mod)
    
    # print(len(data_OUT))
    scalar_new = scalar + '_nearIntrp2Resipymsh'
    # print(mesh_Resipy)
    # get_array(mesh, name, preference='cell'
              

    if 'time' in kwargs:
        time = kwargs['time']
        mesh_new_attr, name_new_attr = add_attribute_2mesh(data_OUT, 
                                                                mesh_Resipy, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                time=time,
                                                                path=path)
    else:
        mesh_new_attr, name_new_attr = add_attribute_2mesh(data_OUT, 
                                                                mesh_Resipy, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                path=path)
    
    if show == True:

        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_new_attr,scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
    # if type(meshERT) is str:
    #     meshERTpv = pv.read(meshERT)
    
    # if savefig == True:
    
    #     plotter = pv.Plotter(notebook=True)
    #     _ = plotter.add_mesh(mesh_new_attr,show_edges=True)
    #     plotter.view_xz(negative=False)
    #     plotter.show_grid()
    #     plotter.save_graphic(path_CATHY + 'ERT' + str(DA_cnb) + str('.ps'), 
    #                           title='ERT'+ str(DA_cnb), 
    #                           raster=True, 
    #                           painter=True)
            
    
    return mesh_new_attr, scalar_new

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
    meshOUT_interp = meshOUT.interpolate(meshIN, radius=rd, pass_point_data=True)
    
    # plot_2d_interpolation_quality(meshIN,scalar,meshOUT,meshOUT_interp)
    
    
    result = meshOUT_interp.point_data_to_cell_data()
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
    
    # fig = plt.figure()
    # ax1 = plt.subplot(131)
    # # print(max(meshIN[scalar]))
    # # print(min(meshIN[scalar]))
    
    # # meshIN.points[:,0].min()
    # # meshOUT.points[:,0].min()
    # # meshOUT.points[:,0].max()
    # # meshOUT.points[:,1].min()
    # # meshOUT.points[:,1].max()
    
    # cm = plt.cm.get_cmap('RdYlBu')
    # # result = meshOUT.interpolate(meshIN, radius=rd, pass_point_data=True)
    # sc = ax1.scatter(meshIN.points[:,0],meshIN.points[:,1],c=meshIN[scalar],label='meshIN[scalar]',
    #             s=35, cmap=cm)
    # plt.colorbar(sc)
    cm = plt.cm.get_cmap('RdYlBu')
    # # fig = plt.figure()
    
    fig, axs = plt.subplots(2)


    # ax2 = plt.subplot(121)
    sc = axs[0].scatter(meshIN.points[:,0],meshIN.points[:,1],c=meshIN[scalar],label='meshIN[scalar]',
                s=35, cmap=cm)
    #,vmin=min(meshIN[scalar]), vmax=max(meshIN[scalar])
    axs[0].set_xlim([min(meshIN.points[:,0]),max(meshIN.points[:,0])])
    axs[0].set_ylim([min(meshIN.points[:,1]),max(meshIN.points[:,1])])
    
    # fig = plt.figure()
    # ax3 = plt.subplot(122)
    sc = axs[1].scatter(meshOUT.points[:,0],meshOUT.points[:,1],c=result[scalar],label='result[scalar]',
                s=35, cmap=cm)

    axs[1].set_xlim([min(meshIN.points[:,0]),max(meshIN.points[:,0])])
    axs[1].set_ylim([min(meshIN.points[:,1]),max(meshIN.points[:,1])])
    
    import matplotlib.ticker as ticker
    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    plt.colorbar(sc,format=ticker.FuncFormatter(fmt))
    # plt.tight_layout()

    
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
    # print(savename)

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





def trace_mesh_pg(meshIN,meshOUT,method='spline', **kwargs):
    '''
    Interpolate CATHY mesh (structured) into pygimli mesh (structured) using pygimli meshtools
    # Specify interpolation method 'linear, 'spline', 'harmonic'

    '''
    meshIN = pg.load(meshIN)
    meshOUT = pg.load(meshOUT)
    out_data = pg.interpolate(meshIN['ER_converted_CATHYmsh*'], meshOUT, method=method)
    
    return out_data
    