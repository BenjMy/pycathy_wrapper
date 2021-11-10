#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 18:03:04 2021

@author: ben
"""

import os 
import numpy as np
import pyvista as pv

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
    
    # print('mean(meshIN.get_array(scalar))')
    # print(np.mean(meshIN.get_array(scalar)))

    in_nodes_mod = np.array(meshIN.points)
    if 'in_nodes_mod' in kwargs:
        in_nodes_mod = kwargs['in_nodes_mod']

    out_data = []
    closest_idx = []
    # for pos in meshOUT.points:
        
        # closest_idx, closest = find_nearest_nodes(pos,in_nodes_mod,
        #                                              threshold)

    if type(meshOUT) is str:
        # print('read vtk file using pv')
        meshOUT = pv.read(meshOUT)
    # print('lenScalar' + str(len(meshIN.get_array(scalar))))
    # print(meshOUT)

    cellOUT_centers = meshOUT.cell_centers()

    for pos in cellOUT_centers.points:

        closest_idx, closest = find_nearest_node(pos,in_nodes_mod,
                                                     threshold)

        if 'nan' in closest:
            out_data.append('nan')
        else:
            out_data.append(meshIN.get_array(scalar)[closest_idx])
    
    out_data = np.hstack(out_data)

    
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



def add_attribute_2mesh(data, mesh, name='ER_pred', overwrite=True, **kwargs):
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
        
        
    mesh.add_field_data(data, name)
    
    meshname = name + '.vtk'


    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']

    print('meshpath=' + path)
    
    if 'time' in kwargs:
        time = kwargs['time']
        meshname = name  + str(time) +'.vtk'

        mesh.save(path + meshname)
    else:
        mesh.save(path + name + '.vtk')
                
    # if overwrite==True:
    #     mesh.save(full_path)

    return mesh




def mshParse(self, gmsh_mesh=None):
    '''
    Forked from Resipy 'meshtool' lib
    

    import and parse gmsh mesh
    Parameters
    ----------
    gmsh_mesh : TYPE, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.
    ValueError
        DESCRIPTION.

    Returns
    -------
    mesh_dict : TYPE
        DESCRIPTION.
    con_matrix : TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''

    if not isinstance(gmsh_mesh, str):
        raise Exception("expected a string argument for msh_parser")

    fid = open(gmsh_mesh, "r")  # Open text file
    # Idea: Read Mesh format lines $MeshFormat until $Nodes
    dump = fid.readlines()
    fid.close()
    if len(dump) <= 1:
        raise Exception("Target file is empty!!...aborting!")
    # check the file is a mesh format
    if (
        dump[0].strip() != "$MeshFormat"
    ):  # removes formating strings, checks if the file is a gmsh file
        raise Exception("Unrecognised file type...aborting!")
    mesh_format = dump[1].strip()
    formats = ["2.2 0 8", "4 0 8"]
    if (
        mesh_format not in formats
    ):  # warn people that the code was developed with a different mesh format in mind
        warnings.warn("Mesh file format unrecognised ... some errors may occur!\n")

    if mesh_format == "2.2 0 8":
        gmshV = 3  # assume its gmsh version 3.06
    else:
        gmshV = 4  # assume V4 and above

    # read gmsh .msh mesh # format 2

    # find where the nodes start
    for i, line in enumerate(dump):
        if line.find("$Nodes") != -1:  # node flag
            line = dump[i + 1].split()
            no_nodes = int(line[-1])
            node_start = i + 2
        elif line.find("$EndNodes") != -1:
            node_end = i
            break  # stop the loop, should find the nodes start before the end

    # read in number of nodes - at line 5
    # allocate lists for node numbers and coordinates
    node_num = [0] * no_nodes
    nodex = [0] * no_nodes
    nodey = [0] * no_nodes
    nodez = [0] * no_nodes
    # read in node information
    gmshV = 3
    if gmshV == 3:
        for i in range(node_start, node_end):
            line_info = dump[i].split()
            # convert string info into floats
            line = [float(k) for k in line_info]
            node_idx = int(line[0])
            node_num[node_idx - 1] = node_idx
            nodex[node_idx - 1] = line[1]
            nodey[node_idx - 1] = line[2]
            nodez[node_idx - 1] = line[3]
    else:
        c = 0
        i = node_start
        while c < no_nodes:
            line_info = dump[i].split()
            line = [int(k) for k in line_info]  # entity line
            ent = line[-1]  # number of nodes to follow
            ent_start = i  # starting line
            for j in range(ent):
                # print(i+j)
                line_info = dump[ent_start + 1 + j].split()
                # convert string info into floats
                line = [float(k) for k in line_info]
                node_idx = int(line[0])
                node_num[node_idx - 1] = node_idx
                nodex[node_idx - 1] = line[1]
                nodey[node_idx - 1] = line[2]
                nodez[node_idx - 1] = line[3]
                c += 1
                i += 1
            i += 1

    #### read in elements
    # find where the elements start
    for i, line in enumerate(dump):
        if line.find("$Elements") != -1:
            element_start = i + 2
        if line.find("$EndElements") != -1:
            element_end = i
            break  # stop the loop, should find the nodes start before the end

    # engage for loop - this time we want to filter out elements which are not needed
    nat_elm_num = []  # native element number to gmsh
    elm_type = []  # element type
    phys_entity = []  # defines the physical entity type the element is assocaited with

    ignored_elements = 0  # count the number of ignored elements

    # determine element type
    for i in range(element_start, element_end):
        line = dump[i].split()
        print(line)
        if gmshV == 3:
            elm_type.append(int(line[1]))
        else:
            elm_type.append(len(line) - 1)

    if gmshV == 3:
        lookin4 = [2, 4, 6]  # number of cell nodes we are looking for
    else:
        lookin4 = [3, 4, 6]
    if lookin4[0] in elm_type:  # then its triangles
        npere = 3  # number of nodes per element
    if lookin4[1] in elm_type:  # forget triangles its tetrahedra
        npere = 4
    if lookin4[2] in elm_type:  # forget triangles its prisms
        npere = 6

    if npere == 3:
        # stream('Triangle')
        con_matrix = [[], [], []]
        cell_type = [5]
    elif npere == 4:
        # stream('Tetrahedra')
        con_matrix = [[], [], [], []]
        cell_type = [10]
    elif npere == 6:
        # stream('Prism')
        con_matrix = [[], [], [], [], [], []]
        cell_type = [13]
    else:
        raise ValueError(
            "Cannot parse mesh becuase the relevant cell types cannot be found"
        )

    phys = 0
    if gmshV == 3:  # parse elements for older gmsh file type
        for i in range(element_start, element_end):
            splitted = dump[i].split()
            line = [int(k) for k in splitted]
            elmV = line[1]
            if npere == 3:
                elmV += 1  # in this format the flag for triangle elements is 2 so add 1
            # convert string info into ints and cache data
            if npere == elmV:  # if line contains number of element vertices we want
                nat_elm_num.append(line[0])
                elm_type.append(line[1])
                phys_entity.append(line[4])
                for j in range(npere):
                    con_matrix[j].append(line[5 + j] - 1)
            else:
                ignored_elements += 1
    else:  # newer gmsh parser
        c = 0
        i = element_start
        while i < element_end:
            # read in flag line
            splitted = dump[i].split()
            line = [int(k) for k in splitted]
            nef = line[3]  # number of elements to follow
            elmV = line[2]  # number of element vertices
            phys = line[0]  # physical entity
            if npere == 3:
                elmV += 1  # in this format the flag for triangle elements is 2 so add 1
            if npere == elmV:
                nef_start = i + 1
                for j in range(nef):
                    splitted = dump[nef_start + j].split()
                    line = [int(k) for k in splitted]
                    nat_elm_num.append(line[0])
                    # elm_type.append(line[1])
                    phys_entity.append(phys)
                    for k in range(npere):
                        con_matrix[k].append(line[1 + k] - 1)
                    c += 1
                    i += 1
                i += 1
            else:
                ignored_elements += nef
                i += nef + 1

    real_no_elements = len(
        nat_elm_num
    )  #'real' number of elements that we actaully want
    if real_no_elements == 0:
        raise Exception(
            "No elements have been imported, please check formatting of .msh file"
        )

    elm_id = [i + 1 for i in range(real_no_elements)]

    mesh_dict = {
        "num_elms": real_no_elements,
        "num_nodes": no_nodes,
        "dump": dump,
        "node_x": nodex,  # x coordinates of nodes
        "node_y": nodey,  # y coordinates of nodes
        "node_z": nodez,  # z coordinates of nodes
        "node_id": node_num,  # node id number
        "elm_id": elm_id,  # element id number
        "parameters": phys_entity,  # the values of the attributes given to each cell
        "parameter_title": "regions",
        "dict_type": "mesh_info",
    }

    return mesh_dict, con_matrix, phys_entity  # return a python dictionary
