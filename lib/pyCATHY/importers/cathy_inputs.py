"""Reader for CATHY input files
""" 

import os
import numpy as np

#%% Meshtool functions



def read_atmbc(filename):
    '''

    Returns
    -------
    None.

    '''
    
    atmbc_file = open(filename, "r")

    
    atmbcfile.close()
    
    

    
    pass   




def read_grid3d(project_name, **kwargs):

    print("reading grid3d")
    grid3d_file = open(os.path.join(project_name, "output/grid3d"), "r")

    # Lines = grid3d_file.readlines()
    nnod, nnod3, nel = np.loadtxt(grid3d_file, max_rows=1)
    grid3d_file.close()

    grid3d_file = open(os.path.join(project_name, "output/grid3d"), "r")
    mesh_tetra = np.loadtxt(grid3d_file, skiprows=1, max_rows=int(nel) - 1)
    grid3d_file.close()

    grid3d_file = open(os.path.join(project_name, "output/grid3d"), "r")
    mesh3d_nodes = np.loadtxt(
        grid3d_file,
        skiprows=1 + int(nel),
        max_rows=1 + int(nel) + int(nnod3) - 1,
    )
    grid3d_file.close()

    xyz_file = open(os.path.join(project_name, "output/xyz"), "r")
    nodes_idxyz = np.loadtxt(xyz_file, skiprows=1)
    xyz_file.close()

    # self.xmesh = mesh_tetra[:,0]
    # self.ymesh = mesh_tetra[:,1]
    # self.zmesh = mesh_tetra[:,2]
    # return mesh3d_nodes

    grid = {
        "nnod": nnod,  # number of surface nodes
        "nnod3": nnod3,  # number of volume nodes
        "nel": nel,
        "mesh3d_nodes": mesh3d_nodes,
        "mesh_tetra": mesh_tetra,
        "nodes_idxyz": nodes_idxyz,
    }

    return grid