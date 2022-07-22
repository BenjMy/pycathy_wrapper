# -*- coding: utf-8 -*-
"""Reader for CATHY input files""" 

import os
import numpy as np
import pandas as pd
from pyCATHY.plotters import cathy_plots as pltCT 

def read_atmbc(filename, grid=[], show=False, **kwargs):
    '''
    read atmbc file 
    
    
    Returns
    -------
    atmbc dataframe.

    '''
    
    # filename = './input/atmbc'
    atmbc_file = open(filename, "r")
    lines = atmbc_file.readlines()
    # lines = [ligne.strip() for ligne in lines] # avoid retour chariot
    atmbc_file.close()

    # read header flags
    Line0_header = [int(word) for word in lines[0].split() if word.isdigit()]
    HSPATM = Line0_header[0]
    IETO = Line0_header[1]


    t = []
    value = []
    tstep_idx = []
    # loop over lines
    for i, ll in enumerate(lines):
        if i>0:
            if 'TIME'.casefold() in ll:
                tstep_idx.append(i)
                splt= ll.split()
                t.append(float(splt[0]))

            # two cases (according to file formatting): 
            # numerical values + flag of value
            # -----------------------------------------------------------------
            # 1/ value (numeric) + 'VALUE'
            elif 'VALUE'.casefold() in ll:
                
                value.append(float(ll.split()[0]))

            # only numerical values
            # 2/ value (numeric)
            # -----------------------------------------------------------------
            else:
                # take care of the additionnal empty lines (retour chariot)
                # in the end of the file
                value.append(float(ll.split()[0]))

    if HSPATM != 0: #homogeneous on all surf mesh nodes
              
        d_atmbc = []
        d_atmbc= np.vstack([t,value])
        cols_atmbc = ['time', 'value']
        df_atmbc =  pd.DataFrame(d_atmbc.T,  columns=cols_atmbc)
        
    else:
        
        d_atmbc = []
        for i in range(len(t)):
           d_atmbc.append(np.column_stack((t[i]*np.ones(int(grid["nnod3"])),
                      value[i]*np.ones(int(grid["nnod3"])),
                      grid['nodes_idxyz'])))
        d_atmbc= np.vstack(d_atmbc)    
        cols_atmbc = ['time', 'value', 'nodeidx', 'x', 'y', 'z']
        df_atmbc =  pd.DataFrame(d_atmbc,  columns=cols_atmbc)
    
    if show==True:
        pltCT.show_atmbc(t,value,IETO=IETO)   
        
        try:
            pltCT.show_atmbc_3d(df_atmbc)
        except:
            print('no vtk file found')



    return df_atmbc, HSPATM, IETO




# def read_parm(filename, **kwargs):
#     '''
#     read parm file 
    
    
#     Returns
#     -------
#     parm dict

#     '''
    
#     # filename = './input/atmbc'
#     parm_file = open(filename, "r")
#     lines = parm_file.readlines()
#     # lines = [ligne.strip() for ligne in lines] # avoid retour chariot
#     parm_file.close()


#     header_fmt_parm = [3, 3, 2, 4, 4, 3, 3, 2, 4, 3, 3, 4, 4, 4, 2, 1, 2]
#     counth = 0

#     # with open(
#     #     os.path.join(self.workdir, self.project_name, self.input_dirname, "parm"),
#     #     "w+",
#     # ) as parmfile:
#     with open(file2write, "w+") as parmfile:
        
            
#     # read header flags
#     # Line0_header = [int(word) for word in lines[0].split() if word.isdigit()]
#     # HSPATM = Line0_header[0]
#     # IETO = Line0_header[1]



#     return dict_parm




def read_grid3d(project_name, **kwargs):

    print("reading grid3d")
    grid3d_file = open(os.path.join(project_name, "output/grid3d"), "r")

    # Lines = grid3d_file.readlines()
    try:
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

    except:
        pass