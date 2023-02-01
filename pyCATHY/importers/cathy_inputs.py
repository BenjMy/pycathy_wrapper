# -*- coding: utf-8 -*-
"""Reader for CATHY input files""" 

import os
import numpy as np
import pandas as pd
from pyCATHY.plotters import cathy_plots as pltCT 
import re

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
            # test_char = re.search('[a-zA-Z]', ll)
            # if test_char is not None:
            if 'TIME'.casefold() in ll.casefold():
                tstep_idx.append(i)
                splt= ll.split()
                t.append(float(splt[0]))

            # two cases (according to file formatting): 
            # numerical values + flag of value
            # -----------------------------------------------------------------
            # 1/ value (numeric) + 'VALUE'
            elif 'VALUE'.casefold() in ll.casefold():
                value.append(float(ll.split()[0]))
    
            elif '\n' in ll:
                # take care of the additionnal empty lines (retour chariot)
                # in the end of the file
                pass
            
            else:
                # only numerical values
                # 2/ value (numeric)
                # -----------------------------------------------------------------
                value.append(float(ll.split()[0]))
    
    if len(value) != len(t):
        raise ValueError('Number of values does not match number of times (check flags TIME, VALUE)')

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

def read_dem_parameters(dem_parametersfile):
    
    header_fmt = [1, 1, 1, 1, 3, 3] # nb of values per lines
    dem_parameters = ["delta_x","delta_y","factor",
                      "dostep",
                      "nzone",
                      "nstr",
                      "n1",
                      "ivert",
                      "isp",
                      "base",
                      "zratio(i),i=1,nstr",
                      ]
    counth = 0
    count_parm = 0
    dem_parameters_dict = {}
    zratio = []
    with open(dem_parametersfile,"r") as f:
        for line_value in f:  # iterate over each line
            if counth<len(header_fmt):
                if header_fmt[counth]==1:
                    dem_parameters_dict[dem_parameters[count_parm]] = float(line_value)
                    count_parm += 1
                elif header_fmt[counth]==3:
                    for j in range(3):
                        dem_parameters_dict[dem_parameters[count_parm+j]]= float(line_value.split()[j])
                    count_parm += 3
            elif counth==len(header_fmt):
                for j in range(len(line_value.split())):
                    zratio.append(float(line_value.split()[j]))
                dem_parameters_dict[dem_parameters[count_parm]]= zratio
                count_parm += len(line_value.split())
            counth += 1
                
        return dem_parameters_dict
                    
                    
def read_soil(soilfile,dem_parm,MAXVEG):
    '''
    Read soil file


    This is an example of 3 layers and  2 zones; The inner reading cycle is by zone, and the outer one by layers, i.e.:

    PERMX_z1_str1 PERMY_z1_str1  PERMZ_z1_str1  ELSTOR_z1_str1 POROS_z1_str1 VGNCELL_z1_str1 VGRMCCELL__z1_str1 VGPSATCELL__z1_str1
    PERMX_z2_str1 PERMY_z2_str1  PERMZ_z2_str1  ELSTOR_z2_str1 POROS_z2_str1 VGNCELL_z2_str1 VGRMCCELL__z2_str1 VGPSATCELL__z2_str1
    PERMX_z1_str2 PERMY_z1_str2  PERMZ_z1_str2  ELSTOR_z1_str2 POROS_z1_str2 VGNCELL_z1_str2 VGRMCCELL__z1_str2 VGPSATCELL__z1_str2
    PERMX_z2_str2 PERMY_z2_str2  PERMZ_z2_str2  ELSTOR_z2_str2 POROS_z2_str2 VGNCELL_z2_str2 VGRMCCELL__z2_str2 VGPSATCELL__z2_str2
    PERMX_z1_str3 PERMY_z1_str3  PERMZ_z1_str3  ELSTOR_z1_str3 POROS_z1_str3 VGNCELL_z1_str3 VGRMCCELL__z1_str3 VGPSATCELL__z1_str3
    PERMX_z2_str3 PERMY_z2_str3  PERMZ_z2_str3  ELSTOR_z2_str3 POROS_z2_str3 VGNCELL_z2_str3 VGRMCCELL__z2_str3  VGPSATCELL__z2_str3
    
            
    Parameters
    ----------
    soilfile : str
        filename (including abs path).
    dem_parm : dict
        dict of dem parameters (from read_dem_parameters)
    Returns
    -------
    dataframe
        df containing soil properties
    dict
        Header of the dem file

    '''
    
    
    # Read FP parameters 
    # --------------------- 
    # To write
    
    # Read SPP parameters 
    # --------------------- 
    soil_header = ['PERMX', 'PERMY','PERMZ',
                   'ELSTOR','POROS',
                   'VGNCELL','VGRMCCELL','VGPSATCELL']
    
    str_hd_soil = {}   
    nb_of_header_lines = 9 + MAXVEG -1 
    with open(
        os.path.join(soilfile), "r"
    ) as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            # if count < 9:
                # str_hd, value_hd = line.split()  # split it by whitespace
                # soilfile[str_hd.replace(":", "")] = value_hd
            count += 1
            
    soil = np.loadtxt(soilfile, skiprows=nb_of_header_lines, 
                         max_rows=count-nb_of_header_lines-1)   
    
    np.shape(soil)
    name_str = []
    name_zone = []
    for dz in range(int(dem_parm['nzone'])):
        for ds in range(int(dem_parm['nstr'])):
            name_str.append(str(ds)) 
            name_zone.append(str(dz))
            
    if len(soil) != len(name_str):
        raise ValueError('Inconsistent number of zones/layers with respect to the number of soil lines')
    df_soil = pd.DataFrame(soil, [name_str,name_zone], soil_header)
    df_soil.index.set_names('str', level=0, inplace=True)
    df_soil.index.set_names('zone', level=1, inplace=True)

    return df_soil, str_hd_soil

def read_raster(filename):
    '''
    Applicable for all DEM type of CATHY outputs such as dtm_Ws1_sf_1, dtm_Ws1_sf_2, ...

    Parameters
    ----------
    filename : str
        filename of the raster type file to read (including abs path).
    Returns
    -------
    np.array([])
        2d array containing the raster info.
    dict
        Header of the raster file

    '''
    str_hd_raster = {}       
    with open(
        os.path.join(filename), "r"
    ) as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_raster[str_hd.replace(":", "")] = value_hd
            count += 1
    raster_mat = np.loadtxt(filename, skiprows=6)
    return raster_mat, str_hd_raster
    
    
    
def read_zone(zonefile):
    '''
    Read zone file and infer soil zones in the DEM

    Parameters
    ----------
    zonefile : str
        filename (including abs path).
    Returns
    -------
    np.array([])
        2d array containing the zones.
    dict
        Header of the zone file

    '''
    str_hd_zones = {}       
    with open(
        os.path.join(zonefile), "r"
    ) as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_zones[str_hd.replace(":", "")] = value_hd
            count += 1
    zones_mat = np.loadtxt(zonefile, skiprows=6)
    return zones_mat, str_hd_zones

def read_dem(demfile,dtm_13file):
    '''
    Read dem file and infer number of altitudes in the DEM

    Parameters
    ----------
    demfile : str
        filename (including abs path).
    dtm_13file : str
        filename (including abs path).
    Returns
    -------
    np.array([])
        2d array containing the elevation.
    dict
        Header of the dem file

    '''
    str_hd_dem = {}       
    with open(
        os.path.join(demfile), "r"
    ) as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_dem[str_hd.replace(":", "")] = value_hd
            count += 1
    dem_file = open(dtm_13file, "r")
    dem_mat = np.loadtxt(dem_file, skiprows=0)
    dem_file.close()
    
    return dem_mat, str_hd_dem


def read_root_map(rootmapfile):
    '''
    Read root map file and infer number of vegetations in the DEM

    Parameters
    ----------
    rootmapfile : str
        filename (including abs path).

    Returns
    -------
    np.array([])
        2d array containing the flag for the vegetation type.
    dict
        Header of the root map file

    '''
    
    str_hd_rootmap = {}       
    with open(
        os.path.join(rootmapfile), "r"
    ) as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_rootmap[str_hd.replace(":", "")] = value_hd
            count += 1
    rootmap_mat = np.loadtxt(rootmapfile, skiprows=6)
    
    return rootmap_mat, str_hd_rootmap
        
def read_grid3d(project_name, **kwargs):

    # print("reading grid3d")
    grid3d_file = open(os.path.join(project_name, "output/grid3d"), "r")

    # Lines = grid3d_file.readlines()
    try:
        nnod, nnod3, nel = np.loadtxt(grid3d_file, max_rows=1)

        grid3d_file.close()
    
        grid3d_file = open(os.path.join(project_name, "output/grid3d"), "r")
        mesh_tetra = np.loadtxt(grid3d_file, skiprows=1, max_rows=int(nel))
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