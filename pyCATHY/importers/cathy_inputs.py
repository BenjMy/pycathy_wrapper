# -*- coding: utf-8 -*-
"""Reader for CATHY input files"""

import os
import re

import numpy as np
import pandas as pd

from pyCATHY.plotters import cathy_plots as pltCT


def read_atmbc(filename, grid=[], show=False, **kwargs):
    """
    read atmbc file


    Returns
    -------
    atmbc dataframe.

    """

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
        if i > 0:
            # test_char = re.search('[a-zA-Z]', ll)
            # if test_char is not None:
            if "TIME".casefold() in ll.casefold():
                tstep_idx.append(i)
                splt = ll.split()
                t.append(float(splt[0]))

            # two cases (according to file formatting):
            # numerical values + flag of value
            # -----------------------------------------------------------------
            # 1/ value (numeric) + 'VALUE'
            elif "VALUE".casefold() in ll.casefold() or "ATMINP".casefold() in ll.casefold():
                value.append(float(ll.split()[0]))

            elif "\n" in ll:
                # take care of the additionnal empty lines (retour chariot)
                # in the end of the file
                pass

            else:
                # only numerical values
                # 2/ value (numeric)
                # -----------------------------------------------------------------
                value.append(float(ll.split()[0]))

    if len(value) != len(t):
        raise ValueError(
            "Number of values does not match number of times (check flags TIME, VALUE)"
        )

    if HSPATM != 0:  # homogeneous on all surf mesh nodes
        d_atmbc = []
        d_atmbc = np.vstack([t, value])
        cols_atmbc = ["time", "value"]
        df_atmbc = pd.DataFrame(d_atmbc.T, columns=cols_atmbc)

    else:
        d_atmbc = []
        for i in range(len(t)):
            d_atmbc.append(
                np.column_stack(
                    (
                        t[i] * np.ones(int(grid["nnod3"])),
                        value[i] * np.ones(int(grid["nnod3"])),
                        grid["nodes_idxyz"],
                    )
                )
            )
        d_atmbc = np.vstack(d_atmbc)
        cols_atmbc = ["time", "value", "nodeidx", "x", "y", "z"]
        df_atmbc = pd.DataFrame(d_atmbc, columns=cols_atmbc)

    if show == True:
        pltCT.show_atmbc(t, value, IETO=IETO)

        try:
            pltCT.show_atmbc_3d(df_atmbc)
        except:
            print("no vtk file found")

    return df_atmbc, HSPATM, IETO


def read_parm(filename, **kwargs):
    '''
    read parm file


    Returns
    -------
    parm dict

    '''


    # parm_header = [
    #     "PERMX",
    #     "PERMY",
    #     "PERMZ",
    #     "ELSTOR",
    #     "POROS",
    #     "VGNCELL",
    #     "VGRMCCELL",
    #     "VGPSATCELL",
    # ]
    
    parm_header = [
        "IPRT1",  # Flag for output of input and coordinate data
        "NCOUT",
        "TRAFLAG",  # Flag for transport
        
        "ISIMGR",  # Flag for type of simulation and type of surface grid
        "PONDH_MIN",  # Minimum ponding head
        "VELREC",
        
        "KSLOPE",
        "TOLKSL",
        
        "PKRL",
        "PKRR",
        "PSEL",
        "PSER",
        
        "PDSE1L",
        "PDSE1R",
        "PDSE2L",
        "PDSE2R",
        
        "ISFONE",
        "ISFCVG",
        "DUPUIT",
        
        "TETAF",
        "LUMP",
        "IOPT",
        
        "NLRELX",
        "OMEGA",
        
        "L2NORM",
        "TOLUNS",
        "TOLSWI",
        "ERNLMX",
        
        "ITUNS",
        "ITUNS1",
        "ITUNS2",
        
        "ISOLV",
        "ITMXCG",
        "TOLCG",
        
        "DELTAT",
        "DTMIN",  # Minimum FLOW3D time step size allowed
        "DTMAX",  # Maximum FLOW3D time step size allowed
        "TMAX", # Time at end of simulation (TMAX is set to 0.0 for steady state problem)

        "DTMAGA",
        "DTMAGM",
        "DTREDS",
        "DTREDM",
        
        "IPRT",
        "VTKF",
        "NPRT",
        "(TIMPRT(I),I=1,NPRT)",
        
        "NUMVP",
        "(NODVP(I),I=1,NUMVP)",  # should be positive (first node is 1)
        
        "NR",
        "CONTR(I),I=1,NR",
        "NUM_QOUT",
        "(ID_QOUT(I),I=1,NUM_QOUT)",
        
    ]

    dict_parm = {}
    parm_file = open(filename, "r")
    lines = parm_file.readlines()
    lines = [ligne.strip() for ligne in lines] # avoid retour chariot
    parm_file.close()
    lines = [ll for ll in lines if len(ll)>0]
    
    lines_s = [ll.split() for ll in lines]
    flat_lines_s = [item for sublist in lines_s for item in sublist]
    lines_s_num, lines_s_names = _search_num_values_in_list(flat_lines_s)
   
    ii = 0
    for lname in lines_s_names:
        if lname == "(TIMPRT(I),I=1,NPRT)" or lname == "(TIMPRT(I).I=1.NPRT)": 
            lines_s_num_TIMPRT = []
            for k in range(dict_parm['NPRT']):
                lines_s_num_TIMPRT.append(lines_s_num[ii+k])
            dict_parm[lname] = lines_s_num_TIMPRT
            ii += dict_parm['NPRT']
        elif lname == "(NODVP(I),I=1,NUMVP)" or lname== "(NODVP(I).I=1.NUMVP)":
            lines_s_num_NODVP = []
            for k in range(dict_parm['NUMVP']):
                lines_s_num_NODVP.append(lines_s_num[ii+k])
            dict_parm[lname] = lines_s_num_NODVP
            ii += dict_parm['NUMVP']
        elif lname == "(ID_QOUT(I).I=1.NUM_QOUT)" or lname == "(ID_QOUT(I),I=1,NUM_QOUT)": 
            lines_s_num_ID_QOUT = []
            if dict_parm['NUM_QOUT'] == 0:
                dict_parm[lname] = 441
                ii +=1
            else:                
                for k in range(dict_parm['NUM_QOUT']):
                    lines_s_num_ID_QOUT.append(lines_s_num[ii+k])
                dict_parm[lname] = lines_s_num_ID_QOUT
                ii += dict_parm['NUM_QOUT']
        else:
            dict_parm[lname] = lines_s_num[ii]
            ii +=1

    return dict_parm


def read_dem_parameters(dem_parametersfile):

    # header_fmt = [1, 1, 1, 1, 3, 3]  # nb of values per lines
    dem_parameters = [
        "delta_x",
        "delta_y",
        "factor",
        "dostep",
        "nzone",
        "nstr",
        "n1",
        "ivert",
        "isp",
        "base",
        "zratio(i),i=1,nstr",
    ]
    
    
    parm_file = open(dem_parametersfile, "r")
    lines = parm_file.readlines()
    lines = [ligne.strip() for ligne in lines] # avoid retour chariot
    parm_file.close()
    lines = [ll for ll in lines if len(ll)>0]
    
    lines_s = [ll.split() for ll in lines]
    flat_lines_s = [item for sublist in lines_s for item in sublist]
    lines_s_num, lines_s_names = _search_num_values_in_list(flat_lines_s)
    

    dict_dem_parm = {}
    ii = 0
    for lname in lines_s_names:
        if lname == "zratio(i),i=1,nstr" or lname == 'zratio(i).i=1.nstr': 
            lines_s_num_zratio = []
            for k in range(dict_dem_parm['nstr']):
                # print(k)
                lines_s_num_zratio.append(lines_s_num[ii+k])
            dict_dem_parm[lname] = lines_s_num_zratio
            ii += dict_dem_parm['nstr']
        else:
            dict_dem_parm[lname] = lines_s_num[ii]
            ii +=1
            
            
    # counth = 0
    # count_parm = 0
    # dem_parameters_dict = {}
    # zratio = []
    # with open(dem_parametersfile, "r") as f:
    #     for line_value in f:  # iterate over each line
    #         if counth < len(header_fmt):
    #             if header_fmt[counth] == 1:
    #                 dem_parameters_dict[dem_parameters[count_parm]] = float(line_value)
    #                 count_parm += 1
    #             elif header_fmt[counth] == 3:
    #                 for j in range(3):
    #                     dem_parameters_dict[dem_parameters[count_parm + j]] = float(
    #                         line_value.split()[j]
    #                     )
    #                 count_parm += 3
    #         elif counth == len(header_fmt):
    #             for j in range(len(line_value.split())):
    #                 zratio.append(float(line_value.split()[j]))
    #             dem_parameters_dict[dem_parameters[count_parm]] = zratio
    #             count_parm += len(line_value.split())
    #         counth += 1

    return dict_dem_parm


def read_soil(soilfile, dem_parm, MAXVEG):
    """
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

    """

    # Read FP parameters
    # ---------------------
    # To write

    # Read SPP parameters
    # ---------------------
    soil_header = [
        "PERMX",
        "PERMY",
        "PERMZ",
        "ELSTOR",
        "POROS",
        "VGNCELL",
        "VGRMCCELL",
        "VGPSATCELL",
    ]

    str_hd_soil = {}
    nb_of_header_lines = 9 + MAXVEG - 1
    with open(os.path.join(soilfile), "r") as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            # if count < 9:
            # str_hd, value_hd = line.split()  # split it by whitespace
            # soilfile[str_hd.replace(":", "")] = value_hd
            count += 1

    soil = np.loadtxt(
        soilfile, skiprows=nb_of_header_lines, max_rows=count - nb_of_header_lines - 1
    )

    np.shape(soil)
    name_str = []
    name_zone = []
    for dz in range(int(dem_parm["nzone"])):
        for ds in range(int(dem_parm["nstr"])):
            name_str.append(str(ds))
            name_zone.append(str(dz))

    if len(soil) != len(name_str):
        raise ValueError(
            "Inconsistent number of zones/layers with respect to the number of soil lines"
        )
    df_soil = pd.DataFrame(soil, [name_str, name_zone], soil_header)
    df_soil.index.set_names("str", level=0, inplace=True)
    df_soil.index.set_names("zone", level=1, inplace=True)

    return df_soil, str_hd_soil


def read_raster(filename):
    """
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

    """
    str_hd_raster = {}
    with open(os.path.join(filename), "r") as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_raster[str_hd.replace(":", "")] = value_hd
            count += 1
    raster_mat = np.loadtxt(filename, skiprows=6)
    return raster_mat, str_hd_raster


def read_zone(zonefile):
    """
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

    """
    str_hd_zones = {}
    with open(os.path.join(zonefile), "r") as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_zones[str_hd.replace(":", "")] = value_hd
            count += 1
    zones_mat = np.loadtxt(zonefile, skiprows=6)
    return zones_mat, str_hd_zones


def read_dem(demfile, dtm_13file):
    """
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

    """
    str_hd_dem = {}
    with open(os.path.join(demfile), "r") as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_dem[str_hd.replace(":", "")] = value_hd
            count += 1
    dem_file = open(dtm_13file, "r")
    dem_mat = np.loadtxt(dem_file, skiprows=0)
    dem_file.close()

    # import matplotlib.pyplot as plt
    # plt.imshow(dem_mat)

    return dem_mat, str_hd_dem


def read_root_map(rootmapfile):
    """
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

    """

    str_hd_rootmap = {}
    with open(os.path.join(rootmapfile), "r") as f:  # open the file for reading
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

def _search_num_values_in_list(flat_lines_s):
    
    lines_s_num = []
    lines_s_names = []
    for lli in flat_lines_s:
        if lli.isdigit():
            lines_s_num.append(int(lli))
        elif lli.replace('-', '', 1).isdigit():
            lines_s_num.append(int(lli))
        elif lli.replace('-', '', 1).replace('.', '', 1).isdigit():
            lines_s_num.append(float(lli))
        elif lli.replace('.', '', 1).isdigit():
            lines_s_num.append(float(lli))
        elif lli.replace('.', '', 1).replace('+', '', 1).isdigit():
            lines_s_num.append(float(lli.replace('+', '')))
        elif lli.replace('.', '', 1).replace('e+', '', 1).isdigit():
            lines_s_num.append(float(lli))
        elif lli.replace('.', '', 1).replace('e-', '', 1).isdigit():
            lines_s_num.append(float(lli))
        else:
            lines_s_names.append(lli)
            
    return lines_s_num, lines_s_names