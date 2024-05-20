"""Reader for CATHY outputs files
"""

import numpy as np
import pandas as pd
   
    
def read_grid3d(grid3dfile, **kwargs):
    nnod, nnod3, nel = np.loadtxt(grid3dfile, max_rows=1)
    grid3d_file = open(grid3dfile, "r")
    mesh_tetra = np.loadtxt(grid3d_file, skiprows=1, max_rows=int(nel))
    grid3d_file.close()

    grid3d_file = open(grid3dfile, "r")
    mesh3d_nodes = np.loadtxt(
        grid3d_file,
        skiprows=1 + int(nel),
        max_rows=1 + int(nel) + int(nnod3) - 1,
    )
    grid3d_file.close()
    grid = {
        "nnod": nnod,  # number of surface nodes
        "nnod3": nnod3,  # number of volume nodes
        "nel": nel,
        "mesh3d_nodes": mesh3d_nodes,
        "mesh_tetra": mesh_tetra,
    }
    return grid


def read_fort777(filename, **kwargs):
    """
    0  0.00000000E+00     NSTEP   TIME
    SURFACE NODE              X              Y      ACT. ETRA

    Parameters
    ----------
    filename : str
        Path to the fort.777 file.
    Returns
    -------
    xyz_df : pd.DataFrame()
        Dataframe containing the coordinates X,Y and Z (elevation) of the mesh nodes.

    """
    fort777_file = open(filename, "r")
    lines = fort777_file.readlines()

    nstep = []
    surf_node = []
    idx_nstep = []
    time_nstep = []

    # loop over lines of file and identified NSTEP and SURFACE NODE line nb
    # ------------------------------------------------------------------------
    for i, ll in enumerate(lines):
        if "NSTEP" in ll:
            nstep.append(i)
            splt = ll.split(" ")
            wes = [string for string in splt if string != ""]
            idx_nstep.append(wes[0])
            time_nstep.append(float(wes[1]))

        if "SURFACE NODE" in ll:
            surf_node.append([i, wes[1]])

    nsurfnodes = abs(surf_node[0][0] - surf_node[1][0]) - 2

    # build a numpy array
    # ------------------------------------------------------------------------
    df_fort777 = []
    for i, sn_line in enumerate(surf_node):
        tt = np.ones([nsurfnodes]) * time_nstep[i]
        df_fort777.append(
                            np.vstack(
                                        [
                                            tt.T,
                                            np.loadtxt(filename, skiprows=sn_line[0] + 1, max_rows=nsurfnodes).T,
                                        ]
                                    ).T
        )

    df_fort777_stack = np.vstack(np.reshape(df_fort777, [np.shape(df_fort777)[0] * np.shape(df_fort777)[1], 5]))
    # fort777_file collumns information
    # -------------------------------------------------------------------------
    cols_fort777 = ["time_sec","SURFACE NODE","X","Y","ACT. ETRA"]
    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_fort777 = pd.DataFrame(df_fort777_stack, columns=cols_fort777)
    df_fort777["time"] = pd.to_timedelta(df_fort777["time_sec"], unit="s")

    return df_fort777
    
    
def read_xyz(filename, **kwargs):
    """
    Output of grid in 2d as XYZ values coordinates of each nodes

    Parameters
    ----------
    filename : str
        DESCRIPTION.
    Returns
    -------
    xyz_df : pd.DataFrame()
        Dataframe containing the coordinates X,Y and Z (elevation) of the mesh nodes.

    """

    xyz = np.loadtxt(filename,skiprows=1)
    xyz_df = pd.DataFrame(xyz, columns=['id','x','y','z'])

    return xyz_df


def read_wtdepth(filename, **kwargs):
    """
    Output of water table depth variations with times

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    CUMFLOWVOL : TYPE
        DESCRIPTION.

    """

    # wtdepth_df = pd.read_csv(filename, skiprows=1, header=None, sep='\t')
    wtdepth = np.loadtxt(filename)
    wtdepth_df = pd.DataFrame(wtdepth, columns=["time (s)", "water table depth (m)"])

    return wtdepth_df


def read_cumflowvol(filename, **kwargs):
    """
    Output of cumulative flow volumes VSFTOT, VNDTOT, VNNTOT, VNUDTOT, and VTOT

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    CUMFLOWVOL : TYPE
        DESCRIPTION.

    """

    cumflowvol_file = open(filename, "r")
    Lines = cumflowvol_file.readlines()
    count = len(Lines)
    cumflowvol_file.close()

    nstep = count - 3  # Number of timesteps

    cumflowvol_file = open(filename, "r")
    CUMFLOWVOL = np.loadtxt(cumflowvol_file, skiprows=8, max_rows=8 + nstep)
    cumflowvol_file.close()

    return CUMFLOWVOL


def read_vp(filename):
    """
    Vertical profile output in fixed nodes

    Returns
    -------
    None.

    """

    vp_file = open(filename, "r")
    lines = vp_file.readlines()

    nstep = []
    surf_node = []
    idx_s_node = []
    idx_nstep = []
    time_nstep = []

    # loop over lines of file and identified NSTEP and SURFACE NODE line nb
    # ------------------------------------------------------------------------
    for i, ll in enumerate(lines):
        if "NSTEP" in ll:
            nstep.append(i)
            splt = ll.split(" ")
            wes = [string for string in splt if string != ""]
            idx_nstep.append(wes[0])
            time_nstep.append(wes[1])

        if "SURFACE NODE =" in ll:
            surf_node.append([i, wes[1]])

            # here we exctract only numeric value types
            idx_s_node.append([int(s) for s in ll.split() if s.isdigit()])

    # check if surf node is not empty
    # ------------------------------------------------------------------------
    # TO IMPLEMENT
    # try:
    #      surf_node
    # except ValueError:
    #      print("Oops!  That was no valid number.  Try again...")

    nstrat = abs(surf_node[0][0] - surf_node[1][0]) - 3

    # build a numpy array
    # ------------------------------------------------------------------------
    dvp = []
    for i, sn_line in enumerate(surf_node):
        tt = np.ones([nstrat]) * float(sn_line[1])
        nn = np.ones([nstrat]) * idx_s_node[i]
        str_nb = np.arange(1, nstrat + 1)
        dvp.append(
            np.vstack(
                [
                    tt.T,
                    nn.T,
                    str_nb,
                    np.loadtxt(filename, skiprows=sn_line[0] + 2, max_rows=nstrat).T,
                ]
            ).T
        )

    dvp_stack = np.vstack(np.reshape(dvp, [np.shape(dvp)[0] * np.shape(dvp)[1], 9]))

    # Vp collumns information
    # -------------------------------------------------------------------------
    cols_vp = ["time", "node", "str_nb", "z", "PH", "SW", "CKRW", "zeros", "QTRANIE"]
    # PH: pressure head
    # SW: Soil Water
    # CKRW: Relative hydraulic conductivity output at all nodes
    # QTRANIE Root‐zone water uptake at current time level always positive, changes sign in BCPIC
    # dict_vp = {'PH': 'pressure head'}

    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_vp = pd.DataFrame(dvp_stack, columns=cols_vp)
    df_vp["time"] = pd.to_timedelta(df_vp["time"], unit="s")

    return df_vp


def read_hgraph(filename):
    """
    IOUT41 hgraph Surface runoff hydrograph: plot the computed discharge at the outlet (streamflow)

    Returns
    -------
    hgraph data numpy array.

    """

    hgraph_file = open(filename, "r")
    lines = hgraph_file.readlines()
    hgraph_file.close()

    nstep = len(lines) - 2

    hgraph_file = open(filename, "r")
    hgraph = np.loadtxt(hgraph_file, skiprows=2, usecols=range(5))
    hgraph_file.close()

    # hgraph collumns information
    # -------------------------------------------------------------------------
    cols_hgraph = ["time", "streamflow", "Unnamed0", "Unnamed1", "Unnamed2"]

    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_hgraph = pd.DataFrame(hgraph, columns=cols_hgraph)

    return df_hgraph


def read_dtcoupling(filename):
    """
    CPU, time stepping, iteration and other diagnostics of the surface and subsurface modules at each time step

    Returns
    -------
    None.

    """

    dtcoupling_file = open(filename, "r")
    lines = dtcoupling_file.readlines()
    dtcoupling_file.close()

    nstep = len(lines) - 31

    dtcoupling_file = open(filename, "r")
    dtcoupling = np.loadtxt(dtcoupling_file, skiprows=2, max_rows=2 + nstep)
    dtcoupling_file.close()

    df_dtcoupling = pd.DataFrame(dtcoupling)

    df_dtcoupling.columns = [
        "Step",
        "Deltat",
        "Time",
        "Back",
        " NL-l",
        "NL-a",
        "Sdt-l",
        "Sdt-a",
        "Atmpot-vf",
        "Atmpot-v",
        " Atmpot-r",
        "Atmpot-d",
        "Atmact-vf",
        "Atmact-v",
        "Atmact-r",
        "Atmact-d",
        "Horton",
        "Dunne",
        "Ponded",
        "Satur",
        "CPU-sub",
        "CPU-surf",
    ]

    return df_dtcoupling


def read_sw(filename, **kwargs):
    """
    Water Saturation output at all nodes

    Returns
    -------
    None.

    """

    sw_file = open(filename, "r")
    lines = sw_file.readlines()
    sw_file.close()

    # nstep = len(lines)-2

    idx = []
    step_i = []
    time_i = []

    # loop over lines of file and identified NSTEP and SURFACE NODE line nb
    # ------------------------------------------------------------------------
    for i, ll in enumerate(lines):
        if "TIME" in ll:
            # print('TIME')
            idx.append(i)
            splt = ll.split()
            step_i.append(splt[0])
            time_i.append(float(splt[1]))
    idx.append(i + 1)

    sw_sub = []
    for j, ind in enumerate(idx):
        for i, ll in enumerate(lines):
            if i > idx[j] and i < idx[j + 1]:
                splt = ll.split()
                sw_sub.append([float(k) for k in splt])

    sw_sub = np.hstack(sw_sub)
    d_sw_t = np.reshape(
        sw_sub, [len(idx) - 1, int(np.shape(sw_sub)[0] / (len(idx) - 1))]
    )

    df_sw_t = pd.DataFrame(d_sw_t, time_i)
    df_sw_t.index.names = ['Time']
    df_sw_t = df_sw_t[~df_sw_t.index.duplicated(keep='first')]


    # if 'delta_t' in kwargs:
    #     d_sw_t['time'] = pd.to_timedelta(df_hgsfdeth['time'],unit='s')

    # df_hgsfdeth.pivot_table(values='NET SEEPFACE VOL',index='time').plot(ax=ax[0],
    #                                                               ylabel='NET SEEPFACE VOL',
    #                                                               xlabel='time (s)')

    return df_sw_t, df_sw_t.index


def read_psi(filename):
    """
    Pressure head output at all nodes

    Returns
    -------
    None.

    """
    psi_file = open(filename, "r")
    lines = psi_file.readlines()
    psi_file.close()
    idx = []
    step_i = []
    time_i = []

    # loop over lines of file and identified NSTEP and SURFACE NODE line nb
    # ------------------------------------------------------------------------
    for i, ll in enumerate(lines):
        if "TIME" in ll:
            idx.append(i)
            splt = ll.split()
            step_i.append(splt[0])
            time_i.append(splt[1])
    idx.append(i + 1)

    psi_sub = []
    for j, ind in enumerate(idx):
        for i, ll in enumerate(lines):
            if i > idx[j] and i < idx[j + 1]:
                splt = ll.split()
                psi_sub.append([float(k) for k in splt])

    psi_sub = np.hstack(psi_sub)
    d_psi_t = np.reshape(
        psi_sub, [len(idx) - 1, int(np.shape(psi_sub)[0] / (len(idx) - 1))]
    )
    
    # time_i = time_i.astype(float)  # Ensure the index is float
    df_psi_t = pd.DataFrame(d_psi_t, time_i)
    df_psi_t.index = df_psi_t.index.astype(float)
    df_psi_t.index.names = ['Time']
    df_psi_t = df_psi_t[~df_psi_t.index.duplicated(keep='first')]

    return df_psi_t


def read_hgsfdet(filename):
    """
    Detailed seepage face hydrograph output (Incoming and outgoing flows at the seepage face)

    Returns
    -------
    None.

    """

    # hgsfdet_file = open(filename, "r")
    # lines = hgsfdet_file.readlines()
    # hgsfdet_file.close()

    # nstep = len(lines)-2

    hgsfdet_file = open(filename, "r")
    hgsfdet = np.loadtxt(hgsfdet_file, skiprows=5, usecols=range(5))
    hgsfdet_file.close()

    # hgraph collumns information
    # -------------------------------------------------------------------------
    cols_hgsfdet = ["NSTEP", "DELTAT", "time", "NET SEEPFACE VOL", "NET SEEPFACE FLX"]

    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_hgsfdet = pd.DataFrame(hgsfdet, columns=cols_hgsfdet)

    return df_hgsfdet


def read_mbeconv(filename):
    """
    Mass balance and convergence behaviour at each time step
    (REL. MBE (%) should be as small as possible)

    Returns
    -------
    None.

    """
    mbeconv_file = open(filename, "r")
    lines = mbeconv_file.readlines()
    mbeconv_file.close()
    nstep = len(lines) - 2
    mbeconv_file = open(filename, "r")
    mbeconv = np.loadtxt(mbeconv_file, skiprows=2, usecols=range(19))
    mbeconv_file.close()

    # hgraph collumns information
    # -------------------------------------------------------------------------
    cols_mbeconv = [
        "NSTEP",
        "DELTAT",
        "TIME",
        "NLIN",
        "AVG.LIN",
        "STORE1",
        "STORE2",
        "DSTORE",
        "CUM.DSTORE",
        "VIN",
        "CUM. VIN",
        "VOUT",
        "CUM.VOUT",
        "VIN+VOUT",
        "CUM. VIN",
        "M.BAL.ERR",
        "REL. MBE",
        "CUM. MBE",
        "CUM.",
    ]
    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_mbeconv = pd.DataFrame(mbeconv, columns=cols_mbeconv)

    # plt.plot(df_mbeconv['TIME'],df_mbeconv['CUM. MBE'])
    # plt.xlabel('Time (hours)')
    # plt.ylabel('Relative mass balance error (%)')

    return df_mbeconv


def read_netris(filename, **kwargs):

    net_ris_file = open(filename, "r")
    lines = net_ris_file.readlines()
    net_ris_file.close()
    
    netris = []
    for i, ll in enumerate(lines[6:]):
        
        if "INDCELL_WITH_LAKES" in ll:
            break
        lmat = []
        splt = ll.split()
        lmat.append([float(k) for k in splt])
        netris.append(lmat)
    
    netris = np.vstack(netris)
    return netris


def read_cq(filename, **kwargs):
    """
    Reads cell discharge data from a file.

    Parameters
    ----------
    filename : str
        The name of the file to read the data from.

    Returns
    -------
    np.ndarray
        A 2D NumPy array containing the cell discharge data.
    """

    # Open the file
    with open(filename, "r") as cq_file:
        lines = cq_file.readlines()
    
    # Initialize lists to store indices, step numbers, and times
    idx = []
    step_i = []
    time_i = []

    # Loop over lines of the file to identify steps
    for i, line in enumerate(lines):
        if "TIME" in line:
            idx.append(i)
            splt = line.split()
            step_i.append(splt[0])
            time_i.append(float(splt[1]))

    # Initialize a NumPy array to store data
    cq_sub = np.ones([len(idx), 5]) * -99

    # Create a list to store data lines
    lmat = []

    # Loop over lines of the file to extract data
    for line in lines:
        splt = line.split()
        try:
            # Convert data to float and round to 4 decimal places
            # lmat.append([np.round(float(k), 4) for k in splt])
            lmat.append([float(k) for k in splt])
        except ValueError:
            # Skip lines that cannot be converted to float
            pass
    
    # Convert the list of lists to a NumPy array
    lmat = np.vstack(lmat)

    # Reshape the data array
    d_cq_t = np.reshape(lmat, [len(idx), np.diff(idx)[0] - 2, 5])

    # Extract cell discharge data
    cell_discharge = d_cq_t[:, :, 4]
    cell_discharge = cell_discharge.T
    cell_discharge = cell_discharge[:, :-1]
    
    return cell_discharge