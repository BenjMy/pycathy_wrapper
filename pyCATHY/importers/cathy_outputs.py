"""Reader for CATHY outputs files
"""

import numpy as np
import pandas as pd
from pathlib import Path

    
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


def load_spatial_file_fast(filename: str | Path, prop: str) -> pd.DataFrame:
    def is_number(s):
        try: float(s); return True
        except: return False

    lines = Path(filename).read_text().splitlines()

    # time_nstep, surf_nodes = [], []
    # for i, line in enumerate(lines):
    #     if "NSTEP" in line:
    #         time_nstep.append(next(float(x) for x in line.split() if is_number(x)))
    #     elif "SURFACE NODE" in line:
    #         surf_nodes.append(i)

    # spatial_out_file = open(filename, "r")
    # lines = spatial_out_file.readlines()

    nstep = []
    surf_nodes = []
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
            # surf_nodes.append([i, wes[1]])
            surf_nodes.append(i)
    nsurfnodes = surf_nodes[1] - surf_nodes[0] - 2

    # Instead of np.loadtxt per block: parse directly from lines
    blocks = []
    for t, start in zip(time_nstep, surf_nodes):
        block_str = "\n".join(lines[start+1:start+1+nsurfnodes])
        data = np.fromstring(block_str, sep=" ").reshape(-1, 4)
        tt = np.full((data.shape[0], 1), t)
        blocks.append(np.hstack([tt, data]))

    arr = np.vstack(blocks)
    cols = ["time_sec","SURFACE NODE","X","Y",prop]
    df = pd.DataFrame(arr, columns=cols)
    df["time"] = pd.to_timedelta(df["time_sec"], unit="s")
    return df.drop_duplicates(subset=["time","X","Y"])



def read_spatial_format(filename,prop=None):
    
    spatial_out_file = open(filename, "r")
    lines = spatial_out_file.readlines()

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
    df_spatial = []
    for i, sn_line in enumerate(surf_node):
        tt = np.ones([nsurfnodes]) * time_nstep[i]
        df_spatial.append(
                            np.vstack(
                                        [
                                            tt.T,
                                            np.loadtxt(filename, skiprows=sn_line[0] + 1, 
                                                       max_rows=nsurfnodes).T,
                                        ]
                                    ).T
        )

    df_spatial_stack = np.vstack(np.reshape(df_spatial, 
                                            [np.shape(df_spatial)[0] * np.shape(df_spatial)[1], 5]
                                            )
                                 )
    # fort777_file collumns information
    # -------------------------------------------------------------------------
    colsnames = ["time_sec","SURFACE NODE","X","Y",prop]
    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_spatial = pd.DataFrame(df_spatial_stack, 
                              columns=colsnames
                              )
    df_spatial["time"] = pd.to_timedelta(df_spatial["time_sec"], unit="s")
    df_spatial_multiindex = df_spatial.set_index(['time', 'X', 'Y'])
    df_spatial_unique = df_spatial_multiindex[~df_spatial_multiindex.index.duplicated(keep='first')]
    df_spatial_unique = df_spatial_unique.reset_index()
    
        
    return df_spatial_unique
    
def read_recharge(filename, **kwargs):   
    df_recharge = read_spatial_format(filename,prop='recharge')
    return df_recharge

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
        Dataframe containing time_sec,SURFACE NODE,X,Y,ACT. ETRA.

    """
    
    # df_fort777 = read_spatial_format(filename,prop='ACT. ETRA')
    df_fort777 = load_spatial_file_fast(filename,prop='ACT. ETRA')

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
        "Atmpot-r",
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
    # Step 1: Read the entire file
    with open(filename, "r") as sw_file:
        lines = sw_file.readlines()

    # Step 2: Find the indices and extract step and time information
    step_i, time_i, idx = [], [], []
    for i, line in enumerate(lines):
        if "TIME" in line:
            parts = line.split()
            step_i.append(parts[0])
            time_i.append(float(parts[1]))  # Directly store as float
            idx.append(i)

    # Step 3: Collect all numerical data between identified indices
    num_data = []
    
    if (idx[-1]-idx[-2])==len(lines):
        idx.append(len(lines))  # Add final boundary for last segment
    # else:
    #     print('Incomplete dataset! check simulation results')
        
    num_data = []
    
    for j in range(len(idx) - 1):
        # extract lines for the block
        block_lines = lines[idx[j]+1 : idx[j+1]]
        
        # split all lines and flatten
        block_values = np.fromiter(
            (float(value) for line in block_lines for value in line.split()),
            dtype=float
        )
        num_data.append(block_values)
    
    # concatenate all blocks into a single array
    num_data = np.concatenate(num_data)

    # Step 4: Convert list to numpy array and reshape
    num_array = np.array(num_data, dtype=float)
    num_rows = len(idx) - 1
    num_cols = num_array.size // num_rows
    d_sw_t = num_array.reshape(num_rows, num_cols)

    # Step 5: Create DataFrame
    df_sw_t = pd.DataFrame(d_sw_t, index=time_i[:len(idx)-1])
    df_sw_t.index.name = 'Time'

    # Step 6: Remove duplicate indices, if any
    df_sw_t = df_sw_t[~df_sw_t.index.duplicated(keep='first')]
    
    return df_sw_t, df_sw_t.index


def read_psi(filename):
    """
    Pressure head output at all nodes

    Returns
    -------
    None.

    """
    
    # Step 1: Read the entire file
    with open(filename, "r") as psi_file:
        lines = psi_file.readlines()
    
    # Step 2: Find the indices and extract step and time information
    step_i, time_i, idx = [], [], []
    for i, line in enumerate(lines):
        if "TIME" in line:
            parts = line.split()
            step_i.append(parts[0])
            time_i.append(float(parts[1]))  # Directly store as float
            idx.append(i)
    
    len(step_i)
    # Step 3: Collect all numerical data between identified indices
    num_data = []
    
    if (idx[-1]-idx[-2])==len(lines):
        idx.append(len(lines))  # Add final boundary for last segment
    # else:
    #     print('Incomplete dataset! check simulation results')
    
    # for j in range(len(idx) - 1):
    #     block_data = [line.split() for line in lines[idx[j] + 1:idx[j + 1]]]
    #     num_data.extend(float(value) for line in block_data for value in line)

    num_data = []
    
    for j in range(len(idx) - 1):
        # extract lines for the block
        block_lines = lines[idx[j]+1 : idx[j+1]]
        
        # split all lines and flatten
        block_values = np.fromiter(
            (float(value) for line in block_lines for value in line.split()),
            dtype=float
        )
        num_data.append(block_values)
    
    # concatenate all blocks into a single array
    num_data = np.concatenate(num_data)

    # Step 4: Convert list to numpy array and reshape
    num_array = np.array(num_data, dtype=float)
    num_rows = len(idx) - 1
    num_cols = num_array.size // num_rows
    d_psi_t = num_array.reshape(num_rows, num_cols)
    
    # Step 5: Create DataFrame
    df_psi_t = pd.DataFrame(d_psi_t, index=time_i[:len(idx)-1])
    df_psi_t.index.name = 'Time'
    
    # Step 6: Remove duplicate indices, if any
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


def read_hgatmsf(filename):
    
    # Processes HGATMSF and creates graph:
    # Time (day) - Overland flow (m3/day)
    # Time (day) - Return Flow (m3/day)


    hgatmsf_file = open(filename, "r")
    hgatmsf = np.loadtxt(hgatmsf_file)
    hgatmsf_file.close()

    df_hgatmsf = pd.DataFrame(hgatmsf)

    df_hgatmsf.columns = [
     '#NSTP',      
     'DELTAT',     
     'TIME',   
     'POT. FLUX',
     'ACT. FLUX',
     'OVL. FLUX', # The superficial flow, that is the difference between POT.FLUX and ACT.FLUX;
     'RET. FLUX', # flusso che dal sotterraneo va al superficiale (return flux);  
     'SEEP FLUX',   
     'REC. FLUX',   
     'REC.VOL.',
    ]
    return df_hgatmsf

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