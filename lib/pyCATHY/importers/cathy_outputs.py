"""Reader for CATHY outputs files
"""

import numpy as np
import pandas as pd

def read_cumflowvol(filename,**kwargs):
    '''
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

    '''
    
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
    '''
    Vertical profile output in fixed nodes

    Returns
    -------
    None.

    '''
    
    
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
        if 'NSTEP' in ll:
            nstep.append(i)
            splt= ll.split(" ")
            wes = [string for string in splt if string != ""]
            idx_nstep.append(wes[0])
            time_nstep.append(wes[1])
    
        if 'SURFACE NODE =' in ll:
            surf_node.append([i,wes[1]])
            
            # here we exctract only numeric value types
            idx_s_node.append([int(s) for s in ll.split() if s.isdigit()]) 
  
    
    # check if surf node is not empty
    # ------------------------------------------------------------------------
    # TO IMPLEMENT
    # try:
    #      surf_node
    # except ValueError:
    #      print("Oops!  That was no valid number.  Try again...")
         
                 
        
    nstrat = abs(surf_node[0][0] - surf_node[1][0])- 3

    
    # build a numpy array
    # ------------------------------------------------------------------------
    dvp = []
    for i, sn_line in enumerate(surf_node):   
        tt = np.ones([nstrat])*float(sn_line[1])
        nn = np.ones([nstrat])*idx_s_node[i]
        str_nb = np.arange(1,nstrat+1)
        dvp.append(np.vstack([tt.T, nn.T, str_nb,
                    np.loadtxt(filename,skiprows=sn_line[0]+2,max_rows=nstrat).T]).T)
        
    dvp_stack = np.vstack(np.reshape(dvp,[np.shape(dvp)[0]*np.shape(dvp)[1],9]))
    
    
    # Vp collumns information
    # -------------------------------------------------------------------------
    cols_vp = ['time', 'node', 'str_nb',  'z', 'PH', 'SW', 'CKRW', 'zeros', 'QTRANIE']
    # PH: pressure head
    # SW: Soil Water
    # CKRW: Relative hydraulic conductivity output at all nodes
    # QTRANIE Root‐zone water uptake at current time level always positive, changes sign in BCPIC
    # dict_vp = {'PH': 'pressure head'}


    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_vp = pd.DataFrame(dvp_stack,  columns=cols_vp)
    df_vp['time'] = pd.to_timedelta(df_vp['time'],unit='s') 

    return   df_vp  


def read_hgraph(filename):
    '''
    IOUT41 hgraph Surface runoff hydrograph: plot the computed discharge at the outlet (streamflow)

    Returns
    -------
    hgraph data numpy array.

    '''
    
    
    hgraph_file = open(filename, "r")
    lines = hgraph_file.readlines()
    hgraph_file.close()
    
    
    nstep = len(lines)-2
    
    hgraph_file = open(filename, "r")
    hgraph = np.loadtxt(hgraph_file, skiprows=2, usecols=range(5))
    hgraph_file.close()

    # hgraph collumns information
    # -------------------------------------------------------------------------
    cols_hgraph = ['time', 'streamflow', '',  '', '']
    
    
    # transform a numpy array into panda df
    # ------------------------------------------------------------------------
    df_hgraph = pd.DataFrame(hgraph,  columns=cols_hgraph)

    return df_hgraph  



def read_dtcoupling(filename):
    '''
    CPU, time stepping, iteration and other diagnostics of the surface and subsurface modules at each time step
    
    Returns
    -------
    None.

    '''
    
    dtcoupling_file = open(filename, "r")
    lines = dtcoupling_file.readlines()
    dtcoupling_file.close()
    
    
    nstep = len(lines)-31
    
    dtcoupling_file = open(filename, "r")
    dtcoupling = np.loadtxt(dtcoupling_file, skiprows=2,max_rows=2+nstep)
    dtcoupling_file.close()
    
    df_dtcoupling = pd.DataFrame(dtcoupling)
    
    
    df_dtcoupling.columns = ['Step','Deltat','Time','Back',' NL-l','NL-a','Sdt-l','Sdt-a',
                  'Atmpot-vf','Atmpot-v',' Atmpot-r','Atmpot-d','Atmact-vf','Atmact-v','Atmact-r','Atmact-d',
                  'Horton','Dunne','Ponded','Satur','CPU-sub','CPU-surf']
    
    return df_dtcoupling   


def read_psi(filename):
    '''
    Pressure head output at all nodes

    Returns
    -------
    None.

    '''

    psi_file = open(filename, "r")
    lines = psi_file.readlines()
    psi_file.close()
    
    # nstep = len(lines)-2
    
    idx=[]
    step_i = []
    time_i = []
    
    # loop over lines of file and identified NSTEP and TIME line nb
    # ------------------------------------------------------------------------
    for i, ll in enumerate(lines):
        if 'TIME' in ll:
            idx.append(i)
            splt= ll.split()
            step_i.append(splt[0])
            time_i.append(splt[1])
    
    
    d_psi_t = []
    for i, ind in enumerate(idx):
        psi_file = open(filename, "r")
        psi_sub = np.loadtxt(psi_file, skiprows=idx[i]+1,max_rows=idx[1]-2)
        d_psi_t.append(psi_sub)
        psi_file.close()
    
    d_psi_t = np.vstack(d_psi_t)

    
    
    return  d_psi_t


def read_sw():
    '''
    Water saturation output at all nodes (SW) for input to TRAN3D and DUAL3D groundwater contaminant transport

    Returns
    -------
    None.

    '''
    
    
    pass   


