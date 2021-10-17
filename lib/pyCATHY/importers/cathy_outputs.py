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
    
def read_vp():
    '''
    Vertical profile output in fixed nodes

    Returns
    -------
    None.

    '''
    filename = 'vp'
    vp_file = open(filename, "r")
    lines = vp_file.readlines()
    
    nstep = []
    surf_node = []
    idx_s_node = []
    idx_nstep = []
    time_nstep = []
    
    for i, ll in enumerate(lines):
        if 'NSTEP' in ll:
            nstep.append(i)
            # idx_nstep.append([int(s) for s in ll.split() if s.isdigit()])
            splt= ll.split(" ")
            wes = [string for string in splt if string != ""]
            idx_nstep.append(wes[0])
            time_nstep.append(wes[1])
    
        if 'SURFACE NODE =' in ll:
            surf_node.append([i,wes[1]])
            idx_s_node.append([int(s) for s in ll.split() if s.isdigit()])
    idx_s_node = np.unique(idx_s_node)
        
    
    nstrat = abs(surf_node[0][0] - surf_node[1][0])- 3
    
    dvp = []
    for i, sn_line in enumerate(surf_node):   
        tt = np.ones([nstrat])*float(sn_line[1])
        nn = np.ones([nstrat])*sn_line[0]
        str_nb = np.arange(1,nstrat+1)
        dvp.append(np.vstack([tt.T, nn.T, str_nb,
                    np.loadtxt(filename,skiprows=sn_line[0]+2,max_rows=nstrat).T]).T)
        
    dvp_stack = np.vstack(np.reshape(dvp,[np.shape(dvp)[0]*np.shape(dvp)[1],9]))
    
    cols_vp = ['time', 'node', 'str_nb',  'z', 'pressure_head', 'SW', 'CKRW', 'zeros', 'QTRANIE']
    
    # Relative hydraulic conductivity output at all nodes
    # Root‐zone water uptake at current time level always positive, changes sign in BCPIC
    
    df_vp = pd.DataFrame(dvp_stack,  columns=cols_vp)
        
    return   df_vp  

def read_psi():
    '''

    Returns
    -------
    None.

    '''
    
    
    pass   


def read_sw():
    '''

    Returns
    -------
    None.

    '''
    
    
    pass   


