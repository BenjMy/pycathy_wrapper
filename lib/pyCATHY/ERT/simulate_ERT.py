"""Class managing ERT data simulation and inversion + petrophysical transformation
"""


import matplotlib.pyplot as plt
import numpy as np
from resipy import R2 # geophysics tools
import pyvista as pv # if not installed : pip install with conda install pyvista (in a conda terminal)
import numpy as np
import os 


        
def create_ERT_survey(pathERT,elecsXYZ,sequence,mesh, **kwargs):
    
    #https://hkex.gitlab.io/resipy/api.html

    # os.chdir(pathERT) 
    
    isExist = os.path.exists(pathERT)

    if not isExist:
      
      # Create a new directory because it does not exist 
      os.makedirs(pathERT) 
    
    ERT = R2(pathERT, typ='R3t') #+ 'ERT_fwdmodel'
    ERT.setTitle('Rhizo_synth')
    
    
    if type(elecsXYZ) is str:
        # elecsXYZ = np.loadtxt(elecsXYZ, delimiter='\t')
        elecsXYZ= np.genfromtxt(elecsXYZ, delimiter=",",skip_header=1)
        
    ERT.setElec(np.c_[elecsXYZ[:,0],elecsXYZ[:,2],elecsXYZ[:,1],elecsXYZ[:,3]])
    
    ERT.importMesh(mesh)
    
    if 'res0' in kwargs:
        ERT.setRefModel(kwargs['res0'])

    ERT.addRegion(np.array([[0.1,0.1],[0.1,0.4],[0.3,0.4],[0.3,0.1],[0.1,0.1]]), 500, iplot=False)
    ERT.importSequence(sequence)
    # -----------------------------------------------
    
    return ERT

def fwd_ERT_survey(ERT,noise,show=False, dump=[]):
    '''
    Fwd ERT model

    Parameters
    ----------
    ERT : TYPE
        DESCRIPTION.
    noise : TYPE
        DESCRIPTION.
    show : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    ERT : TYPE
        DESCRIPTION.

    '''

    # ----------#
    # fwd modelling
    # ---------------------------------------------------------#

    ERT.forward(noise=0.05, iplot=show, dump=dump) # forward modelling with 5 % noise added to the output
    
    return ERT



def invert_ERT_survey(ERT,show=False):

    # ---------------------------------------------------------#
    # inversion
    # ---------------------------------------------------------#
    ERT.param['num_xy_poly'] = 0
    ERT.param['zmin'] = -np.inf
    ERT.param['zmax'] = np.inf
    ERT.param['data_type'] = 1 # using log of resistitivy
    ERT.err = False # if we want to use the error from the error models fitted before
    ERT.param['a_wgt'] = 0.001
    ERT.param['b_wgt'] = 0.05
    
    ERT.invert()
    
    if show==True:
        ERT.showResults(index=0) # show the initial model
        ERT.showResults(index=1, sens=False) # show the inverted model
    


def Archie_ERT(ERT,rFluid,porosity,m,n):
    '''
    Archie transformation
    (only possible on inverted values)

    Parameters
    ----------
    ERT : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    
    # (rho/rFluid/porosity^-m)^(-1/n)
    rFluid=1
    porosity=0.55
    m = 1
    n = 1
    
    str_formula = f"(x['Resistivity']* {rFluid} * {porosity} **(-m))**(-1/ {n})"
    # print(str_formula)
       
    
    # only possible on inverted mesh (meshResults)
    ERT.computeAttribute(str_formula, name='saturation')
    
    pass

