"""Class managing ERT data simulation and inversion + petrophysical transformation
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv # if not installed : pip install with conda install pyvista (in a conda terminal)
import numpy as np
import os 

try: 
    import pygimli as pg
    from pygimli.physics import ert
    import pygimli.meshtools as mt
except ImportError: 
    pygimli = None

try: 
    from resipy import R2 # geophysics tools
except ImportError: 
    resipy = None





# runParallel

def create_ERT_survey_pg(pathERT,sequence,mesh,**kwargs):

    
    noise_level = 5
    if 'noise_level' in kwargs:
        noise_level = kwargs['noise_level']
    
    print('**'*10)
    print(noise_level)
    
    
    # isExist = os.path.exists(pathERT)
    # if not isExist:
    #   # Create a new directory because it does not exist 
    #   os.makedirs(pathERT) 
      
      
    # fname_mesh = './ERT_fwd_DA_sol_rhizo/test_pg/BaseRhizo_Vrte.msh'
    # mesh3d=  mt.readGmsh(mesh, verbose=False)
    
    # pre, ext = os.path.splitext(mesh)
    # print(os.rename(mesh, pre + '.msh'))

    pre, ext = os.path.splitext(mesh)
    
    try:
        mesh3d =  mt.readGmsh(pre + '.msh', verbose=True)
    except:
        try:
            mesh3d =  pg.load(pre + '.bms', verbose=True)
        except:
            raise ValueError('Cannot read '+ mesh + ': valid extensions are .msh and .bms')

    # fname_seq = '/ERT_fwd_DA_sol_rhizo/test_pg/SequenceERT_Rhizo_72Elecs.shm'
    
    # shm = pg.load(fname_seq)
    
    if 'shm' in sequence:
        scheme = pg.load(sequence)
    elif 'txt' in sequence:
        shm = pd.read_csv(sequence, delimiter=' ', header=None)
        shm = shm.to_numpy()
        shm_new = np.subtract(shm,1)
        scheme = pg.DataContainerERT()
        if 'dict_ERT' in kwargs:
            scheme.setSensorPositions(kwargs['dict_ERT']['elecs'])
        for i, elec in enumerate("abmn"):
                scheme[elec] = shm_new[:,i]
    else:
        raise ValueError('Sequence file format not recognized use .shm or .txt as a list of quadrupoles')

        

    
    # scheme["k"] = ert.createGeometricFactors(scheme)
    
    res0 = 1
    if 'res0' in kwargs:
        res0 = kwargs['res0']

    # het = ert.simulate(mesh3d, res=res0, scheme=shm, sr=False, noiseLevel=5,
    #                    calcOnly=True, verbose=True)
    # het.set('k', 1.0/ (hom('u') / hom('i')))
    # het.set('rhoa', het('k') * het('u') / het('i'))
    # het.save('simulated.dat', 'a b m n rhoa k u i')
    
    


        
    if len(res0) != len(mesh3d.cells()):
        raise ValueError('wrong initial resistivity input')
    
    het = ert.simulate(mesh3d, res=res0, scheme=scheme, 
                        calcOnly=False, verbose=False, 
                        noiseLevel=noise_level)


    # het = ert.simulate(mesh3d, res=np.ones(len(mesh3d.cells())), scheme=scheme, 
    #                    calcOnly=False, verbose=False, 
    #                    noiseLevel=noise_level)
    
    # het['rhoa']
    

    return het

        
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

