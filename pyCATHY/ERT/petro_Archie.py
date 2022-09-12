"""Pedophysical transformations"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pyvista as pv
from pyCATHY import meshtools as mt 
from pyCATHY.rhizo_tools import CATHY_2_Resipy, CATHY_2_pg
from pyCATHY.ERT import simulate_ERT as simuERT
import pandas as pd
from pyCATHY.importers import cathy_outputs as out_CT


def get_sw_ens_i(Ens_nb,path_fwd_CATHY,**kwargs):
    ''' Return sw array for a given ensemble (if ensemble exist)'''
    # Get the state vector from kwargs if existing otherwise read output
    # ------------------------------------------------------------------------
    if 'df_sw' in kwargs:
        df_sw = kwargs['df_sw']
    else:
        df_sw = out_CT.read_sw(os.path.join(path_fwd_CATHY,'output/sw'))
        # case of open Loop where df_sw is list of nb of assimilation time length
        # -------------------------------------------------------------------
        if 'time_ass' in kwargs:
            time_ass = kwargs['time_ass']
            df_sw = df_sw[time_ass,:]
            DA_cnb = time_ass
        # case of sequential assimilation
        # -------------------------------------------------------------------
        else:
            df_sw = df_sw[-1,:] # take the last time 
            # (in case the sw file contains intermediate times for plotting observations)
    return df_sw


def get_Archie_ens_i(ArchieParms,Ens_nb):
    #     This should be moved to DA Class
    ''' Return Archie parameter for a given ensemble (if ensemble exist)'''
    ArchieParms2parse = {}
    if len(ArchieParms['rFluid_Archie'])>1:
        for p in ['porosity', 'rFluid_Archie','a_Archie','m_Archie','n_Archie','pert_sigma_Archie']:
            ArchieParms2parse[p] = [ArchieParms[p][Ens_nb]]
    else:
        ArchieParms2parse = ArchieParms
    return ArchieParms2parse

    
def SW_2_ERa_DA(project_name,
                 ArchieParms,
                 porosity,
                 pathERT, meshERT, elecs, sequenceERT,
                 path_fwd_CATHY,
                 **kwargs):
    '''
    This should be moved to DA Class
    Data Assimilation
    -----------------
    Map saturation water from CATHY model to apparent Electrical Resistivities
    1. Import/read CATHY sw file result
    2. Convert SW to ER uisng Archie law
    3. Fwd ERT model to obtain apparent resistivities (predicted observation for DA)
    
    Parameters
    ----------
    project_name : str
        name of the current project.
    ArchieParms : dict
        Archie params.
    porosity : np.array([])
        medium porosity.
    pathERT : str
        path of the ERT forward mesh.
    meshERT : str
        filename of the ERT fwd mesh.
    elecs : TYPE
        electrode positions.
    sequenceERT : TYPE
        ERT sequence.
    path_fwd_CATHY: list of folder containing CATHY simulation outputs
        Use for // simulations
    **kwargs : 
        df_sw : pd df
            saturation mesh values to convert
        savefig : TYPE, optional
            DESCRIPTION. The default is True.
    Returns
    -------
    df_ERT_predicted : pd df
        Dataframe of predicted ER.
    df_Archie : pd df
        Dataframe of the current (given ensemble and given assimilation time) Archie relationship.

    '''
    
    path_CATHY = os.path.join(path_fwd_CATHY,'vtk/')
    
    # Some flag for DA assimilation
    # ------------------------------------------------------------------------
    DA_cnb = None
    Ens_nb = int(os.path.split(path_fwd_CATHY)[-1].split("_",1)[1])-1
    
    data_format = []
    savefig = False

    if 'DA_cnb' in kwargs:
        DA_cnb = kwargs['DA_cnb']
    if 'Ens_nb' in kwargs:
        Ens_nb = kwargs['Ens_nb']
    if 'data_format' in kwargs:
        data_format = kwargs['data_format']
    if 'savefig' in kwargs:
        savefig = kwargs['savefig']

        
    # Get sw array for a given ensemble 
    # ------------------------------------
    df_sw = get_sw_ens_i(Ens_nb,path_fwd_CATHY,**kwargs)

    # Read the input mesh using pyvista
    # ------------------------------------
    if DA_cnb is not None:
        mesh_CATHY_ref = pv.read(os.path.join(path_fwd_CATHY, 'vtk/100.vtk'))
    else:
        mesh_CATHY_ref = pv.read(os.path.join('vtk/100.vtk'))
        
    # Choose archie parameter for a given realisation (from the ensemble)
    # --------------------------------------------------------------------
    ArchieParms2parse = get_Archie_ens_i(ArchieParms,Ens_nb)
       
    # Convert to SW to ER values
    # --------------------------------------------------------------------
    ER_converted_ti = Archie_rho_DA(rFluid_Archie=ArchieParms2parse['rFluid_Archie'], 
                                     sat = [df_sw],
                                     porosity=ArchieParms2parse['porosity'], 
                                     a_Archie=ArchieParms2parse['a_Archie'], 
                                     m_Archie=ArchieParms2parse['m_Archie'],
                                     n_Archie=ArchieParms2parse['n_Archie'],
                                     pert_sigma_Archie=ArchieParms2parse['pert_sigma_Archie']
                                    )

    df_Archie =  pd.DataFrame(columns=['time','ens_nb', 'sw','ER_converted'])
    df_Archie['time'] = DA_cnb*np.ones(len(ER_converted_ti))
    df_Archie['ens_nb'] = Ens_nb*np.ones(len(ER_converted_ti))
    df_Archie['sw'] = df_sw
    df_Archie['ER_converted'] = ER_converted_ti
    df_Archie['porosity'] = ArchieParms2parse['porosity']*np.ones(len(ER_converted_ti))
    
    # add attribute converted to CATHY mesh
    # ------------------------------------------------------------------------
    mesh_CATHY_new_attr, active_attr = mt.add_attribute_2mesh(ER_converted_ti,
                                                                mesh_CATHY_ref,
                                                                'ER_converted' + str(DA_cnb),
                                                                overwrite=True,
                                                                path = path_CATHY)

    if 'pygimli' in data_format:
        # copy attribute to pg mesh
        # ------------------------------------------------------------------------
        

        mesh_geophy_new_attr, scalar_new = CATHY_2_pg(mesh_CATHY_new_attr,
                                                      meshERT,
                                                      scalar='ER_converted'+ str(DA_cnb),
                                                      show=False,
                                                      path= path_CATHY,
                                                      **kwargs
                                                      )
    
        # fwd ERT data
        # ------------------------------------------------------------------------
        # USING PYGIMLI
        # ------------------------------------------------------------------------
        res0 = mesh_geophy_new_attr.get_array(scalar_new)
                    

        ERT_predicted = simuERT.create_ERT_survey_pg(os.path.join(pathERT,
                                                                  project_name,
                                                                  'predicted'), 
                                                     sequence=sequenceERT, 
                                                     mesh=meshERT, 
                                                     res0=res0,
                                                     **kwargs)
        # ERT_predicted
        d = {'a':ERT_predicted['a'], 
              'b':ERT_predicted['b'], 
              'k':ERT_predicted['k'], 
              'm':ERT_predicted['m'], 
              'n':ERT_predicted['n'], 
              'rhoa':ERT_predicted['rhoa'], 
              'valid':ERT_predicted['valid']}
        df_ERT_predicted = pd.DataFrame(data=d)
        

    if savefig:

        # fig = plt.figure()
        # plt.scatter(np.arange(0,len(df_sw)), df_sw)
        # plotname ='suplot'+ str(DA_cnb)
        # plt.savefig(path_CATHY + plotname + '.png', dpi=300)
                
        plotter = pv.Plotter(shape=(3, 1),off_screen=True) # notebook = True
        plotter.subplot(0, 0)
        
        mesh_CATHY_df, name_new_attr_CATHY = mt.add_attribute_2mesh(df_sw[-1], 
                                                                    mesh_CATHY_ref, 
                                                                    'saturation_df', 
                                                                    overwrite=True)
        mesh_CATHY_df.set_active_scalars('saturation_df')
        my_colormap = 'Blues'
        _ = plotter.add_mesh(mesh_CATHY_ref,show_edges=True, cmap=my_colormap)

        plotter.update_scalar_bar_range([0,1]) # max saturation is 1
        plotter.show_grid()
        plotter.view_xz()
       
        plotter.subplot(1, 0)
        mesh_CATHY_new_attr.set_active_scalars(active_attr)
        _ = plotter.add_mesh(mesh_CATHY_new_attr,show_edges=True)
        plotter.show_grid()
        # plotter.view_xz()
       
        plotter.subplot(2, 0)
        mesh_geophy_new_attr.set_active_scalars(scalar_new)
        _ = plotter.add_mesh(mesh_geophy_new_attr,show_edges=True)
       
        if 'pygimli' in data_format:
            plotter.view_xy()
        else:      
            plotter.view_xz()
       
        plotter.show_grid()
        plotname ='suplot_ER'+ str(DA_cnb)       
        plotter.save_graphic(path_CATHY + plotname + str('.svg'), 
                            title='', 
                            raster=True, 
                            painter=True)   
       
        plotter.close()
    

    return df_ERT_predicted, df_Archie
  

def Archie_rho_DA(rFluid_Archie=[], sat=[], porosity=[], a_Archie=[1.0], m_Archie=[2.0], n_Archie=[2.0],
               pert_sigma_Archie=[0]):
    #     This should be moved to DA Class
    '''
    Compute ER values at each mesh nodes. 
    If pert_sigma, add a nornal noise of sigma standart deviation and 0 mean (# See eq. 4.4 thesis Isabelle p.95)

    Parameters
    ----------
    rFluid : list
        Conductivity of the pore fluid for each soil zone
    sat : list
        Water saturation for each soil zone
    porosity : list
        Porosity for each soil zone
    a : list, optional
        Tortuosity factor. The default is 1.0.
    m : list, optional
        Cementation exponent. The default is 2.0.(usually in the range 1.3 -- 2.5 for sandstones)
    n : list, optional
        Saturation exponent. The default is 2.0. 
    pert_sigma_Archie: list, optional
        Normal noise to add to the data

    Returns
    -------
    rho: Converted Apparent Electrical Resistivity (Ohm.m).

    '''
    
    # Loop over soil type
    # -----------------------------------------------
    for i in range(len(rFluid_Archie)): # loop over Archie heterogeneity i.e. soil zones
        # print('!!! shortcut sat[:] is not valid for heteregeneous soil!')
        rho = rFluid_Archie[i] * a_Archie[i] * porosity[i]**(-m_Archie[i]) * sat[i]**(-n_Archie[i])
        sigma = 1/rho
        
        if sat[i].ndim>1:
            for ti in range(np.shape(sat)[0]): # Loop over assimilation times
                for meas_nb in range(np.shape(sat)[1]): # Loop over mesh nodes
                    noise = np.random.normal(0,(pert_sigma_Archie[i]*(sigma[ti,meas_nb]))**2,1) # See eq. 4.4 thesis Isabelle p.95
                    sigma[ti,meas_nb] = sigma[ti,meas_nb] + noise
        else:
            for meas_nb in range(len(sat[i])): # Loop over mesh nodes
                noise = np.random.normal(0,(pert_sigma_Archie[i]*(sigma[meas_nb]))**2,1) # See eq. 4.4 thesis Isabelle p.95
                sigma[meas_nb] = sigma[meas_nb] + noise
            
    
    return (1/sigma)





