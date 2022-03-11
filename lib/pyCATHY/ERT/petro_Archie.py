#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""For the manager look at BERT https://gitlab.com/resistivity-net/bert."""

import numpy as np
import matplotlib.pyplot as plt
import os
import pyvista as pv
from pyCATHY import meshtools as mt 
from pyCATHY.rhizo_tools import CATHY_2_Resipy, CATHY_2_pg
from pyCATHY.ERT import simulate_ERT as simuERT
import pandas as pd
from pyCATHY.importers import cathy_outputs as out_CT




def SW_2_ERa(project_name,
             ArchieParms,
             porosity,
             pathERT, meshERT, elecs, sequenceERT,
             path_fwd_CATHY,
             **kwargs):
    '''
    

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
        DESCRIPTION.
    df_Archie : pd df
        Dataframe of the current (given ensemble and given assimilation time) Archie relationship.

    '''
    
    
    
        
    # Some flag for DA assimilation
    # ------------------------------------------------------------------------
    DA_cnb = None
    if 'DA_cnb' in kwargs:
        DA_cnb = kwargs['DA_cnb']
        
    Ens_nb = int(path_fwd_CATHY[-1])
    if 'Ens_nb' in kwargs:
        Ens_nb = kwargs['Ens_nb']
        
    data_format = []
    if 'data_format' in kwargs:
        data_format = kwargs['data_format']
    
    
    savefig = False
    if 'savefig' in kwargs:
        savefig = kwargs['savefig']
        
        

        
        
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


    
        
    path_CATHY = os.path.join(path_fwd_CATHY,'vtk/')

    
       

    if DA_cnb is not None:
        mesh_CATHY_ref = pv.read(os.path.join(path_fwd_CATHY, 'vtk/100.vtk'))
    else:
        mesh_CATHY_ref = pv.read(os.path.join('vtk/100.vtk'))
    
    
    ER_converted_ti = Archie_rho(rFluid=ArchieParms['rFluid'], 
                                sat = df_sw,
                                porosity=porosity, 
                                a=ArchieParms['a'], 
                                m=ArchieParms['m'],
                                n=ArchieParms['n'])

    df_Archie =  pd.DataFrame(columns=['time','ens_nb', 'sw','ER_converted'])
    df_Archie['time'] = DA_cnb*np.ones(len(ER_converted_ti))
    df_Archie['ens_nb'] = Ens_nb*np.ones(len(ER_converted_ti))
    df_Archie['sw'] = df_sw
    df_Archie['ER_converted'] = ER_converted_ti
    df_Archie['porosity'] = porosity*np.ones(len(ER_converted_ti))
    



    # add attribute converted to CATHY mesh
    # ------------------------------------------------------------------------
    mesh_CATHY_new_attr, active_attr = mt.add_attribute_2mesh(ER_converted_ti,
                                                                mesh_CATHY_ref,
                                                                'ER_converted' + str(DA_cnb),
                                                                overwrite=True,
                                                                path = path_CATHY)

    
    
    if 'pygimli' in data_format:
        # copy attribute to simpeg mesh
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
        d = {'a':ERT_predicted['a'], 
              'b':ERT_predicted['b'], 
              'k':ERT_predicted['k'], 
              'm':ERT_predicted['m'], 
              'n':ERT_predicted['n'], 
              'rhoa':ERT_predicted['rhoa'], 
              'valid':ERT_predicted['valid']}
        df_ERT_predicted = pd.DataFrame(data=d)
        
        # save to dataframe and export file
        # ------------------------------------------------------------------------
        filename = 'ER_predicted.csv'
        
        isExist = os.path.exists(os.path.join(pathERT, project_name))
        if not isExist:
            os.makedirs(os.path.join(pathERT, project_name))
    
        df_ERT_predicted.to_csv(os.path.join(pathERT, project_name, filename))
        # ERT_predicted[]
        # simuERT.invert_ERT_survey(ERT)
        
        
    
    # if savefig is True:
        
        
    #     # fig = plt.figure()
        
    #     # plt.scatter(np.arange(0,len(df_sw)), df_sw)
    #     # plotname ='suplot'+ str(DA_cnb)
    #     # plt.savefig(path_CATHY + plotname + '.png', dpi=300)
                
    #     # plotter0 = pv.Plotter(shape=(1, 1),off_screen=True) # notebook = True
    #     # plotter0.subplot(0, 0)
    #     # mesh_CATHY_df, name_new_attr_CATHY = mt.add_attribute_2mesh(df_sw[-1], 
    #     #                                                             mesh_CATHY_ref, 
    #     #                                                             'saturation_df', 
    #     #                                                             overwrite=True)
    #     # mesh_CATHY_df.set_active_scalars('saturation_df')
    #     # my_colormap = 'Blues'
    #     # _ = plotter0.add_mesh(mesh_CATHY_ref,show_edges=True, cmap=my_colormap)

    #     # plotter0.update_scalar_bar_range([0,1]) # max saturation is 1
    #     # plotter0.show_grid()
    #     # plotter0.view_xz()

    #     # plotname ='suplot'+ str(DA_cnb)

        
    #     # plotter0.save_graphic(path_CATHY + plotname + str('.svg'), 
    #     #                     title='', 
    #     #                     raster=True, 
    #     #                     painter=True)   
       
    #     # plotter0.close()
       
    #     plotter = pv.Plotter(shape=(1, 2),off_screen=True) # notebook = True
    #     plotter.subplot(0, 0)
    #     mesh_CATHY_new_attr.set_active_scalars(active_attr)
    #     _ = plotter.add_mesh(mesh_CATHY_new_attr,show_edges=True)
    #     # plotter.show_grid()
    #     # plotter.view_xz()
       
    #     plotter.subplot(0, 1)
    #     mesh_geophy_new_attr.set_active_scalars(scalar_new)
    #     _ = plotter.add_mesh(mesh_geophy_new_attr,show_edges=True)
       
    #     if 'pygimli' in data_format:
    #         plotter.view_xy()
    #     else:      
    #         plotter.view_xz()
       
    #     plotter.show_grid()


    #     plotname ='suplot_ER'+ str(DA_cnb)

    #     # plotter.save_graphic(plotname + str('.svg'), 
    #     #                     title='', 
    #     #                     raster=True, 
    #     #                     painter=True)   
        
    #     plotter.save_graphic(path_CATHY + plotname + str('.svg'), 
    #                         title='', 
    #                         raster=True, 
    #                         painter=True)   
       
    #     plotter.close()
    

    return df_ERT_predicted, df_Archie
  

def Archie_rho(rFluid, sat, porosity, a=1.0, m=2.0, n=2.0):
    '''
    Compute ER values at each mesh nodes

    Parameters
    ----------
    rFluid : int
        conductivity of the pore fluid.
    sat : np.array
        water saturation.
    porosity : int or np.array
        DESCRIPTION.
    a : float, optional
        tortuosity factor. The default is 1.0.
    m : float, optional
        cementation exponent. The default is 2.0.(usually in the range 1.3 -- 2.5 for sandstones)
    n : float, optional
        saturation exponent. The default is 2.0. 

    Returns
    -------
    TYPE
        Electrical Resistivity (Ohm.m).

    '''

    return rFluid * a * porosity**(-m) * sat**(-n)



def Archie_sat(rho, rFluid, porosity, a=1.0, m=2.0, sat=1.0, n=2.0):
    '''
    
    rho: resistivity
    𝑆𝑤 : water saturation
    𝜙: the porosity of the soil
    𝜎_{𝑤} is the conductivity of the pore fluid
    𝑎, 𝑚, and 𝑛 are empirically derived parameters
    𝑎 is the tortuosity factor
    𝑚 is the cementation exponent
    𝑛 is the saturation exponent
    Returns
    -------
    𝑆𝑤 : water saturation


    '''
    return  (rho/rFluid/porosity**(-m))**((-1/n))

