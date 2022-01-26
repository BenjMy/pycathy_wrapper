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

def SW_2_ERa(
             project_name,
             porosity,
             pathERT, meshERT, elecs, sequenceERT,
             path_fwd_CATHY,
             **kwargs):
    '''
    

    Parameters
    ----------

    path_fwd_CATHY : str
        DESCRIPTION.
    project_name : str
        DESCRIPTION.
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
    **kwargs : 
        df_sw : pd df
            saturation mesh values to convert
        savefig : TYPE, optional
            DESCRIPTION. The default is True.
    Returns
    -------
    df_ERT_predicted : TYPE
        DESCRIPTION.
    df_Archie : TYPE
        DESCRIPTION.

    '''
    
    
    
        
    # Some flag for DA assimilation
    # ------------------------------------------------------------------------
    DA_cnb = []
    if 'DA_cnb' in kwargs:
        DA_cnb = kwargs['DA_cnb']
        
    Ens_nb = []
    if 'Ens_nb' in kwargs:
        Ens_nb = kwargs['Ens_nb']
        
    mesh_time_key = []
    if 'mesh_time_key' in kwargs:
        mesh_time_key = kwargs['mesh_time_key']
        
        
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

        # case of sequential assimilatio
        # -------------------------------------------------------------------
        else:
            df_sw = df_sw[-1,:] # take the last time 
            # (in case the sw file contains intermediate times for plotting observations)


        





    print('---------------------------------------------------')

    print('os.path.join(path_fwd_CATHY/output/sw')
    print(os.path.join(path_fwd_CATHY,'output/sw'))

    print('DA_cnb')
    print(DA_cnb)
    
    print('savefig')
    print(savefig)
    
    print('df_sw')
    print(df_sw)
    
        
    print('np.shape(df_sw)')
    print(np.shape(df_sw))
    
    
    print('df_sw[-1,:]')
    print(df_sw[-1])
    print('---------------------------------------------------')

    
        
    path_CATHY = os.path.join(path_fwd_CATHY,'vtk/')


    if mesh_time_key:
        
        if mesh_time_key<10:
            name_mesh= 'cele20' + str(mesh_time_key) + '.vtk'
        else:
            name_mesh= 'cele2' + str(mesh_time_key) + '.vtk'
            
        mesh_CATHY_ref = pv.read(path_CATHY + name_mesh)       
        name_mesh_backup = name_mesh + 'backup_t_ass'  + str(DA_cnb) +'.vtk'
        
    else:
        mesh_CATHY_ref = pv.read(os.path.join(path_CATHY + 'cele200.vtk'))
        name_mesh_backup = 'cele200_backup_t_ass'  + str(DA_cnb) +'.vtk'
    
    mesh_CATHY_ref.save(path_CATHY + name_mesh_backup)

    #sw2convert = mesh_CATHY['saturation']
    ER_converted_ti = Archie_rho(rFluid=1, 
                                sat = df_sw,
                                porosity=porosity, 
                                a=1.0, m=2.0, n=2.0)

    df_Archie =  pd.DataFrame(columns=['time','ens_nb', 'sw','EC','porosity'])
    df_Archie['time'] = DA_cnb
    df_Archie['ens_nb'] = Ens_nb
    df_Archie['sw'] = df_sw[-1]
    df_Archie['EC'] = ER_converted_ti


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
        mesh_geophy_new_attr, scalar_new = CATHY_2_pg(mesh_CATHY_new_attr,meshERT,scalar='ER_converted'+ str(DA_cnb),
                       show=False, path= os.path.join(path_fwd_CATHY, 'vtk/'))
    
    
        # fwd ERT data
        # ------------------------------------------------------------------------
        
        # USING PYGIMLI
        # ------------------------------------------------------------------------
        res0 = mesh_geophy_new_attr.get_array(scalar_new)
        
    
        ERT_predicted = simuERT.create_ERT_survey_pg(os.path.join(pathERT,project_name,'predicted'), 
                                                    sequence=sequenceERT, 
                                                    mesh=meshERT, 
                                                    res0=res0)
        d = {'a':ERT_predicted['a'], 
              'b':ERT_predicted['b'], 
              'k':ERT_predicted['k'], 
              'm':ERT_predicted['m'], 
              'n':ERT_predicted['n'], 
              'rhoa':ERT_predicted['rhoa'], 
              'valid':ERT_predicted['valid']}
        df_ERT_predicted = pd.DataFrame(data=d)
        
    # pg = True   
    # if simpeg==True:
    #     # copy attribute to simpeg mesh
    #     # ------------------------------------------------------------------------
        
    
    else:
        # copy attribute to resipy mesh
        # ------------------------------------------------------------------------
        mesh_geophy_new_attr, scalar_new = CATHY_2_Resipy(mesh_CATHY_new_attr,meshERT,scalar='ER_converted'+ str(DA_cnb),
                       show=False, path= os.path.join(path_fwd_CATHY, 'vtk/'))

        # fwd ERT data
        # ------------------------------------------------------------------------
        
        # USING RESIPY
        # ------------------------------------------------------------------------
        res0 = mesh_geophy_new_attr.get_array(scalar_new)
        # res0 = mesh_geophy_new_attr[scalar_new]
        # mesh_Resipy_new_attr.cell_data[scalar_new] 
        # mesh_Resipy_new_attr.array_names
    
        # ERT.mesh
        # res = ERT.mesh.df['res0']
        ERT_predicted = simuERT.create_ERT_survey(os.path.join(pathERT,project_name,'predicted'), 
                                                    elecs, 
                                                    sequenceERT, 
                                                    meshERT, 
                                                    res0=res0)
        
        ERT_predicted = simuERT.fwd_ERT_survey(ERT_predicted, noise=10)
        df_ERT_predicted = ERT_predicted.surveys[0].df






        # # save to dataframe and export file
        # # ------------------------------------------------------------------------
        # filename = 'ER_predicted.csv'
        
        # isExist = os.path.exists(os.path.join(pathERT, project_name))
        # if not isExist:
        #     os.makedirs(os.path.join(pathERT, project_name))
    
        # df_ERT_predicted.to_csv(os.path.join(pathERT, project_name, filename))
        # # ERT_predicted[]
        # # simuERT.invert_ERT_survey(ERT)
        
        
    
    if savefig is True:
        
        
        fig = plt.figure()
        
        plt.scatter(np.arange(0,len(df_sw)), df_sw)
        plotname ='suplot'+ str(DA_cnb)
        plt.savefig(path_CATHY + plotname + '.png', dpi=300)
                
        # plotter0 = pv.Plotter(shape=(1, 1),off_screen=True) # notebook = True
        # plotter0.subplot(0, 0)
        # mesh_CATHY_df, name_new_attr_CATHY = mt.add_attribute_2mesh(df_sw[-1], 
        #                                                             mesh_CATHY_ref, 
        #                                                             'saturation_df', 
        #                                                             overwrite=True)
        # mesh_CATHY_df.set_active_scalars('saturation_df')
        # my_colormap = 'Blues'
        # _ = plotter0.add_mesh(mesh_CATHY_ref,show_edges=True, cmap=my_colormap)

        # plotter0.update_scalar_bar_range([0,1]) # max saturation is 1
        # plotter0.show_grid()
        # plotter0.view_xz()

        # plotname ='suplot'+ str(DA_cnb)

        
        # plotter0.save_graphic(path_CATHY + plotname + str('.svg'), 
        #                     title='', 
        #                     raster=True, 
        #                     painter=True)   
       
        # plotter0.close()
       
        plotter = pv.Plotter(shape=(1, 2),off_screen=True) # notebook = True
        plotter.subplot(0, 0)
        mesh_CATHY_new_attr.set_active_scalars(active_attr)
        _ = plotter.add_mesh(mesh_CATHY_new_attr,show_edges=True)
        # plotter.show_grid()
        plotter.view_xz()
       
        plotter.subplot(0, 1)
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

