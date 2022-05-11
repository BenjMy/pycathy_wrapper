# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import numpy as np

# from pyCATHY import cathy_tools as CT
from pyCATHY.plotters import cathy_plots as pltCT
from pyCATHY import meshtools as mt
from pyCATHY.DA import perturbate

import matplotlib.pyplot as plt


import os
import pyvista as pv
# import pygimli as pg
import pandas as pd

class rhizotron(object):
    '''rhizotron class'''
    
    def __init__(self):
        
        # survey
        # ---------------------------
        self.dates_rhizo = []

        # Petrophysics 
        # ---------------------------
        self.Archie_rhizo = []
        self.VGN_rhizo = []
        
        # Potential evapotranspiration
        # ---------------------------
        self.weight_rhizo = []
        self.day_temperature = []
        self.relative_humidity = []
        self.timeOnLight = []
        self.timeOffLight = []
        self.photoperiod = []


        pass


    # Petrophysical relationships
    # -------------------------------------------------------------------------
    def Archie_rhizo(self):
        pass
    
    def VGN_rhizo(self):
        pass
    
    
    

        
    # ETp estimation
    # -------------------------------------------------------------------------
    def estimate_daily_ETp(self):
        self.set_photoperiod()
        self.set_global_radiation()
        self.set_relative_humidity()
        self.set_day_temperature()
        self.estimated_ETp_PenmanMonteith()
        self.import_porometer_data()
        
    def set_photoperiod(self):
        self.timeOnLight = 7 #%
        self.timeOffLight = 19 #%
        self.photoperiod = self.timeOffLight - self.timeOnLight # in hours
    def set_global_radiation(self):
        self.global_radiation = [300,350] # micromol - PAR m^-2 s^-1 (equivalent to ~150–175 W m^2)
    def set_relative_humidity(self):
        self.relative_humidity = 60 #%
    def set_day_temperature(self):
        self.day_temperature = [20,25] # ° Celcius
    def estimated_ETp_PenmanMonteith(self):
        # empiricall values from Guarrigues : 10.1007/s11104-004-7903-0
        # see https://arww.razi.ac.ir/article_792_b4a7c50314adca55aadfe74b33ef1584.pdf 
        # for equations of PenmanMonteith
        self.ETp_mm_day = 3.6 #mm.day^-1       
    def import_porometer_data(self):
        # porometer
        pass

    def estimate_hourly_ETp(self):
        self.df_weight['ETp_mm_hours'] = np.zeros(len(self.df_weight['abs_date_new']))
        for d in self.df_weight['diurn']:
            if d:
                self.df_weight['ETp_mm_sec'] = (self.ETp_mm_day/self.photoperiod) * (1/86400)
                
    
    # weight sensors
    # -------------------------------------------------------------------------
    
    def filter_date_range(self, date0=None,date1=None):
        if date0 is None:
            date0 = self.df_weight['abs_date_new'][0]
            date1 = self.df_weight['abs_date_new'].iloc[-1]
            # delta = date1 - date0
        
        #mask = (self.df_weight['abs_date_new'] > date0) & (self.df_weight['abs_date_new'] <= date1)
        #self.df_weight = self.df_weight['weight (kg)'].loc[mask]
            
            
    def concat_weight_dates(self,sheet_name, sheet_id, initial_date):
        list_data_scale = []
        # sheet_id = sheet_id
        
        for sn in sheet_name:
            url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sn}'
            data_scale_i = pd.read_csv(url, decimal=',')
            
            data_scale_i['abs_date'] = pd.to_datetime(data_scale_i['date'] + " " + data_scale_i['time'],
                                                infer_datetime_format=True)             
            # datestr = data_scale_i['date'] + " " + data_scale_i['time']
            data_scale_i['abs_date'] = pd.to_datetime(data_scale_i['date'] + " " + data_scale_i['time'],
                                            format="%d/%m/%Y %H:%M") 
            if initial_date is None:
                initial_date = data_scale_i['abs_date'][0]
            
            dates = []
            for i in range(len(data_scale_i['time'])):
                dates.append(initial_date + pd.Timedelta(seconds=data_scale_i['sec'][i]))
                
            date_time = [t.strftime("%Y-%m-%d %H:%M:%S") for t in dates]
            data_scale_i['abs_date_new']=date_time
            data_scale_i['abs_date_new']= pd.to_datetime(data_scale_i['abs_date_new'])
            list_data_scale.append(data_scale_i)
            
            self.df_weight = pd.concat(list_data_scale, axis=0, ignore_index=False)
            
    
    def add_diurn_flag(self):
        diurn=[]
        for d in self.df_weight['abs_date_new']:
            if int(d.strftime('%H'))>=self.timeOnLight and int(d.strftime('%H'))<=self.timeOffLight:
                diurn.append(True)
            else:
                diurn.append(False)
        
        self.df_weight['diurn'] = diurn  


    def set_irr_patern(self):
        # date, position
        
        # -------------------------------+
        # This is a test implementation
        self.drip_loc_df = pd.DataFrame()
        my_dict = {'2022-04-06 17:41:00':'hole1','2022-04-10 11:55:57':'hole2'}
        self.drip_loc_df = pd.DataFrame(list(my_dict.items()),columns = ['date','drip_loc'])
        
        self.dict_drip_loc_xyz = {
                                        'hole1': [0.1,0,0],
                                        'hole2': [0.4,0,0],            
                                }
        
    def check_nb_irrigation(self, irr):
        # date, position
        true_count = sum(irr)
        if true_count == len(self.drip_loc_df):
            valid_nb = True
        else:
            valid_nb = False
        return valid_nb
    
    
    def add_irr_flag(self, threshold=5e-2):
        self.set_irr_patern()
        irr=[]
        valid_nb = False
        while valid_nb is False:
            for d in self.df_weight['diff']:
                if d>threshold:
                    irr.append(True)
                else:
                    irr.append(False)
            valid_nb  = self.check_nb_irrigation(irr)
            
            # -------------------------------+
            # This is a test implementation
            valid_nb = True
        self.df_weight['irr_flag'] = irr  


    def units_ET_m_s(self):   
        #rhizo_surface_m2 = 0.03*0.5
        single_leaf_surface_area = 5e-3*5e-3
        nb_of_leaves = 10
        leaf_surf_area = single_leaf_surface_area*nb_of_leaves

        #elapsed_time_seconds= self.df_weight['weight (kg)'].elapsed.dt.total_seconds() 
        ET = (self.df_weight['weight (kg)'])/leaf_surf_area
    
    def estimate_uncertainties_weight(self):
        # Add error estimates on measurements
        pass
        
    def process_weight_sensor(self):
        '''
        Discriminate between irrigation and evapotranspiration from scale data
        Assuming that plant growth is negligeable
        '''            
        self.df_weight['cumweight'] = self.df_weight['weight (kg)'].cumsum()
        self.df_weight['rolling'] = self.df_weight['weight (kg)'].rolling(12).mean()
        self.df_weight['rolling_std'] = self.df_weight['weight (kg)'].rolling(12).std()
        self.df_weight['diff'] = self.df_weight['weight (kg)'].diff()
        
        # df_weight_pos_diff = self.df_weight[self.df_weight['weight_diff'] >=0]
        
        self.add_diurn_flag()
        self.add_irr_flag()
        self.units_ET_m_s()
        self.estimate_uncertainties_weight()

    
    def import_weight_sensor_data(self,initial_date=None):
        '''
        Import from google datasheet directly
        '''
        sheet_id = '1GqDKU8fC0Pw_DQpUD-MBLpN-AXWrAzLuO1hABxxYIhc'
        sheet_name = ['apr_6-11','apr_11-13','apr_13-21','apr_21-26','apr_28-may_2']
        self.concat_weight_dates(sheet_name, sheet_id, initial_date)
        self.filter_date_range()
        self.process_weight_sensor()

        
    def add_irr_patern(self):
        
        
        masse_vol_water = 1e3 # kg.m^3
        drip_hole_surf = 1e-2*1e-2
        elapsed_time_seconds = self.df_weight['sec'].diff() 
        
        self.df_weight['ETp_irr_mm_s']  = np.zeros(len(self.df_weight['sec']))
        
        drip_loc = np.empty(len(self.df_weight['sec']))
        drip_loc[:] = np.nan
        self.df_weight['drip_loc']  = drip_loc

        for i, d in enumerate(self.df_weight['irr_flag']):
            if d:
                self.df_weight['ETp_irr_mm_s'].iloc[i] = ( 
                                                    (self.df_weight['diff'].iloc[i] / masse_vol_water) *
                                                    (drip_hole_surf/elapsed_time_seconds.iloc[i]) * 1e3
                                                )
                
                # -------------------------------+
                # This is a test implementation
                if (i % 2) == 0: 
                    self.df_weight['drip_loc'].iloc[i] = self.drip_loc_df['drip_loc'][0]
                else:
                    self.df_weight['drip_loc'].iloc[i] = self.drip_loc_df['drip_loc'][1]  
        pass

    def net_flux(self):
        self.df_weight['ETP_net'] = self.df_weight['ETp_irr_mm_s'] - self.df_weight['ETp_mm_sec']

    def prepare_atmbc_from_weight(self, simu):
        self.estimate_daily_ETp()
        self.import_weight_sensor_data()
        self.estimate_hourly_ETp()
        self.add_irr_patern()
        self.net_flux()
        
        
    def set_defaults(self,simu):
        
        # MESH 
        # *****
        
        dem_synth=np.ones([9,10])*1e-3
        dem_synth[-1,-1]=0.99e-3
        zb = np.arange(0,0.5+0.1,0.05)
        nstr=len(zb)-1
        zr=list((np.ones(len(zb)))/(nstr))
        simu.update_prepo_inputs(delta_x=0.05,delta_y=0.003, DEM=dem_synth,
                                  N=np.shape(dem_synth)[1],
                                  M=np.shape(dem_synth)[0],
                                  N_celle=np.shape(dem_synth)[0]*np.shape(dem_synth)[1],
                                  nstr=nstr, zratio=zr,  base=max(zb),dr=0.05,show=True)
        simu.run_preprocessor(verbose=True,
                                 KeepOutlet=False) #,show=True)
        simu.update_parm()
        simu.update_veg_map()
        simu.update_soil()
        
        self.make_mesh(simu)
        
        
        pass
    
    
    def make_mesh(self,simu):
        simu.run_processor(IPRT1=3,
                              TRAFLAG=0,
                              verbose=False)
        self.grid = simu.grid
        
        pass
    
    
    def update_atmbc_PRD(self,simu):
        '''   
        '''
        irr = []
        ET = []
        
   
        x_min = max(self.grid["nodes_idxyz"][:, 1]) / 2
        x_max = max(self.grid["nodes_idxyz"][:, 1])
        y_min = min(self.grid["nodes_idxyz"][:, 2])
        y_max = max(self.grid["nodes_idxyz"][:, 2])
    
        # nodes_xyz = grid["nodes_idxyz"]
    
        xmesh = self.grid["nodes_idxyz"][:, 1]
        ymesh = self.grid["nodes_idxyz"][:, 2]
        zmesh = self.grid["nodes_idxyz"][:, 3]
    
        HSPATM = 0
        IETO = 0
    
        # flux = 1.11111E-06
        totarea = 1  # 0.0
        no_irr = 0

        
        with open(os.path.join(simu.workdir, simu.project_name, "input/atmbc"), "w+") as atmbcfile:

            # write no irrigation during the first delta_drying0 days
            atmbcfile.write(
                str(HSPATM) + "\t" + str(IETO) + "\t" + "HSPATM" + "\t" + "IETO" + "\n"
            )
            atmbcfile.write(str(0) + "\t" + "TIME" + "\n")

            count_nmax = 0
            for i in range(int(self.grid["nnod3"])):
                if round(zmesh[i], 3) == round(max(zmesh), 3):
                    count_nmax += 1

            # loop over assimilation times
            # --------------------------------------
            for d in self.df_weight['abs_date_new']:


                # loop over the nodes of the mesh
                # --------------------------------------
    
                count = 0
                for i in range(int(self.grid3d["nnod3"])):
                    if round(zmesh[i], 3) == round(max(zmesh), 3):
                        count += 1
                        if self.dict_drip_loc_xyz[0]:
                            pass
                            


        atmbcfile.close()
            
        
        return irr, ET
    
    
    
    
    # def update_atmbc_PRD(self, simu):

        # simu.grid3d
        # self.update_atmbc_PRD(simu.workdir,
        #                       simu.project_name,
        #                       simu.grid)


        # self.root_surface_area()
        # self.set_drippers()
        # self.create_infitration()
        # self.update_rhizo_inputs()
        # self.perturbate()

    

def atmbc_PRD_synth(self,workdir,project_name,dict_PRD,show=False,**kwargs):
    '''
    Synthetic PRD
    

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.
    dict_PRD : TYPE
        DESCRIPTION.
    show : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    t_irr : TYPE
        DESCRIPTION.
    t_atmbc : TYPE
        DESCRIPTION.
    v_atmbc : TYPE
        DESCRIPTION.

    '''

    
    grid = dict_PRD['grid']
    time_drying0 = dict_PRD['time_drying0']
    irr_days = dict_PRD['irr_days']
    flux = dict_PRD['flux']
    lg_PRD = dict_PRD['lg_PRD']

    
    x_min = max(grid["nodes_idxyz"][:, 1]) / 2
    x_max = max(grid["nodes_idxyz"][:, 1])
    y_min = min(grid["nodes_idxyz"][:, 2])
    y_max = max(grid["nodes_idxyz"][:, 2])

    nodes_xyz = grid["nodes_idxyz"]

    xmesh = grid["nodes_idxyz"][:, 1]
    ymesh = grid["nodes_idxyz"][:, 2]
    zmesh = grid["nodes_idxyz"][:, 3]

    HSPATM = 0
    IETO = 0

    # flux = 1.11111E-06
    totarea = 1  # 0.0

    sec_h = 3600
    days2sec = 24 * 3600



    time_drying0 = (time_drying0 * sec_h)  # days2sec # let the system dry out during 30days
    # irr_days = 1 # irrigate on left side during 10days
    # lg_PRD = 2 # 2 hours length of irrigation during the morning
    no_irr = 0
    t_irr = time_drying0  # initiate t_irr
    t_atmbc = [0]
    v_atmbc = [0]

    with open(os.path.join(workdir, project_name, "input/atmbc"), "w+") as atmbcfile:

        # write no irrigation during the first delta_drying0 days
        atmbcfile.write(
            str(HSPATM) + "\t" + str(IETO) + "\t" + "HSPATM" + "\t" + "IETO" + "\n"
        )
        atmbcfile.write(str(0) + "\t" + "TIME" + "\n")

        count_nmax = 0
        for i in range(int(grid["nnod3"])):
            if round(zmesh[i], 3) == round(max(zmesh), 3):
                count_nmax += 1
        print(count_nmax)

        # loop over the nodes of the mesh
        count = 0
        for i in range(int(grid["nnod3"])):
            if round(zmesh[i], 3) == round(max(zmesh), 3):
                # if i==max(range(int(grid['nnod3']))):
                if count == count_nmax:
                    # atmbcfile.write('0'.format('%f'))
                    atmbcfile.write(str(no_irr))
                else:
                    # atmbcfile.write('0'.format('%f')+ "\n")
                    atmbcfile.write(str(no_irr) + "\n")
            # atmbcfile.write('0'.format('%f')+ "\n")
            # atmbcfile.write(str(no_irr))

        # write irrigation once a day during 2h lasting 10 irr_days
        for k in range(irr_days):
            atmbcfile.write(str(t_irr) + "\t" + "TIME" + "\n")

            # loop over the nodes of the mesh
            count = 0
            for i in range(int(grid["nnod3"])):
                if round(zmesh[i], 3) == round(max(zmesh), 3):
                    count += 1
                    # loop over right or left side (odd or not)
                    if (k % 2) == 0:
                        if ((xmesh[i] > x_min) and (xmesh[i] < x_max)) and (
                            (ymesh[i] > y_min) and (ymesh[i] < y_max)
                        ):
                            # if count==count_nmax:
                            #     atmbcfile.write(str(flux/totarea))
                            # else:
                            atmbcfile.write(str(flux / totarea) + "\n")
                        else:
                            # if count==count_nmax:
                            #     atmbcfile.write(str(no_irr))
                            # else:
                            atmbcfile.write(str(no_irr) + "\n")

                    else:
                        if ((xmesh[i] > x_min) and (xmesh[i] < x_max)) and (
                            (ymesh[i] > y_min) and (ymesh[i] < y_max)
                        ):
                            # if count==count_nmax:
                            #     atmbcfile.write(str(no_irr))
                            # else:
                            atmbcfile.write(str(no_irr) + "\n")
                        else:
                            # if count==count_nmax:
                            #     atmbcfile.write(str(flux/totarea))
                            # else:
                            atmbcfile.write(str(flux / totarea) + "\n")

            t_atmbc.append(t_irr)
            t_irr += lg_PRD * sec_h
            v_atmbc.append(flux / totarea)

            # break irrigation during the rest of the day i.e. 24h -2h = 22h
            # atmbcfile.write("\n")
            atmbcfile.write(str(t_irr) + "\t" + "TIME" + "\n")
            t_atmbc.append(t_irr)
            # break during the resting hours
            count = 0
            for i in range(int(grid["nnod3"])):
                if round(zmesh[i], 3) == round(max(zmesh), 3):
                    count += 1
                    # if count==count_nmax:
                    #     atmbcfile.write(str(no_irr))
                    # else:
                    atmbcfile.write(str(no_irr) + "\n")

                    # atmbcfile.write(str(no_irr) +  "\n")

            v_atmbc.append(no_irr)
            t_irr += (24 - lg_PRD) * sec_h
            # t_atmbc.append(t_irr)
            # atmbcfile.write(str(no_irr) +  "\n")

        atmbcfile.close()

    if show == True:
        x_units = "sec"
        for key, value in kwargs.items():
            if key == "x_units":
                x_units = value
          
        pltCT.show_atmbc(
            t_atmbc, [v_atmbc, np.zeros(len(v_atmbc))], x_units=x_units)

    return t_irr, t_atmbc, v_atmbc
    
    
def ETp_synth(workdir,project_name,dict_ETp,show=False,**kwargs):
    
    '''
    https://plantmethods.biomedcentral.com/articles/10.1186/s13007-019-0409-9
    https://journals.ashs.org/hortsci/view/journals/hortsci/46/12/article-p1677.xml
    https://www.researchgate.net/publication/321734175_Root_growth_water_uptake_and_sap_flow_of_winter_wheat_in_response_to_different_soil_water_availability/figures?lo=1
    https://en.bio-protocol.org/e3190
    '''
    pass

    
def perturbate_rhizo(cathyDA,simu_DA,scenario,prj_name,NENS):
    
    
    list_pert = perturbate.perturbate(simu_DA, scenario, NENS)
        
    for dp in list_pert:
        
        np.random.seed(1)
        
        # need to call perturbate_var as many times as variable to perturbate
        # return a dict merging all variable perturbate to parse into prepare_DA
        parm_per = cathyDA.perturbate_parm(parm=dp, 
                                            type_parm = dp['type_parm'], # can also be VAN GENUCHTEN PARAMETERS
                                            mean =  dp['mean'],
                                            sd =  dp['sd'],
                                            sampling_type =  dp['sampling_type'],
                                            ensemble_size =  dp['ensemble_size'], # size of the ensemble
                                            per_type= dp['per_type'],
                                            show= dp['show'],
                                            savefig= os.path.join(prj_name,
                                                                  prj_name + dp['savefig'])
                                            )
        
    return parm_per
        
        

def update_rhizo_inputs_from_solution(simu_DA, nb_of_days,solution,**kwargs):
    '''
    Update rhizotron simulation input parameters 
    '''
     

    # # MESH 
    # # *****
    
    # dem_synth=np.ones([9,10])*1e-3
    # dem_synth[-1,-1]=0.99e-3
    # zb = np.arange(0,0.5+0.1,0.05)
    # nstr=len(zb)-1
    # zr=list((np.ones(len(zb)))/(nstr))
    # simu_DA.update_prepo_inputs(delta_x=0.05,delta_y=0.003, DEM=dem_synth,
    #                           N=np.shape(dem_synth)[1],
    #                           M=np.shape(dem_synth)[0],
    #                           N_celle=np.shape(dem_synth)[0]*np.shape(dem_synth)[1],
    #                           nstr=nstr, zratio=zr,  base=max(zb),dr=0.05,show=True)
    # simu_DA.run_preprocessor(verbose=True,
    #                          KeepOutlet=False) #,show=True)
    # simu_DA.update_parm()
    # simu_DA.update_veg_map()
    # simu_DA.update_soil()
    
    
    rhizotron.set_defaults(simu_DA)
    
    rhizotron.make_mesh(simu_DA)
         
    
    if 'tobs' in kwargs:
        time_of_interest=kwargs['tobs']       
    else:
        time_of_interest = list(np.arange(0,nb_of_days*3600,3600))
    
    
    # INITIAL CONDITIONS
    # ******************
    # Unbiased scenario
    # -------------------------------
    simu_DA.update_ic(INDP=solution['INDP'],
                      IPOND=solution['IPOND'],
                      WTPOSITION=solution['WTPOSITION'],
                      pressure_head_ini=solution['pressure_head_ini']
                      )
    
    # biased scenario
    # -------------------------------
    # ATMBC CONDITIONS
    # ******************
    simu_DA.update_atmbc(HSPATM=1,IETO=1,TIME=time_of_interest,
                      VALUE=np.ones(len(time_of_interest))*solution['forcing'],
                      show=True,x_units='hours') # just read new file and update parm and cathyH
    
    
    # simu_DA.update_atmbc(HSPATM=1,IETO=1,TIME=time_of_interest,
    #                   VALUE=[np.zeros(len(time_of_interest)),np.ones(len(time_of_interest))*flux],
    #                   show=True,x_units='hours',diff=True) # just read new file and update parm and cathyH
    

    """### 3- Boundary conditions"""
    
    simu_DA.update_nansfdirbc(noflow=True,TIME=time_of_interest) # Dirichlet uniform by default
    simu_DA.update_nansfneubc(TIME=time_of_interest) # Neumann (rain)
    simu_DA.update_sfbc(TIME=time_of_interest)
    
    
    """### 4- Soil and roots inputs"""
    
    # Unbiased scenario
    # -------------------------------
    PERMX = PERMY = PERMZ = solution['PERMX']
    ELSTOR = solution['ELSTOR']
    POROS = solution['POROS']
    VGNCELL = solution['VGNCELL']
    VGRMCCELL = solution['VGRMCCELL']
    VGPSATCELL =  solution['VGPSATCELL']
    
    # biased scenario
    # -------------------------------
    # if 'Ks' in kwargs:
    #     PERMX = PERMY = PERMZ = kwargs['Ks']

    
    SoilPhysProp = {'PERMX':PERMX,'PERMY':PERMY,'PERMZ':PERMZ,
                    'ELSTOR':ELSTOR,'POROS':POROS,
                    'VGNCELL':VGNCELL,'VGRMCCELL':VGRMCCELL,'VGPSATCELL':VGPSATCELL}
    
    """We could also used the predefined Van-Genuchten function
    (located in cathy_utils) 
    to generate a set of soil physical properties"""
    
    #SoilPhysProp = utils.VanG()
    
    """The trick to defining a heterogeneous density of roots 
    into the 3d medium from a single plant, 
    we defined different vegetation areas which are independent 
    with varying root depth and Feddes parameters."""
    
    PCANA=solution['PCANA']       
    PCREF=solution['PCREF']
    PCWLT=solution['PCWLT']
    ZROOT=solution['ZROOT']
    PZ = solution['PZ']
    OMGC = solution['OMGC']
    
    # if 'ZROOT' in kwargs:
    #     ZROOT=[kwargs['ZROOT']]
        
    # if 'PCREF' in kwargs:
    #     PCREF=[kwargs['PCREF']]     
        
    FeddesParam = {'PCANA':PCANA,'PCREF':PCREF,'PCWLT':PCWLT,
                   'ZROOT':ZROOT,'PZ':PZ,
                   'OMGC':OMGC}
    
    veg_map = np.ones([int(simu_DA.hapin['M']),int(simu_DA.hapin['N'])])
    
    if isinstance(solution['ZROOT'], float):
        solution['ZROOT'] = [solution['ZROOT']]
    
    if len(solution['ZROOT'])>1:
        veg_map[:,[0,9]]=1
        veg_map[:,[1,8]]=2
        veg_map[:,[2,7]]=3
        veg_map[:,[3,6]]=4
        veg_map[:,[4,5]]=5
        
        
    simu_DA.update_veg_map(indice_veg=veg_map, show=True)
    
    """The choice of PMIN conditionne the switching condition"""
    
    pmin = solution['pmin']
    
    
    
    soil_zone_map = np.ones([int(simu_DA.hapin['M']),int(simu_DA.hapin['N'])])
    
    if isinstance(PERMX, float):
        PERMX = [PERMX]
        
    if len(PERMX)>1:
        # soil_zone_map[list(np.arange(5,10)),:]=2
        # soil_zone_map[list(np.arange(5,9)),:]=2
        soil_zone_map=solution['soil_zone_map']
        simu_DA.update_zone(soil_zone_map)

    simu_DA.update_soil(  PMIN=pmin,
                          SPP=SoilPhysProp,
                          FP=FeddesParam,
                          verbose=True
                          
                        )
    
    simu_DA.update_parm(TIMPRTi=time_of_interest,
                        TMAX=time_of_interest[-1],
                        IPRT1=2)
    

    return simu_DA


#%% meshes interpolation

def plot_scatter_mesh_nodes(mesh_OUT,in_nodes_mod_m):
    

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(
                in_nodes_mod_m[:,0],
                in_nodes_mod_m[:,1],
                in_nodes_mod_m[:,2]
                )
    ax.scatter(
                np.array(mesh_OUT.points)[:,0],
                np.array(mesh_OUT.points)[:,1],
                np.array(mesh_OUT.points)[:,2]
                )
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()
    
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(
                in_nodes_mod_m[:,0],
                in_nodes_mod_m[:,1],
                in_nodes_mod_m[:,2]
                )
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(
                np.array(mesh_OUT.points)[:,0],
                np.array(mesh_OUT.points)[:,1],
                np.array(mesh_OUT.points)[:,2]
                )
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    
    pass

def CATHY_2_Simpeg(mesh_CATHY,mesh_Simpeg,scalar='saturation',show=False,**kwargs):
    pass

def CATHY_2_pg(mesh_CATHY,mesh_pg,scalar='saturation',show=False,**kwargs):
    '''
    THIS SHOULD BE MOVED TO MESHTOOLS
    
    Convert CATHY mesh attribute to pygimli
    Need to flip axis because convention for CATHY and pygimli are different

    Parameters
    ----------
    mesh_CATHY : TYPE
        DESCRIPTION.
    mesh_pg : TYPE
        DESCRIPTION.
    scalar : TYPE, optional
        DESCRIPTION. The default is 'saturation'.
    show : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    mesh_new_attr : TYPE
        DESCRIPTION.
    scalar_new : TYPE
        DESCRIPTION.

    '''
            
    
    
    mesh_IN = mesh_CATHY
    if type(mesh_CATHY) is str:
        mesh_IN_tmp = pv.read(mesh_CATHY)

    if type(mesh_pg) is str:
        mesh_OUT = pv.read(mesh_pg)
    
    
    # flip y and z axis as CATHY and pg have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_IN.points)
    # in_nodes_mod_pg = np.array(mesh_pg.points)
    idx = np.array([0, 2, 1])
    in_nodes_mod_m = in_nodes_mod[:, idx]
    in_nodes_mod_m = in_nodes_mod[:, idx]
    in_nodes_mod_m[:,2] = -np.flipud(in_nodes_mod_m[:,2])
    in_nodes_mod_m[:,1] = -np.flipud(in_nodes_mod_m[:,1])

    if 'dict_ERT' in kwargs:
        if 'mesh_nodes_modif' in kwargs['dict_ERT'].keys():
            in_nodes_mod_m = kwargs['dict_ERT']['mesh_nodes_modif']
        

    # plot_scatter_mesh_nodes(mesh_OUT,in_nodes_mod_m)
    
    
    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']


    data_OUT, warm_0 = mt.trace_mesh(mesh_IN,mesh_OUT,
                            scalar=scalar,
                            threshold=1e-1,
                            in_nodes_mod=in_nodes_mod_m)
    
    if len(warm_0)>0:
        print(warm_0)
        
    
    scalar_new = scalar + '_nearIntrp2_pg_msh' 
    if 'time' in kwargs:
        time = kwargs['time']
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_OUT, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                time=time,
                                                                path=path)
    else:
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_OUT, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                path=path)
    
    if show == True:

        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=False)
        p.add_mesh(mesh_new_attr,scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
        
        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_IN, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
            
    
    return mesh_new_attr, scalar_new


def CATHY_2_Resipy(mesh_CATHY,mesh_Resipy,scalar='saturation',show=False,**kwargs):
 

    # flip y and z axis as CATHY and Resipy have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    in_nodes_mod[:,2] = -np.flipud(in_nodes_mod[:,2])
    in_nodes_mod[:,1] = -np.flipud(in_nodes_mod[:,1])


    # check with a plot positon of the nodes for both meshes
    # ------------------------------------------------------------------------
    # p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
    # p.add_mesh(mesh_CATHY)
    # _ = p.add_points(np.array(mesh_CATHY.points), render_points_as_spheres=True,
    #                         color='red', point_size=20)
    # _ = p.add_points(in_nodes_mod, render_points_as_spheres=True,
    #                         color='blue', point_size=20)
    # _ = p.show_bounds(grid='front', all_edges=True, font_size=50)
    # cpos = p.show(True)
    
    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']

    data_OUT = mt.trace_mesh(mesh_CATHY,mesh_Resipy,
                            scalar=scalar,
                            threshold=1e-1,
                            in_nodes_mod=in_nodes_mod)
    
    # print(len(data_OUT))
    scalar_new = scalar + '_nearIntrp2Resipymsh'
    # print(mesh_Resipy)
    # get_array(mesh, name, preference='cell'
              

    if 'time' in kwargs:
        time = kwargs['time']
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_Resipy, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                time=time,
                                                                path=path)
    else:
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_Resipy, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                path=path)
    
    if show == True:

        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_new_attr,scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
    # if type(meshERT) is str:
    #     meshERTpv = pv.read(meshERT)
    
    # if savefig == True:
    
    #     plotter = pv.Plotter(notebook=True)
    #     _ = plotter.add_mesh(mesh_new_attr,show_edges=True)
    #     plotter.view_xz(negative=False)
    #     plotter.show_grid()
    #     plotter.save_graphic(path_CATHY + 'ERT' + str(DA_cnb) + str('.ps'), 
    #                           title='ERT'+ str(DA_cnb), 
    #                           raster=True, 
    #                           painter=True)
            
    
    return mesh_new_attr, scalar_new





# -------------------------------------------------------------------#
#%% Infitration DATA

def create_infitration(self, dirfiles):
    self.set_drippers(dirfiles)

    pass

def set_drippers(self, dirfiles, drip_pos="drippers.txt"):

    print(os.getcwd())
    if isinstance(drip_pos, str):
        self.drippers = np.loadtxt(
            os.path.join(self.project_name, dirfiles, drip_pos),
            skiprows=1,
            delimiter=",",
        )
    else:
        self.drippers = drip_pos

    # check drippers position against DEM
    self.hapin
    mesh_x_max = float(self.hapin["xllcorner"]) + float(
        self.hapin["delta_x"]
    ) * float(self.hapin["N"])
    mesh_y_max = float(self.hapin["yllcorner"]) + float(
        self.hapin["delta_y"]
    ) * float(self.hapin["M"])

    mesh_x = []
    for xx in range(int(self.hapin["N"])):
        mesh_x.append(
            float(self.hapin["xllcorner"]) + float(self.hapin["delta_x"]) * xx
        )

    mesh_y = []
    for yy in range(int(self.hapin["M"])):
        mesh_y.append(
            float(self.hapin["yllcorner"]) + float(self.hapin["delta_y"]) * yy
        )

    print(mesh_x)
    print(mesh_y)

    print(mesh_x_max)
    print(max(self.drippers[:, 0]))



    if mesh_x_max < max(self.drippers[:, 0]):
        print(
            "Error: max mesh_x="
            + str(mesh_x_max)
            + "; max dripper x pos="
            + str(max(self.drippers[:, 0]))
        )

    if mesh_y_max < max(self.drippers[:, 1]):
        print(
            "Error: max mesh_x="
            + str(mesh_y_max)
            + "; max dripper y pos="
            + str(max(self.drippers[:, 1]))
        )

    # for dd in self.drippers:
    #     dd==

    #  (A==B).all()

    self.drippers_nodes = []

    pass
