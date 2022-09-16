""" Managing Data Assimilation process pertubation. 
    Take a scenario dictionnary as input describing which and how parameters are perturbated
    Check consistency of the distribution
    Prepare for DA class
"""

def check_distribution(parm2check):
    if parm2check['type_parm']=='porosity':
        if parm2check['per_nom']<0:
            raise ValueError 
    
    if parm2check['type_parm']=='thetar_VG':
        if parm2check['per_nom']>1:
            print('! warming value to high for residual water VG')

    
    
    
    # References for Archie paarmeters range
    # https://www.sciencedirect.com/science/article/pii/S1018363918306123
    # a : TYPE, optional
    #     Tortuosity factor. The default is [1.0].
    # m : TYPE, optional
    #     Cementation exponent. The default is [2.0]. (usually in the range 1.3 -- 2.5 for sandstones)
    # n : TYPE, optional
    #     Saturation exponent. The default is [2.0].
            
    
def check4bounds(scenario, index, clip_min,clip_max,het_size=1,het_nb=None):
    
    try: # try since per_bounds is not mandatory
        if het_size==1:
            sc_bounds = scenario['per_bounds'][index]
        else:
            sc_bounds = scenario['per_bounds'][index][het_nb]
        clip_min = sc_bounds['min']
        clip_max = sc_bounds['max']
    except:
        pass
    
    return clip_min, clip_max
    

    
def perturbate(simu_DA,scenario,NENS):
    ''' Write a list of dictionaries, each containing all the informations on how to 
    perturbate the parameters based on the scenario to consider'''   
    
    list_pert = []

    #%% Initial and boundary conditions parameters
    # ------------------------------------------------------------------------
    if 'ic' in scenario['per_name']:      
        
        index = scenario['per_name'].index('ic')

        ic = {
              'type_parm': 'ic',
              'nominal':  scenario['per_nom'][index], #nominal value
              'mean':  scenario['per_mean'][index],
              'sd':  scenario['per_sigma'][index],
              'units': 'pressure head $(m)$', # units
              'sampling_type': 'normal',
              'ensemble_size': NENS, # size of the ensemble
              'per_type': scenario['per_type'][index],
              'savefig': 'ic.png',
              }    
        list_pert.append(ic)     
        
    if 'atmbc' in scenario['per_name']:
        index = scenario['per_name'].index('atmbc')

        atmbc_nom = scenario['per_nom'] #1e-4
        atmbc_mean = scenario['per_mean'] # sampling mean
        atmbc_sd = scenario['per_sigma']  #1.53
        
        atmbc = {'nominal': atmbc_nom, #nominal value 45
                'units': '$mm.h^{-1}$',
                'type_parm': 'atmbc',
                'nominal': atmbc_nom, #nominal value
                'mean':atmbc_nom,
                'sd':atmbc_sd,
                'atmbc_units': 'atmospheric forcing $mm.h^{-1}$', # units
                'sampling_type': 'lognormal',
                'ensemble_size':NENS, # size of the ensemble
                'per_type': 'multiplicative',
                'time_variable': True,
                # 'data2assimilate': kwargs['data'],
                'data2assimilate':simu_DA.atmbc,
                'time_decorrelation_len': 277200*10,
                'savefig': 'atmbc.png',
                }
        
        list_pert.append(atmbc)
    
    
    #%% Petro-Pedophysical parameters 
    # ------------------------------------------------------------------------
    if 'thetar_VG' in scenario['per_name']:   
        index = scenario['per_name'].index('thetar_VG')
        
        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]
        
        
        clip_min = 0
        clip_max = 1
        
        for nz in range(len(simu_DA.soil_SPP['SPP_map']['PERMX'])):
            if len(simu_DA.soil_SPP['SPP_map']['PERMX'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]
                
            clip_min,clip_max = check4bounds(scenario,index,clip_min,clip_max,
                                             het_size=len(simu_DA.soil_SPP['SPP_map']['PERMX']),
                                             het_nb=nz,
                                             )

            thetar_VG = {
                  'type_parm': 'thetar_VG'+ str(nz),
                  'nominal':  scenario_nom, #nominal value
                  'mean':  scenario_mean,
                  'sd':  scenario_sd,
                  'units': '$m^{3}/m^{3}$', # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'thetar_VG'+ str(nz) + '.png',
                  'surf_zones_param': nz,
                  'clip_min':clip_min,
                  'clip_max':clip_max,
                  }    
            list_pert.append(thetar_VG)
            
           
    if 'alpha_VG' in scenario['per_name']:
        index = scenario['per_name'].index('alpha_VG')
        
        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]
        
        for nz in range(len(simu_DA.soil_SPP['SPP_map']['PERMX'])):
            if len(simu_DA.soil_SPP['SPP_map']['PERMX'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]

            alpha_VG = {
                  'type_parm': 'alpha_VG'+ str(nz),
                  'nominal':  scenario_nom, #nominal value
                  'mean':  scenario_mean,
                  'sd':  scenario_sd,
                  'units': '$cm^{-1}$', # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'alpha_VG'+ str(nz) + '.png',
                  'surf_zones_param': nz
                  }    
            list_pert.append(alpha_VG)

    if 'VGPSATCELL' in scenario['per_name']:
        index = scenario['per_name'].index('VGPSATCELL')
        
        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]
        
        
        # clip_min = 0
        # clip_max = 1
        
        for nz in range(len(simu_DA.soil_SPP['SPP_map']['PERMX'])):
            if len(simu_DA.soil_SPP['SPP_map']['PERMX'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]

            VGPSATCELL_VG = {
                  'type_parm': 'VGPSATCELL_VG'+ str(nz),
                  'nominal':  scenario_nom, #nominal value
                  'mean':  scenario_mean,
                  'sd':  scenario_sd,
                  'units': '$cm^{-1}$', # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'VGPSATCELL_VG'+ str(nz) + '.png',
                  'surf_zones_param': nz
                  }    
            list_pert.append(VGPSATCELL_VG)
                

    if 'n_VG' in scenario['per_name']:   
        index = scenario['per_name'].index('n_VG')
        
        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]
        
        for nz in range(len(simu_DA.soil_SPP['SPP_map']['PERMX'])):
            if len(simu_DA.soil_SPP['SPP_map']['PERMX'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]

            n_VG = {
                  'type_parm': 'n_VG'+ str(nz),
                  'nominal':  scenario_nom, #nominal value
                  'mean':  scenario_mean,
                  'sd':  scenario_sd,
                  'units': '$(-)$', # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'n_VG'+ str(nz) + '.png',
                  'surf_zones_param': nz
                  }    
            list_pert.append(n_VG)
        
    
    if 'VGP' in scenario['per_name']:   
        # - 'PERMX' (NSTR, NZONE): saturated hydraulic conductivity - xx
        # - 'PERMY' (NSTR, NZONE): saturated hydraulic conductivity - yy
        # - 'PERMZ' (NSTR, NZONE): saturated hydraulic conductivity - zz
        # - 'ELSTOR' (NSTR, NZONE): specific storage
        # - 'POROS'  (NSTR, NZONE): porosity (moisture content at saturation) = \thetaS

        # retention curves parameters VGN, VGRMC, and VGPSAT
        # - 'VGNCELL' (NSTR, NZONE): van Genuchten curve exponent  = n
        # - 'VGRMCCELL' (NSTR, NZONE): residual moisture content = \thetaR
        # - 'VGPSATCELL' (NSTR, NZONE): van Genuchten curve exponent --> 
        #                               VGPSAT == -1/alpha (with alpha expressed in [L-1]);
        # ['ks','ss','phi','thetar','alpha','n'])
        # ['PERMX','ELSTOR','POROS','VGRMCCELL','VGPSATCELL','VGNCELL'])
    
    
        # need to account for transformation of the parameters see Boto
        
        VGP_units = ['','','','','','']
        
        for i, p in enumerate(['ks','ss','phi','thetar','alpha','r']):
            index = scenario['per_name'].index('VGP')
            
            if p not in scenario['per_nom'][index]:
                pass
            else:
                p_VGP = {
                      'type_parm': p,
                      'nominal':  scenario['per_nom'][index][p], #nominal value
                      'mean':  scenario['per_mean'][index][p],
                      'sd':  scenario['per_sigma'][index][p],
                      'units': VGP_units[i], # units
                      'sampling_type': 'normal',
                      'ensemble_size': NENS, # size of the ensemble
                      'per_type': scenario['per_type'][index],
                      'savefig': '_VGP' + p + '.png',
                      }    
                list_pert.append(p_VGP)        

        
    Archie_2pert = [ele for ele in scenario['per_name'] if('Archie' in ele)]
    if len(Archie_2pert)>0: 
        # Initialement, chacun des membres de l’ensemble aurait des paramètres d’Archie
        # légèrement différents, qui seraient utilisés dans la fonction d’observation afin de simuler les
        # données à partir de l’état a priori. Ces paramètres seraient ensuite modifiés à l’étape de
        # mise à jour, de la même façon que les variables de teneur en eau.
        # rFluid=[1.0],a=[1.0],m=[2.0],n=[2.0]
        # simu_DA.Archie_parms
        # literature example:
        # ---------------------------------------------------------------------
        # https://www.sciencedirect.com/science/article/pii/S0169772220302680#tf0005
        # c = (mean = 1.6,sd = 0.5,min = 0.0, max = 2.0),
        # m = (mean = 2.5,sd = 0.8,min = 0.0, max = 3.0)
        
    
        archie_units = ['','','','']
        
        for i, p in enumerate(Archie_2pert):
            index = scenario['per_name'].index(p)
            
            p_archie = {
                  'type_parm': p,
                  'nominal':  scenario['per_nom'][index], #nominal value
                  'mean':  scenario['per_mean'][index],
                  'sd':  scenario['per_sigma'][index],
                  'units': archie_units[i], # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': '_Archie' + p + '.png',
                  }    
            list_pert.append(p_archie)
        
        
    #%% SOIL parameters
    # ------------------------------------------------------------------------

    if 'ks' in scenario['per_name']:
        index = scenario['per_name'].index('ks')
        

        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]
        
        for nz in range(len(simu_DA.soil_SPP['SPP_map']['PERMX'])):
            if len(simu_DA.soil_SPP['SPP_map']['PERMX'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]

            ks = {
                  'type_parm': 'ks'+ str(nz),
                  'nominal':  scenario_nom, #nominal value
                  'mean':  scenario_mean,
                  'sd':  scenario_sd,
                  'units': '$m.s^{-1}$', # units
                  'sampling_type': 'lognormal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'ks'+ str(nz) + '.png',
                  'surf_zones_param': nz
                  }    
            list_pert.append(ks)
            

    if 'porosity' in scenario['per_name']:
        index = scenario['per_name'].index('porosity')
        

        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]
        
        clip_min = 0.2 # minimum soil porosity
        clip_max = 0.7 # maximum soil porosity
                
        for nz in range(len(simu_DA.soil_SPP['SPP_map']['POROS'])):
            print('zone nb:' + str(nz))
            
            clip_min,clip_max = check4bounds(scenario,index,clip_min,clip_max,
                                             het_size=len(simu_DA.soil_SPP['SPP_map']['POROS']),
                                             het_nb=nz,
                                             )
            
            if len(simu_DA.soil_SPP['SPP_map']['POROS'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]
                

            porosity = {
                  'type_parm': 'porosity'+ str(nz),
                  'nominal':  scenario_nom, #nominal value
                  'mean':  scenario_mean,
                  'sd':  scenario_sd,
                  'units': '', # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'porosity'+ str(nz) + '.png',
                  'surf_zones_param': nz,
                  'clip_min':clip_min,
                  'clip_max':clip_max,
                  } 
            
            check_distribution(porosity)

            list_pert.append(porosity)

    #%% Plant parameters
    # ------------------------------------------------------------------------

        
    if 'PCREF' in scenario['per_name']:
        index = scenario['per_name'].index('PCREF')
        
        for nz in range(len(simu_DA.soil['PCREF'])):

            PCREF = {
                  'type_parm': 'PCREF'+ str(nz),
                  'nominal':  scenario['per_nom'][index], #nominal value
                  'mean':  scenario['per_mean'][index],
                  'sd':  scenario['per_sigma'][index],
                  'units': '$m$', # units
                  'sampling_type': 'normal',
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'PCREF.png',
                  'surf_zones_param': nz
                  # 'clip_min':scenario['per_bounds'][index]['min'],
                  # 'clip_max':scenario['per_bounds'][index]['max'],
                  }    
            list_pert.append(PCREF)
                   
        
        
    if 'ZROOT' in scenario['per_name']:
        index = scenario['per_name'].index('ZROOT')


        scenario_nom = scenario['per_nom'][index]
        scenario_mean = scenario['per_mean'][index]
        scenario_sd = scenario['per_sigma'][index]

        if 'sampling_type' in scenario:
            scenario_sampling = scenario['sampling_type'][index]
        else:
            scenario_sampling = 'lognormal'

        clip_min = 0
        clip_max = None


        for nz in range(len(simu_DA.soil['ZROOT'])):
            
            clip_min,clip_max = check4bounds(scenario,index,clip_min,clip_max,
                                             het_size=len(simu_DA.soil['ZROOT']),
                                             het_nb=nz,
                                             )
            
            if len(simu_DA.soil['ZROOT'])>1:
                scenario_nom = scenario['per_nom'][index][nz]
                scenario_mean = scenario['per_mean'][index][nz]
                scenario_sd = scenario['per_sigma'][index][nz]
            
            ZROOT = {
                  'type_parm': 'ZROOT'+ str(nz),
                  'nominal': scenario_nom, #nominal value
                  'mean': scenario_mean,
                  'sd':  scenario_sd,
                  'units': '$m.s^{-1}$', # units
                  'sampling_type': scenario_sampling,
                  'ensemble_size': NENS, # size of the ensemble
                  'per_type': scenario['per_type'][index],
                  'savefig': 'ZROOT' + str(nz) + '.png',
                  'surf_zones_param': nz,
                  'clip_min':clip_min,
                  'clip_max':clip_max,
                      
                  }    
            list_pert.append(ZROOT)
            # nb_surf_nodes = 110


        
    return list_pert


