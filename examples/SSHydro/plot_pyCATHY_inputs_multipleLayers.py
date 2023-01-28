"""
Input plots
===========

This example shows how to use pyCATHY object to plot inputs of the hydrological model.

*Estimated time to run the notebook = 5min*

"""

# map_prop_veg ?
# map_prop2zone



#%% Import packages

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import numpy as np

#%% add CATHY object

path2prj ='weil_exemple_inputs_multipleLayers' # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj)


#%%
simu.update_parm()
simu.update_zone()
simu.update_soil()
df_soil, _ = simu.read_inputs('soil')
df = simu.read_inputs('soil')
simu.show_input(prop='soil', yprop='PERMX')

#%% Add layers

import matplotlib.pyplot as plt
import pandas as pd


def get_layer_depths(simu):
    for li in range(simu.dem_parameters["nstr"]):
        get_layer_depth(simu,li)
    
def get_layer_depth(simu,li):
    
    dempar = simu.dem_parameters['zratio(i),i=1,nstr'].split('\t')
    dempar_ratio = [float(d) for d in dempar]
    
    # layeri_top = np.round(dempar_ratio[li]*simu.dem_parameters["base"]/simu.dem_parameters['nstr'],2)
    # layeri_top = dempar_ratio[li]*simu.dem_parameters["base"]/simu.dem_parameters['nstr']
    if li==0:
        layeri_top = 0
    else:
        # layeri_top = np.cumsum(dempar_ratio[0:li])[-1]*(simu.dem_parameters["base"]/(simu.dem_parameters['nstr']-1))
        layeri_top = np.cumsum(dempar_ratio[0:li])[-1]*(simu.dem_parameters["base"])

        # layeri_top = np.cumsum(dempar_ratio[0:li])[-1] 

    if (li+1)<len(dempar_ratio):
        # layeri_bottom = np.round(dempar_ratio[li+1]*simu.dem_parameters["base"]/simu.dem_parameters['nstr'],2)
        # layeri_bottom = layeri_top + dempar_ratio[li+1]*simu.dem_parameters["base"]/(simu.dem_parameters['nstr']-1)
        # layeri_bottom = np.cumsum(dempar_ratio[0:li+1])[-1]*(simu.dem_parameters["base"]/(simu.dem_parameters['nstr']-1))
        layeri_bottom = np.cumsum(dempar_ratio[0:li+1])[-1]*(simu.dem_parameters["base"])
    else:
        layeri_bottom = simu.dem_parameters["base"]

    print(layeri_top,layeri_bottom)
    return layeri_top, layeri_bottom

def zone3d(simu):
    
    # define zone in 3dimension - duplicate number of layer (nstr) times
    # ---------------------------------------------------------------
    zones3d = [simu.zone]*simu.dem_parameters["nstr"]
    
    fig, axs = plt.subplots(int(simu.dem_parameters["nstr"]/2)+1,2, sharex=True,sharey=(True),
                            constrained_layout=False)
    plt.tight_layout()
    axs=axs.ravel()
    for li in range(simu.dem_parameters["nstr"]):
        layeri_top, layeri_bottom = get_layer_depth(simu,li)
        layer_str = '[' + str(f'{layeri_top:.2E}') + '-' + str(f'{layeri_bottom:.2E}') + ']'
        
        
        pmesh = cplt.show_raster(zones3d[li], prop=layer_str, cmap='jet',
                         ax=axs[li],vmin=0,vmax=1)
    
    plt.colorbar(pmesh, ax=axs[:], location='right', shrink=0.6, cmap='jet')
    plt.tight_layout()

    return zones3d


def create_layers_inzones3d(simu,zones3d,layers_names,layers_depths=[[0,1e99]]):
    
    # Loop over layers and zones and change flag if depth conditions is not respected
    # ----------------------------------------------------------------------------
    zones3d_layered = np.ones(np.shape(zones3d))
    for li in range(simu.dem_parameters["nstr"]):
        
        layeri_top, layeri_bottom = get_layer_depth(simu,li)
        # layer_str = '[' + str(layeri_top) + '-' + str(layeri_bottom) + ']'
        
        for zi in range(len(layers_names)):
            print('layers %i analysis' %zi)
            # print(layers_depths[zi][0])
            # print(layeri_top)
            # print(layers_depths[zi][1])
            # print(layeri_bottom)
            zi = 0
            
            # print(layers_depths[zi][0]<=layeri_top)
            # if depth of zone i is sup to mesh layers depth
            #---------------------------------------------------------------------
            # if (layers_depths[zi][0]>=layeri_top) & (layers_depths[zi][1]<layeri_bottom):     
            if (layers_depths[zi][0]<=layeri_top) & (layers_depths[zi][1]>=layeri_bottom):     
            # if (depths_ordered[zi][1]>=layeri_bottom):  
                print('conds ok --> zi:' + str(zi))
                print(layers_depths[zi])
                print(layeri_top,layeri_bottom )

                zones3d_layered[li][zones3d[li]==zi+1]=zi+2
                print('replacing ' + str(np.sum(zones3d[li]==zi+1)) + ' values by' + str(zi+1))
                
                if 10.5*layers_depths[zi][1]<layeri_bottom:
                    raise ValueError('Required layer is finer than mesh layers - refine mesh')
            else:
                print('conds not ok')
                print(layers_depths[zi])
                print(layeri_top,layeri_bottom )
                
    return zones3d_layered

        
def plot_zones3d_layered(simu,zones3d_layered):
    
    fig, axs = plt.subplots(int(simu.dem_parameters["nstr"]/2)+1,2, sharex=False,sharey=(False),
                            constrained_layout=False)
    axs=axs.flat
    for li in range(simu.dem_parameters["nstr"]):
        
        layeri_top, layeri_bottom = get_layer_depth(simu,li)
        layer_str = '[' + str(layeri_top) + '-' + str(layeri_bottom) + ']'
        
        pmesh = cplt.show_raster(zones3d_layered[li], prop=layer_str, #, cmap='jet',
                         ax=axs[li], vmin=0, vmax=2)
        
    plt.colorbar(pmesh, ax=axs[:], location='right', shrink=0.6, cmap='jet')
    plt.tight_layout()


def het_soil_layers_mapping_generic(simu,propertie_names,SPP,layers_names,layers_depths):
   
    # extend to 3d the zone raster file
    # -----------------------------------
    zones3d = zone3d(simu)

    # insert layers flag into the 3d the zone raster file
    # -----------------------------------
    zones3d_layered = create_layers_inzones3d(simu,zones3d,layers_names,layers_depths)
            
    # plot 3d zones files layered
    # ------------------------------------------
    plot_zones3d_layered(simu,zones3d_layered)
    
    # np.shape(zones3d_layered)
    # np.shape(zones3d_axis_swap)


    # Loop over soil properties names
    # -------------------------------------------------
    index_raster = np.arange(0,simu.hapin['N']*simu.hapin['M'])
    zones3d_axis_swap = np.swapaxes(zones3d_layered,0,2)
    prop_df = [] # properties dataframe
    layers_id = [ 'L' + str(i+1) for i in range(simu.dem_parameters["nstr"])]
    layers_id = np.flipud(layers_id)
    
    for i, p in enumerate(propertie_names):
        
        p_df = np.zeros(np.shape(zones3d_axis_swap))
          
        # Loop over soil layers and assign value of soil properties
        # -------------------------------------------------
        for k, lname in enumerate(layers_names):
            p_df[zones3d_axis_swap == k+1] = SPP[k][i]
                        
        df_tmp = pd.DataFrame(
                                np.vstack(p_df),
                                columns= layers_id,
                                index= index_raster,
                                )
        
        prop_df.append(df_tmp)
            
    SoilPhysProp_df_het_layers_p = pd.concat(prop_df, axis=1, keys=propertie_names)
    SoilPhysProp_df_het_layers_p.index.name = 'id raster'
    SoilPhysProp_df_het_layers_p.columns.names = ['soilp','layerid']
        
    
    SPP_map_dict = {}
    for p in propertie_names:
        SPP_map_dict[p] = []
        for li in range(simu.dem_parameters["nstr"]):
            v = SoilPhysProp_df_het_layers_p.iloc[0][p].loc['L'+str(li+1)]
            SPP_map_dict[p].append(v)
            
            
    return SoilPhysProp_df_het_layers_p, SPP_map_dict

#%%

get_layer_depths(simu)


layers_names = ['sub','bedrock']
layers_depths = [[0,1],[1,3]]
propertie_names = simu.soil_SPP['SPP_map'].keys()
p_subsurface_default = [1e-6,1e-6,1e-6,1.00E-05,0.48,1.914e+00,1.296e-01,1.240e+00]

SPP_perlayer = []
for li in range(len(layers_names)):
    SPP_p = []
    for i, p in enumerate(propertie_names):
        if li==1:
            SPP_p.append(p_subsurface_default[i]+1)
        else:
            SPP_p.append(p_subsurface_default[i])
    SPP_perlayer.append(SPP_p)


#%%
SoilPhysProp, SPP_map_dict = het_soil_layers_mapping_generic(simu,
                                                            propertie_names,
                                                            SPP_perlayer,
                                                            layers_names,
                                                            layers_depths,
                                                            )

#%%

simu.update_soil(IVGHU=0,
                  SPP_map=SPP_map_dict,
                  soil_heteregeneous_in_z=True,
                  # FP_map=FP_map,
                  )       
    


#%%
simu.show_input(prop='soil', yprop='PERMX', layer_nb=14)

        