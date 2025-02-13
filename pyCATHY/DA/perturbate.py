""" Managing Data Assimilation process pertubation. 
    Take a scenario dictionnary as input describing which and how parameters are perturbated
    Check consistency of the distribution
    Prepare for DA class
"""

import pandas as pd
import numpy as np
from scipy.stats import stats, qmc, norm
from scipy.stats import truncnorm
import matplotlib.pyplot as plt

from pyCATHY import cathy_utils as utils_CT
from pyCATHY.plotters import cathy_plots as plt_CT
import re 
import os

#%%
def check_distribution(parm2check):
    if parm2check["type_parm"] == "porosity":
        if parm2check["nominal"] < 0:
            raise ValueError

    if parm2check["type_parm"] == "thetar_VG":
        if parm2check["nominal"] > 1:
            print("! warming value to high for residual water VG")

    # References for Archie paarmeters range
    # https://www.sciencedirect.com/science/article/pii/S1018363918306123
    # a : TYPE, optional
    #     Tortuosity factor. The default is [1.0].
    # m : TYPE, optional
    #     Cementation exponent. The default is [2.0]. (usually in the range 1.3 -- 2.5 for sandstones)
    # n : TYPE, optional
    #     Saturation exponent. The default is [2.0].


def check4bounds(scenario, index, clip_min, clip_max, het_size=1, het_nb=None):

    try:  # try since per_bounds is not mandatory
        if het_size == 1:
            sc_bounds = scenario["per_bounds"][index]
        else:
            sc_bounds = scenario["per_bounds"][index][het_nb]
        clip_min = sc_bounds["min"]
        clip_max = sc_bounds["max"]
    except:
        pass

    return clip_min, clip_max

def create_dict_entry(
                         name,scenario,
                         index,nstri, 
                         scenario_nom, scenario_mean, scenario_sd, 
                         clip_min,clip_max,
                         NENS,
                         pert_control_name=None
                      ):
    
    per_type = scenario["per_type"][index]
    # if type(nstri) == int: 
    #     per_type = scenario["per_type"][index][nstri]
        
    dict_parm_entry = {
        "type_parm": name + str(nstri),
        "nominal": scenario_nom,  # nominal value
        "mean": scenario_mean,
        "sd": scenario_sd,
        "units": "",  # units
        "sampling_type": "normal",
        "ensemble_size": NENS,  # size of the ensemble
        "per_type": per_type,
        "savefig": name + str(nstri) + ".png",
        "surf_zones_param": nstri,
        "clip_min": clip_min,
        "clip_max": clip_max,
        "pert_control_name": pert_control_name,
    }
    
    return dict_parm_entry


def perturbate(simu_DA, scenario, NENS, 
               **kwargs):
    """Write a list of dictionaries, each containing all the informations on how to
    perturbate the parameters based on the scenario to consider
    
    
    pertControl = 'Layer'        
    Perturbation per zone is not yet implemented -  Assuming that the 
    dictionnary of perturbated parameters is build per layers i.e. 
    ks0= layer 0, ks1=layer 1, etc...
   
    """    
    if 'nzones' in kwargs: 
        nzones = kwargs.pop('nzones')
    else:
        nzones = len(simu_DA.soil_SPP['SPP_map'].index.get_level_values(0).unique())

    if 'nlayers' in kwargs: 
        nzones = kwargs.pop('nlayers')
    else:
        nlayers = len(simu_DA.soil_SPP['SPP_map'].index.get_level_values(1).unique())
        nlayers_PERMX = len(simu_DA.soil_SPP['SPP_map']['PERMX'].unique())
        
    list_pert = []

    #%% Initial and boundary conditions parameters
    # ------------------------------------------------------------------------
    
    if "WTPOSITION" in scenario["per_name"]:

        index = scenario["per_name"].index("WTPOSITION")

        clip_min = 0
        clip_max = simu_DA.dem_parameters['base']
        
        clip_min, clip_max = check4bounds(
            scenario,
            index,
            clip_min,
            clip_max,
        )
            
        WTPOSITION = {
            "type_parm": "WTPOSITION",
            "nominal": scenario["per_nom"][index],  # nominal value
            "mean": scenario["per_mean"][index],
            "sd": scenario["per_sigma"][index],
            "units": "Water table position $(m)$",  # units
            "sampling_type": "normal",
            "ensemble_size": NENS,  # size of the ensemble
            "per_type": scenario["per_type"][index],
            "savefig": "WTPOSITION.png",
            "clip_min": clip_min,
            "clip_max": clip_max,
        }
        list_pert.append(WTPOSITION)
        
        
    if "ic" in scenario["per_name"]:

        index = scenario["per_name"].index("ic")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        
        clip_min = None
        clip_max = None
        if type(scenario["per_nom"][index]) is list:                   
            for nstri in range(nlayers):
                
                if nlayers > 1:
                    scenario_nom = scenario["per_nom"][index][nstri]
                    scenario_mean = scenario["per_mean"][index][nstri]
                    scenario_sd = scenario["per_sigma"][index][nstri]
    
                ic = create_dict_entry(
                            'ic',scenario,index,nstri, 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            )
                
                list_pert.append(ic)
        else: 
            
            ic = create_dict_entry(
                        'ic',scenario,index,'', 
                        scenario_nom, scenario_mean, scenario_sd, 
                        clip_min,clip_max,
                        NENS,
                        )
            list_pert.append(ic)

    if "atmbc" in scenario["per_name"]:
        index = scenario["per_name"].index("atmbc")

        atmbc_nom = scenario["per_nom"]  # 1e-4
        atmbc_mean = scenario["per_mean"]  # sampling mean
        atmbc_sd = scenario["per_sigma"]  # 1.53
        

        df = simu_DA.atmbc['atmbc_df']
        # Get unique 'time' values
        unique_times = df['time'].unique()

        # type_parm
        # extract value of each node of the mesh for all the atmbc times
        for value_idx in range(len(df) // len(unique_times)):
            # print(value_idx)
            values = []
            for time in unique_times:
                # print(time)
                value = df.loc[df['time'] == time, 'value'].iloc[value_idx]
                values.append(value)
            # print(values)
            atmbc = {
                "type_parm": "atmbc" + str(value_idx),
                "nominal": atmbc_nom[index], 
                "units": "$mm.h^{-1}$",
                "mean": atmbc_nom[index],
                "sd": atmbc_sd[index],
                "atmbc_units": "atmospheric forcing $mm.h^{-1}$",
                "sampling_type": 'normal',
                "ensemble_size": NENS,
                "per_type": scenario["per_type"][index],
                "time_variable": True,
                "data2perturbate": {'time': simu_DA.atmbc['time'], 'VALUE': values},
                "time_decorrelation_len": scenario["time_decorrelation_len"][index],
                "savefig": "atmbc.png"
            }
            
            list_pert.append(atmbc)

    #%% Petro-Pedophysical parameters
    # ------------------------------------------------------------------------
    if "thetar_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("thetar_VG")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        clip_min = 0
        clip_max = 1

        for nstri in range(nlayers):
            if nlayers > 1:
                scenario_nom = scenario["per_nom"][index][nstri]
                scenario_mean = scenario["per_mean"][index][nstri]
                scenario_sd = scenario["per_sigma"][index][nstri]

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil_SPP["SPP_map"]["PERMX"]),
                het_nb=nstri,
            )

            thetar_VG = {
                "type_parm": "thetar_VG" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$m^{3}/m^{3}$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "thetar_VG" + str(nstri) + ".png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(thetar_VG)

    if "alpha_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("alpha_VG")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        for nstri in range(nlayers):
            if nlayers > 1:
                scenario_nom = scenario["per_nom"][index][nstri]
                scenario_mean = scenario["per_mean"][index][nstri]
                scenario_sd = scenario["per_sigma"][index][nstri]

            alpha_VG = {
                "type_parm": "alpha_VG" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$cm^{-1}$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "alpha_VG" + str(nstri) + ".png",
                "surf_zones_param": nstri,
            }
            list_pert.append(alpha_VG)

    if "VGPSATCELL" in scenario["per_name"]:
        index = scenario["per_name"].index("VGPSATCELL")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        # clip_min = 0
        # clip_max = 1

        for nstri in range(nlayers):
            if nlayers > 1:
                scenario_nom = scenario["per_nom"][index][nstri]
                scenario_mean = scenario["per_mean"][index][nstri]
                scenario_sd = scenario["per_sigma"][index][nstri]

            VGPSATCELL_VG = {
                "type_parm": "VGPSATCELL_VG" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$cm^{-1}$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "VGPSATCELL_VG" + str(nstri) + ".png",
                "surf_zones_param": nstri,
            }
            list_pert.append(VGPSATCELL_VG)

    if "n_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("n_VG")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        for nstri in range(nlayers):
            if nlayers > 1:
                scenario_nom = scenario["per_nom"][index][nstri]
                scenario_mean = scenario["per_mean"][index][nstri]
                scenario_sd = scenario["per_sigma"][index][nstri]

            n_VG = {
                "type_parm": "n_VG" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$(-)$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "n_VG" + str(nstri) + ".png",
                "surf_zones_param": nstri,
            }
            list_pert.append(n_VG)

    if "VGP" in scenario["per_name"]:
        # - 'PERMX' (NSTR, nstriONE): saturated hydraulic conductivity - xx
        # - 'PERMY' (NSTR, nstriONE): saturated hydraulic conductivity - yy
        # - 'PERMZ' (NSTR, nstriONE): saturated hydraulic conductivity - zz
        # - 'ELSTOR' (NSTR, nstriONE): specific storage
        # - 'POROS'  (NSTR, nstriONE): porosity (moisture content at saturation) = \thetaS

        # retention curves parameters VGN, VGRMC, and VGPSAT
        # - 'VGNCELL' (NSTR, nstriONE): van Genuchten curve exponent  = n
        # - 'VGRMCCELL' (NSTR, nstriONE): residual moisture content = \thetaR
        # - 'VGPSATCELL' (NSTR, nstriONE): van Genuchten curve exponent -->
        #                               VGPSAT == -1/alpha (with alpha expressed in [L-1]);
        # ['ks','ss','phi','thetar','alpha','n'])
        # ['PERMX','ELSTOR','POROS','VGRMCCELL','VGPSATCELL','VGNCELL'])

        # need to account for transformation of the parameters see Boto

        VGP_units = ["", "", "", "", "", ""]

        for i, p in enumerate(["ks", "ss", "phi", "thetar", "alpha", "r"]):
            index = scenario["per_name"].index("VGP")

            if p not in scenario["per_nom"][index]:
                pass
            else:
                p_VGP = {
                    "type_parm": p,
                    "nominal": scenario["per_nom"][index][p],  # nominal value
                    "mean": scenario["per_mean"][index][p],
                    "sd": scenario["per_sigma"][index][p],
                    "units": VGP_units[i],  # units
                    "sampling_type": "normal",
                    "ensemble_size": NENS,  # size of the ensemble
                    "per_type": scenario["per_type"][index],
                    "savefig": "_VGP" + p + ".png",
                }
                list_pert.append(p_VGP)

    Archie_2pert = [ele for ele in scenario["per_name"] if ("Archie" in ele)]
    if len(Archie_2pert) > 0:
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

        archie_units = ["", "", "", ""]

        for i, p in enumerate(Archie_2pert):
            index = scenario["per_name"].index(p)

            p_archie = {
                "type_parm": p,
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": archie_units[i],  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "_Archie" + p + ".png",
            }
            list_pert.append(p_archie)

    #%% SOIL parameters
    # ------------------------------------------------------------------------


    if "Ks" in scenario["per_name"]:
        index = scenario["per_name"].index("Ks")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        scenario_sampling = "normal"
        if "sampling_type" in scenario:
            scenario_sampling = scenario["sampling_type"][index]
            
        clip_min = None
        clip_max = None

        # Handling different configurations for zones, root maps, and layer perturbations of Ks:
        # If the scenario includes ['zone', 'root_map', 'layers'], it should contain a list of lists of lists.
        # If it includes ['zone', 'root_map'], it should contain a list of lists.
        # And so on...

        pert_control_name = None
        if 'pert_control_name' in scenario:
            pert_control_name = scenario['pert_control_name'][index]
        
        df_SPP, df_FP = simu_DA.read_inputs('soil')
        zone_raster, zone_header = simu_DA.read_inputs('zone')
        veg_raster, veg_header = simu_DA.read_inputs('root_map')
        
        if pert_control_name=='zone': 
            for zonei in range(len(np.unique(zone_raster))):
           
                if len(np.unique(zone_raster)) > 1:
                    scenario_nom = scenario["per_nom"][index][zonei]
                    scenario_mean = scenario["per_mean"][index][zonei]
                    scenario_sd = scenario["per_sigma"][index][zonei]
                Ks = create_dict_entry(
                            'Ks',scenario,index,zonei, 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            pert_control_name=pert_control_name
                            )
                
                list_pert.append(Ks)
        
        elif pert_control_name=='root_map': 
            for vegi in range(len(np.unique(veg_raster))):
           
                if len(np.unique(zone_raster)) > 1:
                    scenario_nom = scenario["per_nom"][index][vegi]
                    scenario_mean = scenario["per_mean"][index][vegi]
                    scenario_sd = scenario["per_sigma"][index][vegi]
                    
                Ks = create_dict_entry(
                            'Ks',scenario,index,vegi, 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            pert_control_name=pert_control_name
                            )
                list_pert.append(Ks)
        
        
        else:
            Ks = create_dict_entry(
                        'Ks',scenario,index,'', 
                        scenario_nom, scenario_mean, scenario_sd, 
                        clip_min,clip_max,
                        NENS,
                        )
            list_pert.append(Ks)
                               
                
    # Check if 'porosity' is a parameter in the scenario and get its index if present
    if "porosity" in scenario["per_name"]:
        index = scenario["per_name"].index("porosity")
    
        # Set perturbation order, defaulting to ['zone', 'root_map', 'layers'] if not specified
        pert_control_name = scenario.get("pert_control_name", ['zone', 'root_map', 'layers'])[index]
    
        # Retrieve nominal, mean, and sigma values for the scenario
        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
    
        # Define clipping bounds for soil porosity
        clip_min, clip_max = 0.2, 0.7
    
        # Perturb by 'zone' if specified
        if 'zone' in pert_control_name:
            for zonei in range(len(np.unique(zone_raster))):
                # Adjust bounds and retrieve parameters for the current zone
                clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                zone_nom, zone_mean, zone_sd = scenario_nom[zonei], scenario_mean[zonei], scenario_sd[zonei]
    
                # Create porosity entry for this zone and validate distribution
                porosity = create_dict_entry(
                    'porosity', scenario, index, zonei, zone_nom, zone_mean, zone_sd, 
                    clip_min, clip_max, NENS
                )
                check_distribution(porosity)
                list_pert.append(porosity)
    
        # Perturb by 'layers' if specified
        elif 'layers' in pert_control_name:
            for nstri in range(nlayers):
                # Adjust bounds and retrieve parameters for the current layer
                clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                layer_nom, layer_mean, layer_sd = scenario_nom[nstri], scenario_mean[nstri], scenario_sd[nstri]
    
                # Create porosity entry for this layer and validate distribution
                porosity = create_dict_entry(
                    'porosity', scenario, index, nstri, layer_nom, layer_mean, layer_sd, 
                    clip_min, clip_max, NENS
                )
                check_distribution(porosity)
                list_pert.append(porosity)
    
        # Default case: apply global perturbation if neither 'zone' nor 'layers' specified
        else:
            clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
            porosity = create_dict_entry(
                'porosity', scenario, index, '', scenario_nom, scenario_mean, scenario_sd, 
                clip_min, clip_max, NENS
            )
            check_distribution(porosity)
            list_pert.append(porosity)
    


    #%% Plant parameters
    # ------------------------------------------------------------------------

    if "PCWLT" in scenario["per_name"]:
        index = scenario["per_name"].index("PCWLT")


        clip_min = -180
        clip_max = None
            
        for nstri in range(len(simu_DA.soil["PCWLT"])):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil["PCWLT"]),
                het_nb=nstri,
            )
            
            PCWLT = {
                "type_parm": "PCWLT" + str(nstri),
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": "$m$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "PCWLT.png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(PCWLT)
            
    if "OMGC" in scenario["per_name"]:
        index = scenario["per_name"].index("OMGC")


        clip_min = 0
        clip_max = 1
            
        for nstri in range(len(simu_DA.soil["OMGC"])):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil["OMGC"]),
                het_nb=nstri,
            )
            
            OMGC = {
                "type_parm": "OMGC" + str(nstri),
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": "$?$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "OMGC.png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(OMGC)
            
            
    if "PZ" in scenario["per_name"]:
        index = scenario["per_name"].index("PZ")


        clip_min = 0
        clip_max = None
            
        for nstri in range(len(simu_DA.soil["PZ"])):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(simu_DA.soil["PZ"]),
                het_nb=nstri,
            )
            
            PZ = {
                "type_parm": "PZ" + str(nstri),
                "nominal": scenario["per_nom"][index],  # nominal value
                "mean": scenario["per_mean"][index],
                "sd": scenario["per_sigma"][index],
                "units": "$?$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "PZ.png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(PZ)
            
    if "PCREF" in scenario["per_name"]:
        index = scenario["per_name"].index("PCREF")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        
        clip_min = -30
        clip_max = -1

        df_SPP, df_FP = simu_DA.read_inputs('soil')

        for nstri in range(len(df_FP.index.unique())):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(df_FP.index.unique()),
                het_nb=nstri,
            )
        
        
            if len(df_FP.index.unique()) > 1:
                scenario_nom = scenario["per_nom"][index][nstri]
                scenario_mean = scenario["per_mean"][index][nstri]
                scenario_sd = scenario["per_sigma"][index][nstri]

            PCREF = {
                "type_parm": "PCREF" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$m$",  # units
                "sampling_type": "normal",
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "PCREF" + str(nstri) + ".png",
                "surf_zones_param": nstri,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(PCREF)

    if "ZROOT" in scenario["per_name"]:
        index = scenario["per_name"].index("ZROOT")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        if "sampling_type" in scenario:
            scenario_sampling = scenario["sampling_type"][index]
        else:
            scenario_sampling = "lognormal"

        clip_min = 0
        clip_max = simu_DA.dem_parameters['base']
        
        df_SPP, df_FP = simu_DA.read_inputs('soil')
        
        for nzonesi in range(len(df_FP.index.unique())):

            clip_min, clip_max = check4bounds(
                scenario,
                index,
                clip_min,
                clip_max,
                het_size=len(df_FP.index.unique()),
                het_nb=nzonesi,
            )

            if len(df_FP.index.unique()) > 1:
                scenario_nom = scenario["per_nom"][index][nzonesi]
                scenario_mean = scenario["per_mean"][index][nzonesi]
                scenario_sd = scenario["per_sigma"][index][nzonesi]

            ZROOT = {
                "type_parm": "ZROOT" + str(nzonesi),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$m$",  # units
                "sampling_type": scenario_sampling,
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "ZROOT" + str(nzonesi) + ".png",
                "surf_zones_param": nzonesi,
                "clip_min": clip_min,
                "clip_max": clip_max,
            }
            list_pert.append(ZROOT)
            # nb_surf_nodes = 110

    return list_pert


#%%


# SAMPLING distribution
# ----------------------

def sampling_dist_trunc(myclip_a, myclip_b, ensemble_size, **kwargs):
    # https://stackoverflow.com/questions/18441779/how-to-specify-upper-and-lower-limits-when-using-numpy-random-normal
    X = truncnorm(
        (myclip_a - kwargs["loc"]) / kwargs["scale"],
        (myclip_b - kwargs["loc"]) / kwargs["scale"],
        loc=kwargs["loc"],
        scale=kwargs["scale"],
    )
    return X.rvs(ensemble_size)


def sampling_dist(sampling_type, mean, sd, ensemble_size, **kwargs):
    # sampling
    np.random.seed(1)
    if sampling_type == "lognormal":
        parm_sampling = np.random.lognormal(np.log(mean), sigma=sd, size=ensemble_size)       
    elif sampling_type == "normal":
        # parm_sampling = np.random.normal(mean, sd, size=ensemble_size)
        parm_sampling = np.random.normal(mean, scale=sd, size=ensemble_size)
    elif sampling_type == "uniform":
        minmax_uni = kwargs["minmax_uni"]
        parm_sampling = np.random.uniform(minmax_uni[0], minmax_uni[1], ensemble_size)
    return parm_sampling


# PERTUBATE distribution
# ----------------------


def perturbate_dist(parm, per_type, parm_sampling, ensemble_size):
    # pertubate
    parm_mat = np.ones(ensemble_size) * parm["nominal"]
    if per_type == None:
        parm_per_array = parm_sampling
    if per_type == "multiplicative":
        parm_per_array = parm_mat * parm_sampling
    elif per_type == "additive":
        parm_per_array = parm_mat + parm_sampling
    return parm_per_array


def Carsel_Parrish_VGN_pert():
    cholesky_diag_mat = np.diag(3)
    pass


def Archie_pert_rules(
    parm, type_parm, ensemble_size, mean, sd, per_type, sampling_type
):
    # a : TYPE, optional
    #     Tortuosity factor. The default is [1.0].
    # m : TYPE, optional
    #     Cementation exponent. The default is [2.0]. (usually in the range 1.3 -- 2.5 for sandstones)
    # n : TYPE, optional
    #     Saturation exponent. The default is [2.0].

    if "rFluid" in type_parm:
        parm_sampling = sampling_dist_trunc(
            myclip_a=0, myclip_b=np.inf, ensemble_size=ensemble_size, loc=mean, scale=sd
        )
    elif "a" in type_parm:
        parm_sampling = sampling_dist_trunc(
            myclip_a=0, myclip_b=2.5, ensemble_size=ensemble_size, loc=mean, scale=sd
        )
    elif "m" in type_parm:
        parm_sampling = sampling_dist_trunc(
            myclip_a=1.3, myclip_b=2.5, ensemble_size=ensemble_size, loc=mean, scale=sd
        )
    elif "n" in type_parm:
        parm_sampling = sampling_dist_trunc(
            myclip_a=2.5, myclip_b=3, ensemble_size=ensemble_size, loc=mean, scale=sd
        )
    else:
        parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size)

    parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    return parm_sampling, parm_per_array


def VG_pert_rules(
    var_per_2add,
    parm,
    type_parm,
    ensemble_size,
    mean,
    sd,
    per_type,
    sampling_type,
    **kwargs,
):
    print(
        "The parameters of the van Genuchten retention curves α,"
        + "n, and θ r are perturbed taking into account their mutual cor-"
        + "relation according to Carsel and Parrish (1988)"
    )

    if "Carsel_Parrish_VGN_pert" in kwargs:
        utils.Carsel_Parrish_1988(soilTexture=None)
        Carsel_Parrish_VGN_pert()
    else:
        if "clip_min" in var_per_2add[type_parm].keys():
            parm_sampling = sampling_dist_trunc(
                myclip_a=var_per_2add[type_parm]["clip_min"],
                myclip_b=var_per_2add[type_parm]["clip_max"],
                ensemble_size=ensemble_size,
                loc=mean,
                scale=sd,
            )
        else:
            parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size)
        parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    return parm_per_array


def atmbc_pert_rules(
        var_per_2add,
    parm,
    type_parm,
    ensemble_size,
    mean,
    sd,
    per_type,
    sampling_type,
    **kwargs,
):

    atmbc_times = parm["data2perturbate"]["time"]
    atmbc_values = parm["data2perturbate"]["VALUE"]
    parm_per_ti = []
    white_noise =  np.random.normal(0, 1, size=(ensemble_size, 
                                                    len(atmbc_times)
                                                    )
                                    )

    Tau = parm["time_decorrelation_len"]

    # Compute perturbed hyetographs
    q = np.zeros_like(white_noise)
    q[:, 0] = white_noise[:, 0]  # Initialize with white noise
    
    deltaT = np.diff(np.r_[atmbc_times,atmbc_times[-1]+86400])
    gamma = np.array([1 - (dt / Tau) for dt in deltaT])
    
    for i in range(1, len(deltaT)):
        q[:, i] = (gamma[i] * q[:, i-1]) + (np.sqrt(1 - gamma[i] ** 2) * white_noise[:, i])

    q_Evenson = q
  
    # Transform using log-normal function
    # m_i, s_i = 1, 0.1  # Example mean and std for transformation
    q_log = np.log(mean) - 0.5 * np.log(((sd / mean) ** 2) + 1) + q_Evenson * np.sqrt(np.log(((sd / mean) ** 2) + 1))
    
    # Back transformation
    q_logback = np.exp(q_log)
    
    # Compute final perturbed forcing
    parm_per_ti = np.multiply(atmbc_values, q_logback)
    
    # Plot perturbed hyetographs
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(ensemble_size):
        ax.step(atmbc_times, parm_per_ti[i, :], where="post", alpha=0.5, label=f"Perturbation {i+1}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Forcing (Perturbed)")
    plt.title("Perturbed Hyetographs")
    # plt.legend()
    # plt.show()

    # var_per_2add[type_parm] = parm_per_ti
    key = "time_variable_perturbation"
    var_per_2add[type_parm][key] = parm_per_ti
    fig.savefig(os.path.join(os.getcwd(), kwargs["savefig"]), dpi=350)


    return white_noise, parm_per_ti


def Johnson1970(self):
    print("not yet implemented - see Botto 2018")


def Evensen2003(q, wk, deltaT, gamma):
    """
    Ensemble Generation of time-variable atmospheric forcing rates.
    Parameters
    ----------
    wk : np.array([])
        is a sequence of white noise drawn from the standard normal distribution.
    deltaT : float
        assimilation interval in sec
    Tau : np.array([])
        the specified time decorrelation length
    time : int
        assimilation time index

    Returns
    -------
    qki : Ensemble of time-variable atmospheric forcing rates at time i

    """
    for i in range(1, len(deltaT)):
        # if Tau < deltaT[i]:
        #     raise ValueError(
        #         "Time decorrelation length is too small; should be at least>=" + str(deltaT)
        #     )
        # else:
            q[:, i] = (gamma[i] * q[:, i-1]) + (np.sqrt(1 - gamma[i] ** 2) * wk[:, i])

    return q


def build_dict_attributes_pert(
    var_per_2add, type_parm, parm_per_array, parm_sampling, **kwargs
):

    key = "ini_perturbation"
    var_per_2add[type_parm][key] = parm_per_array
    key = "sampling"
    var_per_2add[type_parm][key] = parm_sampling
    # Parameter tranformation
    # --------------------------------------------------------------------
    if "transf_type" in kwargs:
        var_per_2add[type_parm]["transf_type"] = kwargs["transf_type"]
        if "transf_bounds" in kwargs:
            var_per_2add[type_parm]["transf_bounds"] = kwargs["transf_bounds"]
    else:
        var_per_2add[type_parm]["transf_type"] = None
    # Parameter spatial extension
    # --------------------------------------------------------------------
    if "surf_zones_param" in kwargs:
        nb_surf_zones = kwargs["surf_zones_param"]
        # parm_per_array = np.tile(parm_per_array,nb_surf_zones)
        var_per_2add[type_parm]["surf_zones_param"] = kwargs["surf_zones_param"]
    return var_per_2add


def perturbate_parm_by_layers(type_parm,
                              ensemble_size,
                              nlayers,
                              mean=[],
                              sd=[],
                              ):
            
        # Define the number of layers and the size of the ensemble
        # num_layers = 5

        # Create a Sobol sampler
        sobol_sampler = qmc.Sobol(d=ensemble_size)
        
        # Generate Sobol sequence samples
        sobol_samples = sobol_sampler.random(n=nlayers)
        
        # Transform the Sobol sequence to follow a normal distribution
        pressure_head_withLayers = np.zeros((ensemble_size, nlayers))
        pressure_head_withLayers = norm.ppf(sobol_samples, 
                                          loc=mean, 
                                          scale=sd
                                          )
        return pressure_head_withLayers

        
def perturbate_parm(
    var_per,
    parm,
    type_parm,
    mean=[],
    sd=[],
    per_type=None,
    sampling_type="lognormal",
    ensemble_size=128,
    seed=True,
    **kwargs,
):
    """
    Perturbate parameter for ensemble generation

    Possible variable to perturbate:
        - initial conditions
        - hyetograph atmbc
        - van Genuchten retention curves parameters
        - Feddes parameters

    Perturbation:
        - parameters with constrainsts (truncated distribution)
        - parameters time dependant
        - other types of parameters

    Returns
    -------
    var_per : dict
        Perturbated variable dictionnary

    """

    savefig = False
    if 'savefig' in kwargs:
       savefig = kwargs["savefig"]
    
    if per_type=='None':
        per_type = None
    

    match_ic_withLayers = re.search(r'ic\d+', type_parm)
    nlayers = None
    if 'nlayers' in kwargs:
        nlayers = kwargs['nlayers']

    var_per_2add = {}

    # copy initiail variable dict and add 'sampling' and 'ini_perturbation' attributes
    # -------------------------------------------------------------------------
    var_per_2add[type_parm] = parm

    key = "sampling_type"
    var_per_2add[type_parm][key] = sampling_type
    key = "sampling_mean"
    var_per_2add[type_parm][key] = mean
    key = "sampling_sd"
    var_per_2add[type_parm][key] = sd
    key = "per_type"
    var_per_2add[type_parm][key] = per_type

    # Contrainsted perturbation (bounded)
    # --------------------------------------------------------------------
    if "Archie" in type_parm:
        parm_sampling, parm_per_array = Archie_pert_rules(
            parm, type_parm, ensemble_size, mean, sd, per_type, sampling_type
        )
    elif "porosity" in type_parm:
        parm_sampling = sampling_dist_trunc(
            myclip_a=0, myclip_b=1, ensemble_size=ensemble_size, loc=mean, scale=sd
        )
        parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)

    elif "zroot".casefold() in type_parm.casefold():
        parm_sampling = sampling_dist_trunc(
            myclip_a=var_per_2add[type_parm]["clip_min"],
            myclip_b=var_per_2add[type_parm]["clip_max"],
            ensemble_size=ensemble_size,
            loc=mean,
            scale=sd,
        )
        parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)

    # check if parameters in part of van Genuchten retention curves
    # ----------------------------------------------------------------------
    # if type_parm in ['Alpha', 'nVG', 'thethaR']: #van Genuchten retention curves
    elif "VG" in type_parm:  # van Genuchten retention curves
        parm_per_array = VG_pert_rules(
            var_per_2add,
            parm,
            type_parm,
            ensemble_size,
            mean,
            sd,
            per_type,
            sampling_type,
            **kwargs,
        )

    # Time dependant perturbation
    # --------------------------------------------------------------------
    elif "atmbc" in type_parm:
    # elif any("atmbc" in item for item in type_parm):
        parm_sampling, parm_per_array = atmbc_pert_rules(    
                                                        var_per_2add,
                                                        parm,
                                                        type_parm,
                                                        ensemble_size,
                                                        mean,
                                                        sd,
                                                        per_type,
                                                        sampling_type,
                                                        **kwargs
                                                        )
        savefig = False

    elif match_ic_withLayers:
        pressure_head_withLayers = perturbate_parm_by_layers(type_parm,
                                                              ensemble_size,
                                                              nlayers,
                                                              mean,
                                                              sd,
                                                              )
        key_root = re.split("(\d+)", type_parm)

        parm_sampling = pressure_head_withLayers[int(key_root[1])]
        parm_per_array = pressure_head_withLayers[int(key_root[1])]
        var_per_2add[type_parm]['ini_ic_withLayers'] = pressure_head_withLayers
    
    
    # For all other types of perturbation
    # --------------------------------------------------------------------
    else:
        if var_per_2add[type_parm]['clip_min'] is not None:
            parm_sampling = sampling_dist_trunc(
                myclip_a=var_per_2add[type_parm]["clip_min"],
                myclip_b=var_per_2add[type_parm]["clip_max"],
                ensemble_size=ensemble_size,
                loc=mean,
                scale=sd,
            )
        else:
            parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size)
        parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)

    # build dictionnary of perturbated variable
    # --------------------------------------------------------------------
    var_per_2add = build_dict_attributes_pert(
        var_per_2add, type_parm, parm_per_array, parm_sampling, **kwargs
    )

    # Add to var perturbated stacked dict
    # ----------------------------------
    var_per = var_per | var_per_2add

    if savefig is not False:
        plt_CT.plot_hist_perturbated_parm(
                                            parm, 
                                            var_per, 
                                            type_parm, 
                                            parm_per_array, 
                                            **kwargs
                                        )

    return var_per



