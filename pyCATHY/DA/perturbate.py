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
                         per_type,
                         pert_control_name=None,
                         transf_type=None
                      ):
    
    # per_type = scenario["per_type"][index]
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
        "transf_type": transf_type,
        "savefig": name + str(nstri) + ".png",
        "surf_zones_param": nstri,
        "clip_min": clip_min,
        "clip_max": clip_max,
        "pert_control_name": pert_control_name,
    }
    
    return dict_parm_entry

def process_perturbation_block(
    scenario, 
    NENS,
    index,
    nlayers=None, 
    zone_raster=None,
    veg_raster=None,
):
    """
    Process perturbation for a single parameter (given by index) from scenario dict.

    Returns a list of perturbation dict entries (does NOT append to external list).

    Parameters:
    - scenario: dict with keys 'per_name', 'per_nom', 'per_mean', 'per_sigma', 'per_type', optional 'per_bounds' and 'pert_control_name'.
    - NENS: int, number of ensemble members.
    - index: int, index of the parameter in scenario dict lists to process.
    - nlayers: int, required for 'layers' mode.
    - zone_raster: np.ndarray, required for 'zone' mode.
    - veg_raster: np.ndarray, required for 'root_map' mode.
    """
    entries = []
    pname = scenario["per_name"][index]
    pert_control_name = scenario.get('pert_control_name', [None]*len(scenario["per_name"]))[index]
    scenario_nom = scenario["per_nom"][index]
    scenario_mean = scenario["per_mean"][index]
    scenario_sd = scenario["per_sigma"][index]
    per_type = scenario["per_type"][index]
    transf_type = scenario["transf_type"][index]

    # Pull bounds dynamically if available
    clip_min, clip_max = None, None
    if "per_bounds" in scenario:
        bounds = scenario["per_bounds"][index]
        if isinstance(bounds, (list, tuple)) and len(bounds) == 2:
            clip_min, clip_max = bounds

    # Handle perturbation by control type
    if isinstance(scenario_nom, list):
        if pert_control_name == 'layers':
            if nlayers is None:
                raise ValueError("nlayers must be provided for 'layers' perturbation control.")
            for layer_idx in range(nlayers):
                entry = create_dict_entry(
                    pname, scenario, index, layer_idx,
                    scenario_nom[layer_idx],
                    scenario_mean[layer_idx],
                    scenario_sd[layer_idx],
                    clip_min, clip_max, NENS,
                    per_type=per_type[layer_idx],
                    pert_control_name=pert_control_name,
                    transf_type = transf_type
                )
                check_distribution(entry)
                entries.append(entry)

        elif pert_control_name == 'zone':
            if zone_raster is None:
                raise ValueError("zone_raster must be provided for 'zone' perturbation control.")
            unique_zones = np.unique(zone_raster)
            for zone_idx in range(len(unique_zones)):
                z_nom = scenario_nom[zone_idx]
                z_mean = scenario_mean[zone_idx]
                z_sd = scenario_sd[zone_idx]
                z_type = per_type[zone_idx]

                entry = create_dict_entry(
                    pname, scenario, index, zone_idx,
                    z_nom, z_mean, z_sd,
                    clip_min, clip_max, NENS,
                    per_type=z_type,
                    pert_control_name=pert_control_name,
                    transf_type = transf_type
                )
                check_distribution(entry)
                entries.append(entry)

        elif pert_control_name == 'root_map':
            if veg_raster is None:
                raise ValueError("veg_raster must be provided for 'root_map' perturbation control.")
            unique_roots = np.unique(veg_raster)
            for veg_idx in range(len(unique_roots)):
                r_nom = scenario_nom[veg_idx]
                r_mean = scenario_mean[veg_idx]
                r_sd = scenario_sd[veg_idx]
                r_type = per_type[veg_idx] if isinstance(per_type, list) else per_type

                entry = create_dict_entry(
                    pname, scenario, index, veg_idx,
                    r_nom, r_mean, r_sd,
                    clip_min, clip_max, NENS,
                    per_type=r_type,
                    pert_control_name=pert_control_name,
                    transf_type = transf_type
                )
                check_distribution(entry)
                entries.append(entry)

        else:
            raise ValueError(
                f"Unknown list-type perturbation control '{pert_control_name}' for parameter '{pname}'. "
                "Expected 'layers', 'zone', or 'root_map'."
            )

    else:
        # Scalar perturbation
        entry = create_dict_entry(
            pname, scenario, index, '',
            scenario_nom, scenario_mean, scenario_sd,
            clip_min, clip_max, NENS,
            per_type=per_type,
            pert_control_name=pert_control_name,
            transf_type = transf_type
        )
        check_distribution(entry)
        entries.append(entry)

    return entries



def perturbate(simu_DA, scenario, NENS, 
               **kwargs):
    """Write a list of dictionaries, each containing all the informations on how to
    perturbate the parameters based on the scenario to consider
    
    
    pertControl = 'Layer'        
    Perturbation per zone is not yet implemented -  Assuming that the 
    dictionnary of perturbated parameters is build per layers i.e. 
    ks0= layer 0, ks1=layer 1, etc...
   
    """  
    
    df_SPP, df_FP = simu_DA.read_inputs('soil')
    zone_raster, zone_header = simu_DA.read_inputs('zone')
    veg_raster, veg_header = simu_DA.read_inputs('root_map')
        
    if 'nzones' in kwargs: 
        nzones = kwargs.pop('nzones')
    else:
        nzones = len(df_SPP.index.get_level_values(0).unique())
    if 'nlayers' in kwargs: 
        nzones = kwargs.pop('nlayers')
    else:
        nlayers = len(df_SPP.index.get_level_values(1).unique())
        nlayers_PERMX = len(df_SPP['PERMX'].unique())

    #%%
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

        pert_control_name = None
        if 'pert_control_name' in scenario:
            pert_control_name = scenario['pert_control_name'][index]
        
        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        
        
        clip_min = None
        clip_max = None
        if type(scenario["per_nom"][index]) is list:    
            if pert_control_name=='layers':
                for nstri in range(nlayers):
                    
                    if nlayers > 1:
                        scenario_nom = scenario["per_nom"][index][nstri]
                        scenario_mean = scenario["per_mean"][index][nstri]
                        scenario_sd = scenario["per_sigma"][index][nstri]
                        scenario_per_type = scenario["per_type"][index][nstri]
                        transf_type = scenario["transf_type"][index][nstri]
        
                    ic = create_dict_entry(
                                'ic',scenario,index,nstri, 
                                scenario_nom, scenario_mean, scenario_sd, 
                                clip_min,clip_max,
                                NENS,
                                per_type=scenario_per_type,
                                transf_type = transf_type
                                )
                    
                    list_pert.append(ic)
                    
            elif pert_control_name=='zone':
                for zonei in range(len(np.unique(zone_raster))):
                    # Adjust bounds and retrieve parameters for the current zone
                    clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                    zone_nom, zone_mean, zone_sd = scenario_nom[zonei], scenario_mean[zonei], scenario_sd[zonei]
                    scenario_per_type = scenario["per_type"][index][zonei]

                    ic = create_dict_entry(
                        'ic', scenario, index, zonei, zone_nom, zone_mean, zone_sd, 
                        clip_min, clip_max, NENS,
                        per_type=scenario_per_type,
                        pert_control_name=pert_control_name,
                        transf_type = scenario["transf_type"][index]
                    )
                    check_distribution(ic)
                    list_pert.append(ic)        
            else:
                raise ValueError(f"List-type perturbation provided, but unknown control name: '{pert_control_name}'. "
                         f"Expected 'layers' or 'zone'.")
        else: 
            
            ic = create_dict_entry(
                        'ic',scenario,index,'', 
                        scenario_nom, scenario_mean, scenario_sd, 
                        clip_min,clip_max,
                        NENS,
                        per_type=scenario["per_type"][index],
                        transf_type = scenario["transf_type"][index]
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
    
    if "SATCELL_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("SATCELL_VG")
        SATCELL_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        list_pert.extend(SATCELL_VG)


    if "thetar_VG" in scenario["per_name"]:
        
        index = scenario["per_name"].index("thetar_VG")
        thetar_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        list_pert.extend(thetar_VG)
        
    if "alpha_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("alpha_VG")          
        alpha_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        list_pert.extend(alpha_VG)

    if "n_VG" in scenario["per_name"]:
        index = scenario["per_name"].index("n_VG")
        n_VG = process_perturbation_block(scenario,
                                                NENS=NENS, 
                                                index=index,
                                                nlayers=nlayers, 
                                                zone_raster=zone_raster, 
                                                veg_raster=veg_raster
                                                )
        
        list_pert.extend(n_VG)


    # if "VGP" in scenario["per_name"]:
    #     # - 'PERMX' (NSTR, nstriONE): saturated hydraulic conductivity - xx
    #     # - 'PERMY' (NSTR, nstriONE): saturated hydraulic conductivity - yy
    #     # - 'PERMZ' (NSTR, nstriONE): saturated hydraulic conductivity - zz
    #     # - 'ELSTOR' (NSTR, nstriONE): specific storage
    #     # - 'POROS'  (NSTR, nstriONE): porosity (moisture content at saturation) = \thetaS

    #     # retention curves parameters VGN, VGRMC, and VGPSAT
    #     # - 'VGNCELL' (NSTR, nstriONE): van Genuchten curve exponent  = n
    #     # - 'VGRMCCELL' (NSTR, nstriONE): residual moisture content = \thetaR
    #     # - 'VGPSATCELL' (NSTR, nstriONE): van Genuchten curve exponent -->
    #     #                               VGPSAT == -1/alpha (with alpha expressed in [L-1]);
    #     # ['ks','ss','phi','thetar','alpha','n'])
    #     # ['PERMX','ELSTOR','POROS','VGRMCCELL','VGPSATCELL','VGNCELL'])

    #     # need to account for transformation of the parameters see Boto

    #     VGP_units = ["", "", "", "", "", ""]

    #     for i, p in enumerate(["ks", "ss", "phi", "thetar", "alpha", "r"]):
    #         index = scenario["per_name"].index("VGP")

    #         if p not in scenario["per_nom"][index]:
    #             pass
    #         else:
    #             p_VGP = {
    #                 "type_parm": p,
    #                 "nominal": scenario["per_nom"][index][p],  # nominal value
    #                 "mean": scenario["per_mean"][index][p],
    #                 "sd": scenario["per_sigma"][index][p],
    #                 "units": VGP_units[i],  # units
    #                 "sampling_type": "normal",
    #                 "ensemble_size": NENS,  # size of the ensemble
    #                 "per_type": scenario["per_type"][index],
    #                 "savefig": "_VGP" + p + ".png",
    #             }
    #             # list_pert.append(p_VGP)

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

        # archie_units = ["", "", "", ""]
        archie_units = {
            "porosity": "-",                   # dimensionless (volume fraction)
            "a": "-",                          # tortuosity factor (dimensionless)
            "m": "-",                          # cementation exponent (dimensionless)
            "n": "-",                          # saturation exponent (dimensionless)
            "rFluid_Archie": "ohm.m",           # resistivity of water
        }


        for i, p in enumerate(Archie_2pert):
            
            index = scenario["per_name"].index(p)      
            pert_control_name = None
            if 'pert_control_name' in scenario:
                pert_control_name = scenario['pert_control_name'][index]
                
            
            scenario_nom = scenario["per_nom"][index]
            scenario_mean = scenario["per_mean"][index]
            scenario_sd = scenario["per_sigma"][index]
            scenario_per_type = scenario["per_type"][index]



            if pert_control_name is not None and pert_control_name == 'layers':
                for nstri in range(nlayers):
                    # Adjust bounds and retrieve parameters for the current layer
                    # clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                    layer_nom, layer_mean, layer_sd = scenario_nom[nstri], scenario_mean[nstri], scenario_sd[nstri]
                    # Create porosity entry for this layer and validate distribution
                    p_archie = create_dict_entry(
                        p, scenario, index, nstri, layer_nom, layer_mean, layer_sd, 
                        clip_min, clip_max, NENS,
                        per_type=scenario_per_type[nstri],
                    )
                    check_distribution(p_archie)
                    list_pert.append(p_archie)
            else:
                # pass
                index = scenario["per_name"].index(p)
        
                p_archie = {
                    "type_parm": p,
                    "nominal": scenario["per_nom"][index],  # nominal value
                    "mean": scenario["per_mean"][index],
                    "sd": scenario["per_sigma"][index],
                    "units": archie_units[p],  # units
                    "sampling_type": "normal",
                    "ensemble_size": NENS,  # size of the ensemble
                    "per_type": scenario["per_type"][index],
                    "savefig": "_Archie" + p + ".png",
                    'pert_control_name':pert_control_name

                }
                list_pert.append(p_archie)

    #%% SOIL parameters
    # ------------------------------------------------------------------------


    if "Ks" in scenario["per_name"]:
        index = scenario["per_name"].index("Ks")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        scenario_per_type = scenario["per_type"][index]

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
           
                zone_nom, zone_mean, zone_sd = scenario_nom[zonei], scenario_mean[zonei], scenario_sd[zonei]
                Ks = create_dict_entry(
                            'Ks',scenario,index,zonei, 
                            zone_nom, zone_mean, zone_sd, 
                            clip_min,clip_max,
                            NENS,
                            per_type=scenario_per_type[zonei],
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
                            pert_control_name=pert_control_name,
                            per_type=scenario_per_type[vegi],
                            )
                list_pert.append(Ks)
                
        # Perturb by 'layers' if specified
        elif pert_control_name == 'layers':
            for nstri in range(nlayers):
                # Adjust bounds and retrieve parameters for the current layer
                clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                layer_nom, layer_mean, layer_sd = scenario_nom[nstri], scenario_mean[nstri], scenario_sd[nstri]
    
                # Create porosity entry for this layer and validate distribution
                Ks = create_dict_entry(
                    'Ks', scenario, index, nstri, layer_nom, layer_mean, layer_sd, 
                    clip_min, clip_max, NENS,
                    per_type=scenario_per_type[nstri],
                    pert_control_name=pert_control_name
                )
                check_distribution(Ks)
                list_pert.append(Ks)
        
        else:
            Ks = create_dict_entry(
                        'Ks',scenario,index,'', 
                        scenario_nom, scenario_mean, scenario_sd, 
                        clip_min,clip_max,
                        NENS,
                        per_type=scenario["per_type"][index],
                        )
            list_pert.append(Ks)
                               
                
    # Check if 'porosity' is a parameter in the scenario and get its index if present
    if "porosity" in scenario["per_name"]:
        index = scenario["per_name"].index("porosity")
    
        # Set perturbation order, defaulting to ['zone', 'root_map', 'layers'] if not specified
        
        
        pert_control_name = None
        if 'pert_control_name' in scenario:
            pert_control_name = scenario['pert_control_name'][index]
        
        
        # Retrieve nominal, mean, and sigma values for the scenario
        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]
        scenario_per_type = scenario["per_type"][index]
    
        # Define clipping bounds for soil porosity
        clip_min, clip_max = 0.2, 0.7
    
        # Perturb by 'zone' if specified
        if pert_control_name=='zone':
            for zonei in range(len(np.unique(zone_raster))):
                # Adjust bounds and retrieve parameters for the current zone
                clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                zone_nom, zone_mean, zone_sd = scenario_nom[zonei], scenario_mean[zonei], scenario_sd[zonei]
    
                # Create porosity entry for this zone and validate distribution
                porosity = create_dict_entry(
                    'porosity', scenario, index, zonei, zone_nom, zone_mean, zone_sd, 
                    clip_min, clip_max, NENS,
                    per_type=scenario_per_type[zonei],
                    pert_control_name=pert_control_name
                )
                check_distribution(porosity)
                list_pert.append(porosity)
    
        # Perturb by 'layers' if specified
        elif pert_control_name=='layers':
            for nstri in range(nlayers):
                # Adjust bounds and retrieve parameters for the current layer
                clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
                layer_nom, layer_mean, layer_sd = scenario_nom[nstri], scenario_mean[nstri], scenario_sd[nstri]
    
                # Create porosity entry for this layer and validate distribution
                porosity = create_dict_entry(
                    'porosity', scenario, index, nstri, layer_nom, layer_mean, layer_sd, 
                    clip_min, clip_max, NENS,
                    per_type=scenario_per_type[nstri],
                    pert_control_name=pert_control_name
                )
                check_distribution(porosity)
                list_pert.append(porosity)
    
        # Default case: apply global perturbation if neither 'zone' nor 'layers' specified
        else:
            clip_min, clip_max = check4bounds(scenario, index, clip_min, clip_max)
            porosity = create_dict_entry(
                'porosity', scenario, index, '', scenario_nom, scenario_mean, scenario_sd, 
                clip_min, clip_max, NENS,
                per_type=scenario["per_type"][index],
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

    # if "Carsel_Parrish_VGN_pert" in kwargs:
    #     utils.Carsel_Parrish_1988(soilTexture=None)
    #     Carsel_Parrish_VGN_pert()
    # else:
    #     if "clip_min" in var_per_2add[type_parm].keys():
    #         parm_sampling = sampling_dist_trunc(
    #             myclip_a=var_per_2add[type_parm]["clip_min"],
    #             myclip_b=var_per_2add[type_parm]["clip_max"],
    #             ensemble_size=ensemble_size,
    #             loc=mean,
    #             scale=sd,
    #         )
    #     else:
    #         parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size)
    #     parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    
    """
    Sample a single van Genuchten physical parameter and transform it to Gaussian space via Johnson transformation.

    Parameters:
    - parm: dict
        Must include:
            - 'nominal': float
            - 'sampling_type': 'lognormal' | 'normal' | 'uniform'
            - 'mean': mean of sampling distribution (in physical space)
            - 'sd': standard deviation of sampling
            - 'per_type': 'additive' | 'multiplicative' | None
            - 'type': 'LN' | 'SB' | 'SU' (Johnson transform type)
            - 'A': Johnson transform A parameter (for SB)
            - 'B': Johnson transform B parameter (for SB)
            - Optional: 'minmax_uni': [min, max] (if uniform)

    - ensemble_size: int
    - seed: int

    Returns:
    - v_array: (ensemble_size,) ndarray of physical parameter samples
    - y_array: (ensemble_size,) ndarray of transformed (Gaussian) samples
    """
    
    print(
        "(NOT YET IMPLEMENTED - The parameters of the van Genuchten retention curves α,"
        + "n, and θ r are perturbed taking into account their mutual cor-"
        + "relation according to Carsel and Parrish (1988)"
    )
    
    np.random.seed(seed=1)

    # --- Sampling ---
    sampling_type = parm["sampling_type"]
    mean = parm["mean"]
    sd = parm["sd"]
    per_type = parm.get("per_type", None)
    transf_type = parm.get("transf_type", None)

    kwargs = {}
    if sampling_type == "uniform" and "minmax_uni" in parm:
        kwargs["minmax_uni"] = parm["minmax_uni"]

    # parm_sampling = sampling_dist(
    #     sampling_type=sampling_type,
    #     mean=mean,
    #     sd=sd,
    #     ensemble_size=ensemble_size,
    #     **kwargs
    # )

    parm_sampling = sampling_dist_trunc(
        myclip_a=parm["clip_min"],
        myclip_b=parm["clip_max"],
        ensemble_size=ensemble_size,
        loc= parm['mean'],
        scale=sd,
    )
            
    # --- Perturbation (physical space) ---
    v_array = perturbate_dist(
        parm=parm,
        per_type=per_type,
        parm_sampling=parm_sampling,
        ensemble_size=ensemble_size
    )

    # # --- Johnson transform (to Gaussian space) ---
    # y_array = Johnson1970_transform(
    #     v_array,
    #     transform_type=transf_type,
    #     A=parm['clip_min'],
    #     B=parm['clip_max']
    # )

    return parm_sampling, v_array

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


# def Johnson1970():
#     print("not yet implemented - see Botto 2018")


def Johnson1970_transform(V, transform_type, A=None, B=None):
    """
    Apply Johnson transformation to a parameter V. ( see Botto 2018)

    Parameters:
    - V: float or ndarray
        Original parameter(s) to transform.
    - transform_type: str
        One of 'LN', 'SB', or 'SU' indicating which Johnson transformation to use.
    - A, B: float
        Bounds required for the 'SB' transformation.

    Returns:
    - Y: float or ndarray
        Transformed variable (normally distributed).
    """
    rootName = 'Johnson_'
    if transform_type == rootName+'LN':
        return np.log(V)
    elif transform_type == rootName+'SB':
        if A is None or B is None:
            raise ValueError("A and B must be provided for the SB transformation.")
        U = (V - A) / (B - V)
        return np.log(U)
    elif transform_type == rootName+'SU':
        return np.arcsinh(V)
    else:
        raise ValueError("Invalid transform_type. Choose from 'LN', 'SB', or 'SU'.")

def inverse_Johnson1970_transform(Y, transform_type, A=None, B=None):
    """
    Inverse Johnson transformation from Y to V.

    Parameters:
    - Y: float or ndarray
        Transformed (normal) variable.
    - transform_type: str
        One of 'LN', 'SB', or 'SU'.
    - A, B: float
        Bounds required for 'SB' transformation.

    Returns:
    - V: float or ndarray
        Original variable value.
    """
    if transform_type == 'LN':
        return np.exp(Y)
    elif transform_type == 'SB':
        if A is None or B is None:
            raise ValueError("A and B must be provided for the SB transformation.")
        U = np.exp(Y)
        return (A + B * U) / (1 + U)
    elif transform_type == 'SU':
        return np.sinh(Y)
    else:
        raise ValueError("Invalid transform_type. Choose from 'LN', 'SB', or 'SU'.")


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
    
    # type_parm2check = 
    # match_withLayers = re.search(rf'{type_parm}\d+', type_parm)
    # match_withLayers = re.search(r'\d', type_parm)
    nlayers = None
    if 'nlayers' in kwargs:
        nlayers = kwargs['nlayers']

    var_per_2add = {}
    
    # pert_control_name

    # copy initiail variable dict and add 'sampling' and 'ini_perturbation' attributes
    # -------------------------------------------------------------------------
    var_per_2add[type_parm] = parm

    key = "sampling_type"
    var_per_2add[type_parm][key] = sampling_type
    key = "sampling_mean"
    var_per_2add[type_parm][key] = mean
    key = "sampling_nominal"
    var_per_2add[type_parm][key] = parm['nominal']
    key = "sampling_sd"
    var_per_2add[type_parm][key] = sd
    key = "per_type"
    var_per_2add[type_parm][key] = per_type

    key_root = re.split("(\d+)", type_parm)
    parm_key_nb = 0
    if len(key_root)>1:
        parm_key_nb = int(key_root[1])
       
        
    # Contrainsted perturbation (bounded)
    # --------------------------------------------------------------------
    if "Archie" in type_parm:
        if parm['pert_control_name']=='layers':
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                        ensemble_size,
                                                        nlayers,
                                                        mean,
                                                        sd,
                                                        )
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
        else:
            parm_sampling, parm_per_array = Archie_pert_rules(
                parm, type_parm, ensemble_size, mean, sd, per_type, sampling_type
            )
    elif "porosity" in type_parm:
        if parm['pert_control_name']=='layers':
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                                  ensemble_size,
                                                                  nlayers,
                                                                  mean,
                                                                  sd,
                                                                  )
            key_root = re.split("(\d+)", type_parm)
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
    
        else:
            parm_sampling = sampling_dist_trunc(
                myclip_a=0, myclip_b=1, ensemble_size=ensemble_size, loc=mean, scale=sd
            )
            # parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
            
            if len(key_root)>1:
                # Case with different zones
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )
            else:
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )  
    
    
    
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
        parm_sampling, parm_per_array = VG_pert_rules(
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

    elif "ic" in type_parm:

        if parm['pert_control_name']=='layers':
            param_withLayers = perturbate_parm_by_layers(type_parm,
                                                                  ensemble_size,
                                                                  nlayers,
                                                                  mean,
                                                                  sd,
                                                                  )
            key_root = re.split("(\d+)", type_parm)
    
            parm_sampling = param_withLayers[parm_key_nb]
            parm_per_array = param_withLayers[parm_key_nb]
            var_per_2add[type_parm][f'ini_{type_parm}_withLayers'] = param_withLayers
        else:
            parm_sampling = sampling_dist(sampling_type,  parm['mean'], sd, ensemble_size)
            
            if len(key_root)>1:
                # Case with different zones
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )
            else:
                parm_per_array = perturbate_dist(parm,
                                                 per_type,
                                                 parm_sampling,
                                                 ensemble_size
                                                 )  
    
    # For all other types of perturbation
    # --------------------------------------------------------------------
    else:
        if var_per_2add[type_parm]['clip_min'] is not None:
            parm_sampling = sampling_dist_trunc(
                myclip_a=var_per_2add[type_parm]["clip_min"],
                myclip_b=var_per_2add[type_parm]["clip_max"],
                ensemble_size=ensemble_size,
                loc= parm['mean'],
                scale=sd,
            )
        else:
            parm_sampling = sampling_dist(sampling_type,  parm['mean'], sd, ensemble_size)
        parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)

    # build dictionnary of perturbated variable
    # --------------------------------------------------------------------
    var_per_2add = build_dict_attributes_pert(
                                                var_per_2add, 
                                                type_parm, 
                                                parm_per_array, 
                                                parm_sampling, 
                                                **kwargs
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



