""" Managing Data Assimilation process pertubation. 
    Take a scenario dictionnary as input describing which and how parameters are perturbated
    Check consistency of the distribution
    Prepare for DA class
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
from pyCATHY import cathy_utils as utils_CT
from pyCATHY.plotters import cathy_plots as plt_CT

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
                         name,scenario,index,nstri, 
                         scenario_nom, scenario_mean, scenario_sd, 
                         clip_min,clip_max,
                         NENS,
                      ):
    dict_parm_entry = {
        "type_parm": name + str(nstri),
        "nominal": scenario_nom,  # nominal value
        "mean": scenario_mean,
        "sd": scenario_sd,
        "units": "",  # units
        "sampling_type": "normal",
        "ensemble_size": NENS,  # size of the ensemble
        "per_type": scenario["per_type"][index],
        "savefig": name + str(nstri) + ".png",
        "surf_zones_param": nstri,
        "clip_min": clip_min,
        "clip_max": clip_max,
    }
    
    return dict_parm_entry


def perturbate(simu_DA, scenario, NENS, pertControl='Layer'):
    """Write a list of dictionaries, each containing all the informations on how to
    perturbate the parameters based on the scenario to consider
    
    
    pertControl = 'Layer'        
    Perturbation per zone is not yet implemented -  Assuming that the 
    dictionnary of perturbated parameters is build per layers i.e. 
    ks0= layer 0, ks1=layer 1, etc...
   
    """

    nzones = len(simu_DA.soil_SPP['SPP_map'].index.get_level_values(0).unique())
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
        
        index = scenario["per_name"].index("ic")

        ic = create_dict_entry(
                    'ic',scenario,index,'', 
                    scenario_nom, scenario_mean, scenario_sd, 
                    clip_min,clip_max,
                    NENS,
                    )
        # check_distribution(ic)
                
        # ic = {
        #     "type_parm": "ic",
        #     "nominal": scenario["per_nom"][index],  # nominal value
        #     "mean": scenario["per_mean"][index],
        #     "sd": scenario["per_sigma"][index],
        #     "units": "pressure head $(m)$",  # units
        #     "sampling_type": "normal",
        #     "ensemble_size": NENS,  # size of the ensemble
        #     "per_type": scenario["per_type"][index],
        #     "savefig": "ic.png",
        # }
        list_pert.append(ic)

    # if "atmbcETp" in scenario["per_name"]:  
    # any("atmbc" in item for item in ['ic', 'atmbcETp'])
    
    
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

    if "ks" in scenario["per_name"]:
        index = scenario["per_name"].index("ks")

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
    
                ks = create_dict_entry(
                            'ks',scenario,index,nstri, 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            )
                
                # ks = {
                #     "type_parm": "ks" + str(nstri),
                #     "nominal": scenario_nom,  # nominal value
                #     "mean": scenario_mean,
                #     "sd": scenario_sd,
                #     "units": "$m.s^{-1}$",  # units
                #     "sampling_type": "lognormal",
                #     "ensemble_size": NENS,  # size of the ensemble
                #     "per_type": scenario["per_type"][index],
                #     "savefig": "ks" + str(nstri) + ".png",
                #     "surf_zones_param": nstri,
                # }
                list_pert.append(ks)
        else:
                ks = create_dict_entry(
                            'ks',scenario,index,'', 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            )
                # check_distribution(ks)
                list_pert.append(ks)
                
                
    if "porosity" in scenario["per_name"]:
        index = scenario["per_name"].index("porosity")

        scenario_nom = scenario["per_nom"][index]
        scenario_mean = scenario["per_mean"][index]
        scenario_sd = scenario["per_sigma"][index]

        clip_min = 0.2  # minimum soil porosity
        clip_max = 0.7  # maximum soil porosity           

        if type(scenario["per_nom"][index]) is list:                   
            for nstri in range(nlayers):
    
                clip_min, clip_max = check4bounds(
                    scenario,
                    index,
                    clip_min,
                    clip_max,
                    het_size=len(simu_DA.soil_SPP["SPP_map"]["POROS"]),
                    het_nb=nstri,
                )
    
                if nlayers > 1:
                    if type(scenario["per_nom"][index]) is list:
                        scenario_nom = scenario["per_nom"][index][nstri]
                        scenario_mean = scenario["per_mean"][index][nstri]
                        scenario_sd = scenario["per_sigma"][index][nstri]
    
                porosity = create_dict_entry(
                            'porosity',scenario,index,nstri, 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            )
                check_distribution(porosity)
                list_pert.append(porosity)
        else:
                clip_min, clip_max = check4bounds(
                    scenario,
                    index,
                    clip_min,
                    clip_max,
                )

                scenario_nom = scenario["per_nom"][index]
                scenario_mean = scenario["per_mean"][index]
                scenario_sd = scenario["per_sigma"][index]
                
                porosity = create_dict_entry(
                            'porosity',scenario,index,'', 
                            scenario_nom, scenario_mean, scenario_sd, 
                            clip_min,clip_max,
                            NENS,
                            )
                
                # porosity = {
                #     "type_parm": "porosity" + str(''),
                #     "nominal": scenario_nom,  # nominal value
                #     "mean": scenario_mean,
                #     "sd": scenario_sd,
                #     "units": "",  # units
                #     "sampling_type": "normal",
                #     "ensemble_size": NENS,  # size of the ensemble
                #     "per_type": scenario["per_type"][index],
                #     "savefig": "porosity" + str('') + ".png",
                #     "surf_zones_param": '',
                #     "clip_min": clip_min,
                #     "clip_max": clip_max,
                # }
                
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
        clip_max = None
        
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

            ZROOT = {
                "type_parm": "ZROOT" + str(nstri),
                "nominal": scenario_nom,  # nominal value
                "mean": scenario_mean,
                "sd": scenario_sd,
                "units": "$m$",  # units
                "sampling_type": scenario_sampling,
                "ensemble_size": NENS,  # size of the ensemble
                "per_type": scenario["per_type"][index],
                "savefig": "ZROOT" + str(nstri) + ".png",
                "surf_zones_param": nstri,
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
    X = stats.truncnorm(
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
        parm_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
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

    parm_sampling = sampling_dist(sampling_type, mean, sd, ensemble_size)
    # parm_sampling = sampling_dist(sampling_type, mean, 1e-6, ensemble_size)
    parm_per_array = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
    var_per_2add[type_parm] = parm
    Tau = parm["time_decorrelation_len"]
    qk_0 = parm_per_array
    atmbc_times = parm["data2perturbate"]["time"]
    atmbc_values = parm["data2perturbate"]["VALUE"]
    # print(per_type)
    # print(atmbc_values)

    # # Create a histogram plot for perturbed sand conductivities
    # plt.figure(figsize=(10, 5))
    # # plt.hist(parm_sampling, bins=20, alpha=0.5, color='blue', label='Atmbc ENS')
    # plt.hist(parm_per_array, bins=20, alpha=0.5, color='red', label='* Atmbc ENS')
    # plt.xlabel('Perturbed atmbc')
    # plt.ylabel('Frequency')
    # # plt.title('Perturbed Saturated Hydraulic Conductivities')
    # plt.legend(loc='upper right')
    # plt.grid(True)
    
    
    parm_per_array_time_variable = []
    for i, t in enumerate(atmbc_times):
        if i == 0:
            # qk_0 = wk0
            parm_per_array_time_variable.append(qk_0)
        else:
            qk_0 = parm_per_array_time_variable[i - 1]
            parm["nominal"] = atmbc_values[i]
            # wk = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
            # wk = perturbate_dist(parm, per_type, parm_sampling, ensemble_size)
            wk = sampling_dist(sampling_type, mean, sd, ensemble_size)
            # plt.plot(wk)
            deltaT = abs(atmbc_times[i] - atmbc_times[i - 1])
            # qk_i = Evensen2003(qk_0, wk, deltaT, Tau)
            qk_i = wk
            # plt.plot(qk_i)
            parm_per_array_time_variable.append(qk_i)

    # print(parm["nominal"])
    # print(parm_per_array_time_variable)
    key = "time_variable_perturbation"
    var_per_2add[type_parm][key] = parm_per_array_time_variable
    
    return parm_sampling, parm_per_array, var_per_2add


def Johnson1970(self):
    print("not yet implemented - see Botto 2018")


def Evensen2003(qk_0, wk, deltaT, Tau):
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
    if Tau < deltaT:
        raise ValueError(
            "Time decorrelation length is too small; should be at least>=" + str(deltaT)
        )
    gamma = 1 - deltaT / Tau
    qki = gamma * qk_0 + np.sqrt(1 - (gamma * gamma)) * wk

    return qki


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
    if kwargs["savefig"]:
       savefig = kwargs["savefig"]
    
    if per_type=='None':
        per_type = None
        
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
        parm_sampling, parm_per_array, var_per_2add = atmbc_pert_rules(    
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



